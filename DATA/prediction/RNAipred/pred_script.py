"""
Sign prediction based on methods from: https://www.nature.com/nmeth/journal/v11/n1/full/nmeth.2733.html
Uses only z-scores
    is positive regulation if it's above +1
    is negative regulation below -1
    is 0 in-between +1 and -1
In a given RNAi screen, if both interacting proteins have nonzero values, then the relationship is classified as either a
positive correlation (both  +1  or  both  −1)  or  a  negative  correlation  (one  is  +1  and another is −1).
For each interacting pair, we computed the total number  of  positive  and  negative  correlations.
Creates gene-phenotype matrix based on data from: http://www.genomernai.org/DownloadAllExperimentsForm
    needed header data: PMID, organism, score type

"""

# Imports
import sqlite3, io, math
import itertools as it
import pandas as pd
import numpy as np
import logging

# Defining constants
DATA_FILE = ['../prediction/RNAipred/files/GenomeRNAi_v16_Homo_sapiens.txt']
BUILDER_SET = 'ARN2_layers.db'
MAPPER_DB = 'mapper.db'
SPECIES_DICT = {
    "Homo sapiens": "taxid:9606"}

# # Initiating logger
# logger = logging.getLogger()
# handler = logging.FileHandler('../../workflow/SLK3.log')
# logger.setLevel(logging.DEBUG)
# handler.setLevel(logging.DEBUG)
# formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
# handler.setFormatter(formatter)
# logger.addHandler(handler)

# Setting up sqlite3 connection
# Read database to tempfile
conn2 = sqlite3.connect(MAPPER_DB)
tempfile = io.StringIO()
for line in conn2.iterdump():
    tempfile.write('%s\n' % line)
conn2.close()
tempfile.seek(0)

# Create a database in memory and import from tempfile
DB = sqlite3.connect(":memory:")
with DB:
    DB.cursor().executescript(tempfile.read())
    DB.cursor().execute("CREATE INDEX mapp_uniprot ON MAPP(uniprot_ac);")
    DB.cursor().execute("CREATE INDEX uniprot_id ON UNIPROT_AC(id);")
    DB.cursor().execute("CREATE INDEX mapp_foreign ON MAPP(foreign_id);")
    DB.cursor().execute("CREATE INDEX uniprot_uniprot ON UNIPROT_AC(uniprot_ac);")


def main(logger):

    header_dict = {}
    data = []
    pheno_dict = {}
    genes = []
    mapped_genes_dict = {}
    phenotypes = []
    score_dict = {}
    score_list = []

    df = pd.DataFrame()
    # Counters
    good_score_count = 0
    null_score_count = 0
    cutoff_count = 0
    notcut_count = 0
    entry_counter = 0
    # If not Z or P score
    not_searched_score_count = 0
    # If valueerror
    error_count = 0
    # If 'np' or 'sp' instead of score
    no_score_count = 0
    with DB:
        c = DB.cursor()
        for org_file in DATA_FILE:
            with open(org_file, encoding="ISO-8859-1") as infile:
                # Splitting the files into entries
                for key, group in it.groupby(infile, lambda line: line.startswith('//')):
                    # If the line is the entry separator, empties lists and dictionaries
                    if key is True:
                        header_dict = {}
                        data = []
                        pheno_dict = {}
                        genes = []
                        mapped_genes_dict = {}
                        phenotypes = []
                        score_dict = {}
                        score_list = []
                        entry_counter += 1
                        # Limit
                        # if entry_counter == 26:
                        #     break
                    # Iterating the entry
                    else:
                        for line in group:
                            line = line.strip()
                            # Getting data from header
                            if line[0] == '#':
                                if len(line.split('=')) <= 1:
                                    pass
                                else:
                                    header_dict[line.split("=")[0]] = line.split('=')[1]
                            # Getting data
                            else:
                                if header_dict['#Score Type'] == 'Z-score' or header_dict['#Score Type'] == 'P-value' \
                                        or header_dict['#Score Type'] == 'Fold inhibition (p-value)':
                                    data = line.split('\t')
                                    if data[0] != '//' and data[5] != 'np' and data[5] != 'sp':
                                        # For entries with multiple scores, handle them as different entries
                                        if '/' in data[5]:  # TODO split
                                            score = data[5].split('/')[0]
                                        elif data[5] == '>= 2' or data[5] == '> 2' or data[5] == '1.5 - 2' or data[5] == '>2'\
                                                or data[5] == '> 4' or data[5] == '2 - 4':
                                            score = '2'
                                        elif data[5] == '<= -2' or data[5] == '< -2' or data[5] == '(-2) - (-1.5)' \
                                                or data[5] == '-2 - -1.2' or data[5] == '-2 - -1.5' or data[5] == '<-2'\
                                                or data[5] == '-4 - -2' or data[5] == '< -4':
                                            score = '-2'
                                        elif data[5] == '-1.2 - -0.75':
                                            score = '-1.2'
                                        else:
                                            score = data[5]

                                        if data[3] != '':
                                            gene = data[3]
                                        else:
                                            gene = data[2]
                                            if "'" in gene:
                                                gene = gene.replace("'", '')
                                        phenotype = data[6]

                                        # List of phenotypes
                                        if phenotype not in phenotypes:
                                            phenotypes.append(phenotype)

                                        #print(gene)
                                        # List of genes
                                        if gene not in genes:
                                            genes.append(gene)
                                        # Gene-phenotype dictionary
                                        if gene not in pheno_dict:
                                            pheno_dict[gene] = phenotype
                                        else:
                                            pheno_dict[gene] = phenotype

                                        # Gene-score dictionary
                                        if header_dict['#Score Type'] == 'Fold inhibition (p-value)':
                                            score_dict[gene] = float(score.split(' ')[1].replace('(', '',).replace(')', ''))
                                            good_score_count += 1

                                        else:
                                            try:
                                                score_dict[gene] = float(score)
                                            except ValueError:
                                                print(score)
                                                error_count += 1
                                                continue
                                    # Counting instances where score is np or sp
                                    else:
                                        no_score_count += 1
                                # Counting instances where the score type wasn't Z or P
                                else:
                                    not_searched_score_count += 1

                        logging.debug('entry:' + str(entry_counter))

                        # Declaring additional entry variables (from header)
                        PMID = 'PMID:' + header_dict['#Pubmed ID']
                        species = SPECIES_DICT[header_dict['#Organism']]
                        score_type = 'score_type:' + header_dict['#Score Type']
                        data_list = [pheno_dict, PMID, species, score_type]

                        # Converting p-values into Z-scores
                        # (P-value - mean of scores from the whole entry) / standard deviation of scores from the whole entry

                        if header_dict['#Score Type'] == 'P-value' \
                                or header_dict['#Score Type'] == 'Fold inhibition (p-value)':
                            for gene in score_dict:
                                score_list.append(score_dict[gene])
                            for gene in score_dict:
                                Z_score = (score_dict[gene] - np.mean(score_list, axis=0)) \
                                                   / np.std(score_list, ddof=0)
                                score_dict[gene] = Z_score

                        # Creating matrix with scores
                        for gene in score_dict:
                            try:
                                # Mapping gene names to uniprot ids based on db created in mapper directory
                                c.execute(
                                    "SELECT UNIPROT_AC.uniprot_ac FROM MAPP LEFT JOIN UNIPROT_AC "
                                    "ON MAPP.uniprot_ac=UNIPROT_AC.id WHERE MAPP.foreign_id='%s'"
                                    % gene
                                )
                                while True:
                                    row = c.fetchone()
                                    if row is None:
                                        break
                                    else:
                                        mapped_genes_dict[gene] = row[0]
                            # TODO what is this???
                            except sqlite3.OperationalError:
                                print('sqliterror' + gene)
                                pass

                            # Standardizing Z scores (see in description)
                            try:
                                score_value = float(score_dict[gene])
                                if score_value >= 1:
                                    score_value = 1
                                elif score_value <= -1:
                                    score_value = -1
                                else:
                                    score_value = 0
                                if gene in mapped_genes_dict:
                                    df.at[mapped_genes_dict[gene], pheno_dict[gene]] = score_value
                            # If score is a valid number
                            except ValueError:
                                error_count += 1
                                pass
    #print(df)

    # SIGN PREDICTION
    conn = sqlite3.connect(BUILDER_SET)
    layer_list = [0, 1, 2, 3, 5, 6, 7]
    with conn:
        for layer in layer_list:
            logging.debug('Adding scores to layer' + str(layer))

            c = conn.cursor()
            c2 = conn.cursor()
            c.execute("SELECT interactor_a_node_name, interactor_b_node_name FROM layer%d"
                      % layer)
            while True:
                edge = c.fetchone()
                if edge is None:
                    break
                else:
                    #print(edge)
                    # Selecting score value for interacting proteins
                    pos_corr = 0
                    neg_corr = 0
                    A_node = edge[0].split(':')[1]
                    B_node = edge[1].split(':')[1]
                    # For each phenotype of interacting genes
                    for index in df:
                        try:
                            A_node_score = df._get_value(A_node, index)
                            B_node_score = df._get_value(B_node, index)
                        except KeyError:
                            error_count += 1
                            continue
                        # Positive correlation
                        if (A_node_score != 0 and B_node_score != 0) and A_node_score == B_node_score:
                            pos_corr += 1
                        # Negative correlation
                        elif (A_node_score != 0 and B_node_score != 0) and A_node_score != B_node_score:
                            neg_corr += 1
                        #  Counting instances where node scores are 0
                        else:
                            null_score_count += 1

                        # Calculating sign score according to formula
                        if pos_corr + neg_corr >= 1:
                            pheno_match_counter = pos_corr + neg_corr
                            weight_factor = math.sqrt(pheno_match_counter)
                            sign_score = ((pos_corr - neg_corr) / pheno_match_counter) * weight_factor
                            #print(sign_score)
                            # Adding scores to our dataset
                            c2.execute("UPDATE layer%d SET confidence_scores = confidence_scores || '%s'"
                                       "WHERE layer%d.interactor_a_node_name = '%s' AND layer%d.interactor_b_node_name = '%s'"
                                       % (layer, '|sign_pred:' + str(sign_score) + '|', layer, edge[0], layer, edge[1]))
                            notcut_count += 1
                        else:
                            cutoff_count += 1

    logging.debug('Sign prediction done')
    logging.debug('good', str(good_score_count))
    logging.debug('null score', str(null_score_count))
    logging.debug('cutoff allowed', str(notcut_count))
    logging.debug('cutoff_cut', str(cutoff_count))
    logging.debug('not P or Z score', str(not_searched_score_count))
    logging.debug('error', str(error_count))
    logging.debug('np or sp', str(no_score_count))
