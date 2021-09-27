"""
 Tissue expression data from: https://bgee.org/?page=download&action=proc_values#id5
"""

# Imports
import logging
import os
import sqlite3
import json
import io

# # Initiating logger
# logger = logging.getLogger()
# handler = logging.FileHandler('../../workflow/SLK3.log')
# logger.setLevel(logging.DEBUG)
# handler.setLevel(logging.DEBUG)
# formatter = logging.Formatter(
#     '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
# handler.setFormatter(formatter)
# logger.addHandler(handler)


def tissue_prediction(logger):

    # INFILE_LIST = ['files/Homo_sapiens_expr_simple_development.tsv',
    #                'files/Caenorhabditis_elegans_expr_simple_development.tsv',
    #                'files/Drosophila_melanogaster_expr_simple_development.tsv',
    #                'files/Danio_rerio_expr_simple_development.tsv']
    INFILE_LIST = ['../prediction/tissue/files//Homo_sapiens_expr_simple_development.tsv',
                    '../prediction/tissue/files/Caenorhabditis_elegans_expr_simple_development.tsv',
                    '../prediction/tissue/files/Drosophila_melanogaster_expr_simple_development.tsv',
                    '../prediction/tissue/files/Danio_rerio_expr_simple_development.tsv']

    # Initiating mapper
    logger.debug('Initiating mapper')
    #map_conn = sqlite3.connect('../../workflow/mapper.db')
    map_conn = sqlite3.connect('mapper.db')

    tempfile = io.StringIO()
    for line in map_conn.iterdump():
        tempfile.write('%s\n' % line)
    map_conn.close()
    tempfile.seek(0)
    # Create a database in memory and import from tempfile
    map_db = sqlite3.connect(":memory:")
    with map_db:
        map_db.cursor().executescript(tempfile.read())
        map_db.cursor().execute("CREATE INDEX mapp_foreign ON MAPP(foreign_id);")
        map_db.cursor().execute("CREATE INDEX mapp_uniprot ON MAPP(uniprot_ac);")
        map_db.cursor().execute("CREATE INDEX mapp_gene_name ON MAPP(gene_name);")
        map_db.cursor().execute("CREATE INDEX uniprot_id ON UNIPROT_AC(id);")

    for file in INFILE_LIST:
        exp_dict = {}
        with map_db:
            c2 = map_db.cursor()
            with open(file) as infile:
                infile.readline()
                for line in infile:
                    line = line.strip().split('\t')
                    genename = line[1].replace('"', '')
                    # creating tissue identifier by concatenating ids and names
                    # we get this: CL:XXXXXXX(tissue name)
                    tissue = line[2] + '(' + line[3].replace('"', '') + ')'
                    is_present = line[6]

                    # we only add the tissue data if its expression is present in the current tissue
                    if is_present == 'present':
                        # mapping ENS gene ids to uniprot ids
                        try:
                            c2.execute(
                                "SELECT UNIPROT_AC.uniprot_ac FROM MAPP LEFT JOIN UNIPROT_AC ON MAPP.uniprot_ac=UNIPROT_AC.id "
                                "WHERE MAPP.gene_name='%s' GROUP BY MAPP.gene_name"
                                % genename
                            )
                            while True:
                                allrows = c2.fetchone()
                                if allrows is None:
                                    break
                                else:
                                    uniprot_id = allrows[0]
                                    if uniprot_id not in exp_dict:
                                        exp_dict[uniprot_id] = [tissue]
                                    else:
                                        # Adding data if it's not already in the dictionary
                                        if tissue not in exp_dict[uniprot_id]:
                                            exp_dict[uniprot_id].append(tissue)
                        except Exception as e:
                            print("### ERROR tissue prediction:")
                            print(e)
                            print("### /")

        # Connecting to builder node table, adding data
        #conn2 = sqlite3.connect('../../workflow/SLK3_layers.db')
        conn2 = sqlite3.connect('SLK3_layers.db')
        with conn2:
            c3 = conn2.cursor()
            c4 = conn2.cursor()

            for key, value in exp_dict.items():
                try:
                    tissues = ['tissue:%s' % t for t in value]
                    topology = ""
                    uniprotac = "Uniprot:%s" % key
                    c3.execute(
                        "SELECT node.topology FROM node WHERE node.name = '%s'" % uniprotac)
                    while True:
                        row = c3.fetchone()
                        if row is None:
                            break
                        else:
                            if row[0] != "":
                                # Appending tissue data to already existing data in topology column
                                topology = "%s|%s" % (row[0], '|'.join(tissues))
                            else:
                                topology = '|'.join(tissues)
                            c4.execute("UPDATE node SET topology = ? WHERE node.name = ?", (
                                topology, uniprotac))

                except Exception as e:
                    print("### ERROR tissue prediction:")
                    print(e)
                    print("### /")


if __name__ == '__main__':
    tissue_prediction(None)
