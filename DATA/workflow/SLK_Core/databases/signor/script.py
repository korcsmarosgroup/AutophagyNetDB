"""
    Creating signor database
        :argument: DB_TYPE: name of the source database
        :argument: DB_DESTINATION: saving location of the created database files
        :argument: CSV_LIST: list of .csv files of each signalink pathway
        :argument: FILENAME_TO_PATHWAY_MAP: dictionary from files name to SLK pathway
        :argument: IS_DIRECT_MAP: dictionary of directness to MI ids
        :argument: EFFECT_MAP: dictionary of effect to MI ids
        :argument: MOLECULAR_MAP: dictionary of molecular background to MI ids

    Important!
    Since the signor db provides it's data in multiple files. This script can be called multiple times to generate
    a single SQL .db files. On the first run the script creates a signore.db files. If the script is called again on
    another signor .tsv files, it can extend the previously created signor.db files by adding 'signor.db' as a fourth argument.
"""

# -*- coding: utf-8 -*-

import csv, sys

from ARNlib.SQLiteDBApi.sqlite_db_api import PsimiSQL

#Defining constants
SQL_SEED = '../../../../../ARNlib/SQLiteDBApi/network-db-seed.sql'
DB_TYPE = 'Signor'
DB_DESTINATION = '../../../layer0/output/signor'
CSV_LIST = ['files/SIGNOR_WNT.csv',
            'files/SIGNOR_TOLLR.csv',
            'files/SIGNOR_TGFb.csv',
            'files/SIGNOR_SAPK_JNK.csv',
            'files/SIGNOR_P38.csv',
            'files/SIGNOR_NOTCH.csv',
            'files/SIGNOR_NFKBNC.csv',
            'files/SIGNOR_NFKBC.csv',
            'files/SIGNOR_MTOR.csv',
            'files/SIGNOR_MCAPO.csv',
            'files/SIGNOR_INSR.csv',
            'files/SIGNOR_IL1R.csv',
            'files/SIGNOR_IAPO.csv',
            'files/SIGNOR_AMPK.csv',
            'files/SIGNOR_BMP.csv',
            'files/SIGNOR_DR.csv',
            'files/SIGNOR_EGF.csv',
            'files/SIGNOR_HPP.csv',
            'files/SIGNOR_Hedgehog.csv' ]

NUMBER_OF_FILES = len(CSV_LIST)
FILENAME_TO_PATHWAY_MAP = {
    'SIGNOR_AMPK.csv': 'Receptor tyrosine kinase',
    'SIGNOR_BMP.csv': 'TGF',
    'SIGNOR_DR.csv': 'TNF pathway',
    'SIGNOR_EGF.csv': 'Receptor tyrosine kinase',
    'SIGNOR_HPP.csv': 'Hippo',
    'SIGNOR_Hedgehog.csv': 'Hedgehog',
    'SIGNOR_IAPO.csv': 'TNF pathway',
    'SIGNOR_IL1R.csv': 'JAK/STAT',
    'SIGNOR_INSR.csv': 'Receptor tyrosine kinase',
    'SIGNOR_MCAPO.csv': 'TNF pathway',
    'SIGNOR_MTOR.csv': 'Receptor tyrosine kinase',
    'SIGNOR_NFKBC.csv': 'Innate immune pathways',
    'SIGNOR_NFKBNC.csv': 'Innate immune pathways',
    'SIGNOR_NOTCH.csv': 'Notch',
    'SIGNOR_P38.csv': 'Receptor tyrosine kinase',
    'SIGNOR_SAPK_JNK.csv': 'Receptor tyrosine kinase',
    'SIGNOR_TGFb.csv': 'TGF',
    'SIGNOR_TOLLR.csv': 'Toll-like receptor',
    'SIGNOR_WNT.csv': 'WNT/Wingless'
}

IS_DIRECT_MAP = {
    "YES": "MI:0407(direct interaction)",
    "NO": "indirect",
    "UNK": "unknown"
}

EFFECT_MAP = {
    'Unknown': 'MI:0190(interaction type)',
    'down-regulates': 'MI:2240(down-regulates)',
    "down-regulates activity": 'MI:2241(down-regulates activity)',
    "down-regulates quantity": 'MI:2242(down-regulates quantity)',
    "down-regulates quantity by destabilization": 'MI:2244(down-regulates quantity by destabilization)',
    "down-regulates quantity by repression": 'MI:2243(down-regulates quantity by repression)',
    'unknown': 'MI:0190(interaction type)',
    'up-regulates': 'MI:2235(up-regulates)',
    "up-regulates activity": 'MI:2236(up-regulates activity)',
    "up-regulates quantity by expression": 'MI:2238(up-regulates quantity by expression)',
    "up-regulates quantity by stabilization": 'MI:2239(up-regulates quantity by stabilization)'
}

MOLECULAR_MAP = {
    'binding' : 'MI:0462(bind)',
    'transcriptional regulation' : 'MI:2247(transcriptional regulation)',
    'phosphorylation' : 'MI:0217(phosphorylation reaction)',
    '' : 'MI:0190(interaction type)',
    'ubiquitination' : 'MI:0220(ubiquitination reaction)',
    'relocalization' : 'MI:2256(relocalization reaction)',
    'dephosphorylation' : 'MI:0203(dephosphorylation reaction)',
    'cleavage' : 'MI:0194(cleavage reaction)',
    'deubiquitination' : 'MI:0204(deubiquitination reaction)',
    'guanine nucleotide exchange factor' : 'MI:2252(guanine nucleotide exchange factor)'
}


def main(logger):
    # Initiating a PsimiSQL class
    db_api = PsimiSQL(SQL_SEED)

    # Making the script user friendly
    file_counter = 1
    print("Started parsing .csv files")

    # Parsing data files
    for csv_file_location in CSV_LIST:

        csv_file_name = csv_file_location.split('/')[-1]

        sys.stdout.write("Parsing '%s' (%d/%d)\r" % (csv_file_name, file_counter, NUMBER_OF_FILES))

        csv_file = csv.reader(open(csv_file_location, encoding="ISO-8859-1"), delimiter = ';', quotechar = '"')

        pathway = FILENAME_TO_PATHWAY_MAP[csv_file_name]

        # Skipping the header
        for cells in csv_file:

                type_a = cells[1].lower()
                type_b = cells[5].lower()

                taxids = cells[12].split(';')[0]

                if type_a == 'protein' and type_b == 'protein' and taxids == '9606':

                    # Dealing with the first node

                    node_a_name = f'Uniprot:{cells[2]}'
                    node_a_taxid = 'taxid:' + taxids
                    node_a_taxid = node_a_taxid

                    node_a_dict = {}

                    # If the node already exists in the db, than only it's pathway will be modified, otherwise it will be added to the db
                    if db_api.get_node(node_a_name,node_a_taxid):
                        node_a_dict = db_api.get_node(node_a_name,node_a_taxid)
                        if not pathway in node_a_dict['pathways']:
                            node_a_dict['pathways'] += '|'+pathway
                            db_api.update_node(node_a_dict)
                    else:
                        node_a_dict = {
                            'name' : node_a_name,
                            'alt_accession' : 'entrez gene/locuslink:'+cells[0],
                            'tax_id' : node_a_taxid,
                            'pathways' : pathway,
                            'aliases' : None,
                            'topology' : ""
                        }
                        db_api.insert_node(node_a_dict)

                    # Doing the same with node b

                    node_b_name = f'Uniprot:{cells[2]}'
                    node_b_taxid = 'taxid:' + taxids
                    node_b_taxid = node_b_taxid

                    node_b_dict = {}

                    # If the node already exists in the db, than only it's pathway will be modified, otherwise it will be added to the db

                    if db_api.get_node(node_b_name,node_b_taxid):
                        node_b_dict = db_api.get_node(node_b_name,node_b_taxid)
                        if not pathway in node_b_dict['pathways']:
                            node_b_dict['pathways'] += '|'+pathway
                            db_api.update_node(node_b_dict)
                    else:
                        node_b_dict = {
                            'name' : node_b_name,
                            'alt_accession' : 'entrez gene/locuslink:'+cells[4],
                            'tax_id' : node_b_taxid,
                            'pathways' : pathway,
                            'aliases' : None,
                            'topology' : ""
                        }
                        db_api.insert_node(node_b_dict)

                    # Getting publication id
                    publication_id = ['pubmed:'+cells[21]]
                    publication_id.append("pubmed:26467481")

                    effect = EFFECT_MAP[cells[8]]

                    molecular_background = MOLECULAR_MAP[cells[9]]

                    inttype_final = effect + '|' + molecular_background

                    is_direct = IS_DIRECT_MAP[cells[22]].strip()
                    if "MI:0407(direct interaction)" in is_direct:
                        is_direct = "true"
                    else:
                        is_direct = "false"

                    # Setting up the interaction type
                    interaction_types = "%s|is_directed:%s|is_direct:%s" \
                                        % (inttype_final, "true", is_direct)

                    edge_dict = {
                        'interaction_detection_method': None,
                        'first_author': None,
                        'publication_ids': "|".join(publication_id),
                        'interaction_types': interaction_types,
                        'source_db': 'Signor',
                        'interaction_identifiers': None,
                        'confidence_scores': None,
                        'layer': "0"
                    }

                    db_api.insert_edge(node_a_dict,node_b_dict,edge_dict)

    print("Parsing files finished!")
    print("Finished parsing Signor. Saving db to %s.db" % (DB_TYPE))
    db_api.save_db_to_file(DB_DESTINATION)


if __name__ == '__main__':
    main()
