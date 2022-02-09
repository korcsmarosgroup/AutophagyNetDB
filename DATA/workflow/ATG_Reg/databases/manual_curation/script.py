"""
  Importing manual curation data from ARN v1
"""

# Imports
import logging
from ARNlib.SQLiteDBApi.sqlite_db_api import PsimiSQL

# Defining constants
SQL_SEED = '../../../../ARNlib/SQLiteDBApi/network-db-seed.sql'
DATA_FILE_LIST = ['files/Autofágia Regulatory Network - TD-nek_2013. 09. 26..txt',
                  'files/Autofágia Regulatory Network - TD-nek_v2.txt']
EXPORT_DB_DESTINATION = '../output/manualcur'
DB_TYPE = 'manual curation'

# Initiating logger
logger = logging.getLogger()
handler = logging.FileHandler('../../SLK3.log')
logger.setLevel(logging.DEBUG)
handler.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)


def get_node_a(id, taxid, psi_mi_to_sql_object):
    """
    This function sets up a node dict and returns it.
    If the node is already in the SQLite database it fetches that node from the db, so it won't be inserted multiple times.
    """

    # Testing if the node is already in the database

    node_dict = psi_mi_to_sql_object.get_node(id, node_tax_id=taxid)

    if not node_dict:
        node_dict = {
            "name": id,
            "tax_id": taxid,
            "alt_accession": None,
            'pathways': None,
            "aliases": None,
            "topology": None
        }

    return node_dict


def get_node_b(id, taxid, psi_mi_to_sql_object):
    """
    This function sets up a node dict and returns it. If the node is already in the SQLite database it fetches that node from the db, so it won't be inserted multiple times.

    """

    # Testing if the node is already in the database
    node_dict = psi_mi_to_sql_object.get_node(id, node_tax_id=taxid)

    if not node_dict:
        node_dict = {
            "name": id,
            "tax_id": taxid,
            "alt_accession": None,
            'pathways': None,
            "aliases": None,
            "topology": None
        }

    return node_dict


def main(logger):
    # Initiating the parser
    db_api = PsimiSQL(SQL_SEED)

    # Parsing data file
    for file in DATA_FILE_LIST:
        print(file)
        with open(file, encoding="ISO-8859-1") as data:
            # Skipping header
            data.readline()
            data.readline()
            data.readline()
            data.readline()
            data.readline()

            for line in data:
                line = line.strip().split('\t')
                if len(line) < 4:    # Probably because of conversion from xlsx to tsv
                    continue
                # Mapping species to taxid
                if line[0] == 'human':
                    taxid_source = 'taxid:9606'
                else:
                    pass
                if line[2] == 'human':
                    taxid_target = 'taxid:9606'
                else:
                    pass

                # Creating the node dicts, if the node is already in the db assigning that to the node dict
                source_dict = get_node_a('Uniprot:' + line[1], taxid_source, db_api)
                target_dict = get_node_b('Uniprot:' + line[3], taxid_target, db_api)

                # Nodes are inserted to the db if they are not in it yet
                if not 'id' in source_dict:
                    db_api.insert_node(source_dict)

                if not 'id' in target_dict:
                    db_api.insert_node(target_dict)

                # Mapping interaction identifiers
                # Directed/undirected
                if line[5] == 'D':
                    directed = 'directed'
                else:
                    directed = 'undirected'
                # Direct/indirect
                if line[7] == 'D':
                    direct = 'MI:0407(directed)'
                else:
                    direct = 'MI:2246(indirect)'
                # Stimulation/inhibition
                if line[8] == 'S' or line[8] == 's':
                    stimulation = 'MI:0840(stimulator)'
                elif line[8] == 'I':
                    stimulation = 'MI:0586(inhibitor)'
                else:
                    pass
                # Molecular background
                molec_map = {
                    'P': 'MI:0217(phosphorylation reaction)',
                    'Acetylation': 'MI:0192(acetylation reaction)',
                    'degradation (ubiquitinilation)': 'MI:0220(ubiquitination reaction)',
                    'autoP': 'MI:0217(phosphorylation reaction)',
                    'csak beköt': 'MI:0462(bind)',
                    'proteolízis': 'MI:0414(enzymatic reaction)',
                    'proteolízis ("delipidálás")': 'MI:0414(enzymatic reaction)',
                    '"proteolízis (""delipidálás"")"': 'MI:0414(enzymatic reaction)',
                    'E2 - kovalens tioészter kötés': 'MI:0195(covalent binding)',
                    'kovalens': 'MI:0195(covalent binding)',
                    'kovalens tioészter kötés': 'MI:0195(covalent binding)',
                    'E1 - kovalens tioészter kötés': 'MI:0195(covalent binding)',
                    'E1-E2 komplex': 'MI:0195(covalent binding)',
                    '': ''
                }

                # Constructing interaction data line
                int_types = '|'.join([stimulation, molec_map[line[9]],
                                      'is_direct:' + 'true', 'is_directed:' + 'true'])

                edge_dict = {
                    'publication_ids': 'pubmed:' + line[4],
                    'layer': '1',
                    'source_db': 'manual curation',
                    'interaction_identifiers': None,
                    'confidence_scores': None,
                    'interaction_detection_method': None,
                    'interaction_types': int_types,
                    'first_author': None
                }

                db_api.insert_edge(source_dict, target_dict, edge_dict)

    # Saving the to a DB_TYPE.db file
    db_api.save_db_to_file(EXPORT_DB_DESTINATION)


if __name__ == '__main__':
    print("Parsing database...")
    main(logger=None)
    print("Parsing database is completed. SQLite database is saved to: " + EXPORT_DB_DESTINATION)
