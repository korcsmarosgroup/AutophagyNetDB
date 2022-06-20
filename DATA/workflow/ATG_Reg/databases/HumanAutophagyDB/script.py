"""
  Parses the mapped Human autophagy database list and finds any edges from third party databases
  if they contain the listed genes
  :param merger = merger database file THIS USES UNIPROT IDS
"""

import sqlite3
from ARNlib.SQLiteDBApi.sqlite_db_api import PsimiSQL

HA_DATA_FILE = 'HAdb_mapped.txt'
ATG_DATA_FILE = 'atg_mapped.txt'
l3_file = 'files/ARN2_merged_L1.csv'
EXPORT_DB_LOCATION = '../../output/HADB'


def get_node_a(name, taxid, psi_mi_to_sql_object):
    """
    This function sets up a node dict and returns it. If the node is already in the SQLite database it fetches that node from the db, so it won't be inserted multiple times
    """

    # Testing if the node is already in the database
    node_dict = psi_mi_to_sql_object.get_node(name, node_tax_id=taxid)

    if not node_dict:
        node_dict = {
            "name" : name,
            "tax_id": taxid,
            "alt_accession": None,
            'pathways': None,
            "aliases": None,
            "topology": None
        }

    return node_dict


def get_node_b(name, taxid, psi_mi_to_sql_object):
    """
    This function sets up a node dict and returns it. If the node is already in the SQLite database it fetches that node from the db, so it won't be inserted multiple times.

    """

    # Testing if the node is already in the database
    node_dict = psi_mi_to_sql_object.get_node(name, node_tax_id=taxid)

    if not node_dict:
        node_dict = {
            "name": name,
            "tax_id": taxid,
            "alt_accession": None,
            'pathways': None,
            "aliases": None,
            "topology": None
        }

    return node_dict

def main(logger):
    db_api = PsimiSQL(SQL_SEED)

    with open(HA_DATA_FILE) as ha_infile:
            ha_infile.readline()
            ATG_PROT_LIST = []
            for line in ha_infile:
                line = line.strip().split('\t')
                uniprot = line[1]
                if uniprot not in ATG_PROT_LIST:
                    ATG_PROT_LIST.append(uniprot)
    with open(ATG_DATA_FILE) as atg_infile:
            atg_infile.readline()
            for line in atg_infile:
                line = line.strip().split('\t')
                uniprot = line[1]
                if uniprot not in ATG_PROT_LIST:
                    ATG_PROT_LIST.append(uniprot)

    myrows = []
    with open(l3_file) as l3_edges:
        l3_edges.readline()
        for line in l3_edges:
            line = line.strip().split(',')
            source_uni = line[3].split(':')[1]
            target_uni = line[4].split(':')[1]
            for item in ATG_PROT_LIST:
                if item == target_uni:
                    myrows.append(line)

    for row in myrows:
        source = row[3]
        target = row[4]
        # Creating the node dicts, if the node is already in the db assigning that to the node dict
        source_dict = get_node_a(source, 'taxid:9606', db_api)
        target_dict = get_node_b(target, 'taxid:9606', db_api)

        # Nodes are inserted to the db if they are not in it yet
        if not 'id' in source_dict:
            db_api.insert_node(source_dict)

        if not 'id' in target_dict:
            db_api.insert_node(target_dict)

        # Pubmed references
        pub_id = row[7]
        interaction_types = row[8]

        edge_dict = {
            'publication_ids': pub_id,
            'layer': '1',
            'source_db': 'HumanAutophagyDB',
            'interaction_identifiers': None,
            'confidence_scores': row[11],
            'interaction_detection_method': None,
            'interaction_types': interaction_types,
            'first_author': None
        }

        db_api.insert_edge(source_dict, target_dict, edge_dict)

        # Saving the to a DB_TYPE.db files
    db_api.save_db_to_file(EXPORT_DB_LOCATION)






