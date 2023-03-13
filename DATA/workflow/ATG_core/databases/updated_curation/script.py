"""
2022 updated manual curation interaction data
"""

# Imports
from ARNlib.SQLiteDBApi.sqlite_db_api import PsimiSQL
import sqlite3, os

# Defining constants
SQL_SEED = '../../../../../ARNlib/SQLiteDBApi/network-db-seed.sql'
DATA_FILE = 'files/autophagynet_updated_curation_TH_LG.csv'
DB_DESTINATION = '../../../all_output/updated_curation.db'


def get_node_a(aid, taxid, psi_mi_to_sql_object):
    """
    This function sets up a node dict and returns it. If the node is already in the SQLite database it fetches that node from the db, so it won't be inserted multiple times.
    """

    # Testing if the node is already in the database
    node_dict = psi_mi_to_sql_object.get_node(aid, node_tax_id=taxid)

    if not node_dict:
        node_dict = {
            "name": aid,
            "tax_id": taxid,
            "alt_accession": None,
            'pathways': None,
            "aliases": None,
            "topology": None
        }

    return node_dict


def get_node_b(bid, taxid, psi_mi_to_sql_object):
    """
    This function sets up a node dict and returns it. If the node is already in the SQLite database it fetches that node from the db, so it won't be inserted multiple times.
    """

    # Testing if the node is already in the database
    node_dict = psi_mi_to_sql_object.get_node(bid, node_tax_id=taxid)

    if not node_dict:
        node_dict = {
            "name": bid,
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
    dest_conn = sqlite3.connect(DB_DESTINATION)

    # Parsing file
    with open(DATA_FILE, encoding='ISO-8859-1') as data:
        data.readline()
        for line in data:
            line = line.strip().split('\t')

            source_genename = line[0].strip()
            target_genename= line[1].strip()

            # Creating the node dicts, if the node is already in the db assigning that to the node dict
            source_dict = get_node_a('HGNC:' + source_genename, 'taxid:9606', db_api)
            target_dict = get_node_b('HGNC:' + target_genename, 'taxid:9606', db_api)

            # Nodes are inserted to the db if they are not in it yet
            if 'id' not in source_dict:
                db_api.insert_node(source_dict)

            if 'id' not in target_dict:
                db_api.insert_node(target_dict)

            # Layer mapping
            layer_map_dict = {
                'Direct regulation of autophagy': 1,
                'Interaction between autophagy proteins': 0,
                'Potentially new regulator proteins/Additional protein regulators': 2
            }
            # Effect mapping
            effect_map_dict = {
                'Unknown': 'unknown',
                "Activation": 'MI:0624(stimulation)',
                'Inhibition': 'MI:0623(inhibition)'

            }

            directedness = line[3].lower()
            effect = effect_map_dict[line[2].strip()]
            # Assembling line
            ident = ('effect:' + effect + '|is_direct:true' + '|is_directed:' + directedness)

            edge_dict = {
                'publication_ids': 'pubmed:' + line[7],
                'layer': layer_map_dict[line[4]],
                'source_db': 'manual curation',
                'interaction_identifiers': ident,
                'confidence_scores': None,
                'interaction_detection_method': None,
                'interaction_types': None,
                'first_author': None
            }

            db_api.insert_edge(source_dict, target_dict, edge_dict)

            # Saving the to a DB_TYPE.db file
        db_api.save_db_to_file(DB_DESTINATION)


# if __name__ == '__main__':
#     main(None)

