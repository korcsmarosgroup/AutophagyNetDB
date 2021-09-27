'''
 :argument: DATA_FILE: from http://www.cell.com/trends/cell-biology/fulltext/S0962-8924(09)00271-2 article: potential scaffold proteins
'''

from SLKlib.SQLiteDBApi.sqlite_db_api import PsimiSQL

# Defining constants
SQL_SEED = '../../../../../SLKlib/SQLiteDBApi/network-db-seed.sql'
DB_TYPE = 'PSP'
EXPORT_DB_LOCATION = '../../output/PSP'
DATA_FILE = 'files/mmc2.txt'


def get_node_a(name, taxid, psi_mi_to_sql_object):
    """
    This function sets up a node dict and returns it. If the node is already in the SQLite database it fetches that node from the db, so it won't be inserted multiple times
    """

    # Testing if the node is already in the database
    node_dict = psi_mi_to_sql_object.get_node(name, node_tax_id=taxid)

    if not node_dict:
        node_dict = {
            "name": 'Uniprot:' + name,
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
            "name": 'Uniprot:' + name,
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

    with open(DATA_FILE) as data:

        # Skipping the header
        data.readline()

        for line in data:
            columns = line.split('\t')
            taxid = 'taxid:9606'

            # Creating the node dicts, if the node is already in the db assigning that to the node dict
            source_dict = get_node_a(columns[0], taxid, db_api)
            target_dict = get_node_b(columns[2], taxid, db_api)

            # Nodes are inserted to the db if they are not in it yet
            if not 'id' in source_dict:
                db_api.insert_node(source_dict)

            if not 'id' in target_dict:
                db_api.insert_node(target_dict)

            interaction_types = "%s|is_directed:%s|is_direct:%s" \
                                % ('MI:0190(interaction type)', "false", 'false')
            edge_dict = {
                'publication_ids': 'pubmed:20005715',
                'layer': '1',
                'source_db': DB_TYPE,
                'interaction_identifiers': None,
                'confidence_scores': None,
                'interaction_detection_method': None,
                'interaction_types': interaction_types,
                'first_author': None
            }

            db_api.insert_edge(source_dict, target_dict, edge_dict)

    # Saving the to a DB_TYPE.db files
    db_api.save_db_to_file(EXPORT_DB_LOCATION)


if __name__ == '__main__':
    print("Parsing database...")
    main(logger=None)
    print("Parsing database is completed. SQLite database is saved to: " +
          EXPORT_DB_LOCATION)
