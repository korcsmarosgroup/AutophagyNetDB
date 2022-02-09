'''
 :argument: DATA_FILE: from SLK 2.0 endocytosis data
'''

from ARNlib.SQLiteDBApi.sqlite_db_api import PsimiSQL

# Defining constants
SQL_SEED = '../../../../../ARNlib/SQLiteDBApi/network-db-seed.sql'
DB_TYPE = 'SignaLink3'
EXPORT_DB_LOCATION = '../../output/signalink3'
DATA_FILE = 'files/SLK3_human_core.csv'


def get_node_a(name, taxid, alt_acc, pathway, topology, psi_mi_to_sql_object):
    """
    This function sets up a node dict and returns it. If the node is already in the SQLite database it fetches that node from the db, so it won't be inserted multiple times
    """

    # Testing if the node is already in the database
    node_dict = psi_mi_to_sql_object.get_node(name, node_tax_id=taxid)

    if not node_dict:
        node_dict = {
            "name" : 'Uniprot:' + name,
            "tax_id": taxid,
            "alt_accession": alt_acc,
            'pathways': pathway,
            "aliases": None,
            "topology": topology
        }

    return node_dict


def get_node_b(name, taxid, alt_acc, pathway, topology, psi_mi_to_sql_object):
    """
    This function sets up a node dict and returns it. If the node is already in the SQLite database it fetches that node from the db, so it won't be inserted multiple times.

    """

    # Testing if the node is already in the database
    node_dict = psi_mi_to_sql_object.get_node(name, node_tax_id=taxid)

    if not node_dict:
        node_dict = {
            "name": 'Uniprot:' + name,
            "tax_id": taxid,
            "alt_accession": alt_acc,
            'pathways': pathway,
            "aliases": None,
            "topology": topology
        }

    return node_dict


def main(logger):
    # Initiating the parser
    db_api = PsimiSQL(SQL_SEED)

    with open(DATA_FILE) as data:

        # Skipping the header
        data.readline()

        for line in data:
            columns = line.strip().split(',')

            # Creating the node dicts, if the node is already in the db assigning that to the node dict
            source_dict = get_node_a('Uniprot:' + columns[1], 'taxid:' + columns[2], columns[0], columns[5], columns[4], db_api)
            target_dict = get_node_b('Uniprot:' + columns[7], 'taxid:' + columns[6], columns[8], columns[11], columns[10], db_api)

            # Nodes are inserted to the db if they are not in it yet
            if not 'id' in source_dict:
                db_api.insert_node(source_dict)

            if not 'id' in target_dict:
                db_api.insert_node(target_dict)

            # Pubmed references
            pub_id = '|pubmed:'.join(columns[16].split('|'))

            # Directedness
            effect = columns[13]

            interaction_types = "%s|is_directed:truw|is_direct:true" \
                                % effect

            edge_dict = {
                'publication_ids': columns[16],
                'layer': '3',
                'source_db': 'SignaLink3',
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
    main(logger = None)
    print("Parsing database is completed. SQLite database is saved to: " + EXPORT_DB_LOCATION)



