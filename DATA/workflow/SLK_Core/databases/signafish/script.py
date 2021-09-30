'''
 :argument: DATA_FILE: from SLK 2.0 endocytosis data
'''

from ARNlib.SQLiteDBApi.sqlite_db_api import PsimiSQL

# Defining constants
SQL_SEED = '../../../../../ARNlib/SQLiteDBApi/network-db-seed.sql'
DB_TYPE = 'SignaFish'
EXPORT_DB_LOCATION = '../../output/signafish'
DATA_FILE = 'files/signafish-L0.csv'


def get_node_a(name, taxid, pathway, topology, psi_mi_to_sql_object):
    """
    This function sets up a node dict and returns it. If the node is already in the SQLite database it fetches that node from the db, so it won't be inserted multiple times
    """

    # Testing if the node is already in the database
    node_dict = psi_mi_to_sql_object.get_node(name, node_tax_id=taxid)

    if not node_dict:
        node_dict = {
            "name" : 'Uniprot:' + name,
            "tax_id": taxid,
            "alt_accession": None,
            'pathways': pathway,
            "aliases": None,
            "topology": topology
        }

    return node_dict


def get_node_b(name, taxid, pathway, topology, psi_mi_to_sql_object):
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
            columns = line.strip().split(';')
            taxid = 'taxid:7955'

            source_path = []
            for item in columns[5].split(','):
                if '(' in item:
                    sp = item.split('(')[0]
                else:
                    sp = item
                source_path.append(sp)

            source_path_final = '|'.join(source_path)

            target_path = []
            for item in columns[11].split(','):
                if '(' in item:
                    tp = item.split('(')[0]
                else:
                    tp = item
                source_path.append(tp)

            target_path_final = '|'.join(target_path)

            # Creating the node dicts, if the node is already in the db assigning that to the node dict
            source_dict = get_node_a(columns[1], taxid, source_path_final,'|'.join(columns[4].replace('d ', 'd').split(',')), db_api)
            target_dict = get_node_b(columns[7], taxid, target_path_final,'|'.join(columns[10].replace('d ', 'd').split(',')), db_api)

            # Nodes are inserted to the db if they are not in it yet
            if not 'id' in source_dict:
                db_api.insert_node(source_dict)

            if not 'id' in target_dict:
                db_api.insert_node(target_dict)

            # Pubmed references
            pub_id = '|pubmed:'.join(columns[16].split('|'))

            # Directedness
            if columns[14] == 'direct':
                isdirect = 'true'
            else:
                isdirect = 'false'

            if columns[13] == 'PPI directed':
                isdirected = 'true'
            else:
                isdirected = 'false'

            # Effect
            if columns[15] == 'stimulation':
                effect = 'MI:0624(stimulation)'
                interaction_types = "%s|is_directed:%s|is_direct:%s" \
                                    % (effect, isdirected, isdirect)
            elif columns[15] == 'inhibition':
                effect = 'MI:0623(inhibition)'

                interaction_types = "%s|is_directed:%s|is_direct:%s" \
                                    % (effect, isdirected, isdirect)
            else:
                interaction_types = "is_directed:%s|is_direct:%s" \
                                    % (isdirected, isdirect)
            edge_dict = {
                'publication_ids': 'pubmed:' + pub_id,
                'layer': '0',
                'source_db': 'SignaFish',
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



