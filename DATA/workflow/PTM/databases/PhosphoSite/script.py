'''
Parses Phosphosite Kinase-substarate data
 :argument: DATA_FILE: http://www.phosphosite.org/staticDownloads.action: Kinase_Substrate_Dataset.gz files
'''

from SLKlib.SQLiteDBApi.sqlite_db_api import PsimiSQL

# Defining constants
SQL_SEED = '../../../../../SLKlib/SQLiteDBApi/network-db-seed.sql'
DB_TYPE = 'PhosphoSite'
EXPORT_DB_LOCATION = '../../output/PhosphoSite'
DATA_FILE = 'files/Kinase_Substrate_Dataset'


def insert_or_get_node_dict(id, idtype, taxid, node_names_to_id, db_api):
    node_dict = {
        "name": idtype.strip() + ':' + id.strip(),
        "tax_id": taxid,
        "alt_accession": None,
        'pathways': None,
        "aliases": None,
        "topology": None
    }

    if node_dict['name'] in node_names_to_id:
        node_dict['id'] = node_names_to_id[node_dict['name']]
    else:
        db_api.insert_unique_node(node_dict)
        node_names_to_id[node_dict['name']] = node_dict['id']

    return node_dict

def main(logger):
    # Initiating the parser
    db_api = PsimiSQL(SQL_SEED)
    node_names_to_id = {}

    with open(DATA_FILE, encoding='ISO-8859-1') as data:

        # Skipping the header
        data.readline()
        data.readline()
        data.readline()
        data.readline()

        for line in data:
            columns = line.split('\t')

            if columns[3] == 'human' and columns[8] == 'human':
                taxid = 'taxid:9606'

                # Creating the node dicts, if the node is already in the db assigning that to the node dict
                source_dict = insert_or_get_node_dict(columns[2], "Uniprot", taxid, node_names_to_id, db_api)
                target_dict = insert_or_get_node_dict(columns[6], "Uniprot", taxid, node_names_to_id, db_api)

                # Nodes are inserted to the db if they are not in it yet
                if not 'id' in source_dict:
                    db_api.insert_node(source_dict)

                if not 'id' in target_dict:
                    db_api.insert_node(target_dict)

                interaction_types = "%s|is_directed:%s|is_direct:%s" \
                                    % ('MI:0217(phosphorylation)', "true", 'false')

                edge_dict = {
                    'publication_ids': 'pubmed:22135298',
                    'layer': '2',
                    'source_db': DB_TYPE,
                    'interaction_identifiers': None,
                    'confidence_scores': None,  # if available
                    'interaction_detection_method': None,  # probably exp type
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
