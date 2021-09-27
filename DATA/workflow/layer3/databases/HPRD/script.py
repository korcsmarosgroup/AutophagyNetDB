"""
Parses HPRD data files
    :argument: INTERACTION_FILE:
    :argument: EXPORT_DB_LOCATION: saving location
    :argument: DB_TYPE: name of the database
"""

# Imports
from SLKlib.SQLiteDBApi.sqlite_db_api import PsimiSQL
import re

# Defining constants
SQL_SEED = '../../../../../../SLKlib/SQLiteDBApi/network-db-seed.sql'
INTERACTION_FILE = '../HPRD/files/HPRD_Release9_062910/BINARY_PROTEIN_PROTEIN_INTERACTIONS.txt'
EXPORT_DB_LOCATION = '../../output/HPRD'
DB_TYPE = 'hprd'


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

    with open(INTERACTION_FILE) as intfile:
        node_names_to_id = {}
        lines = 0

        for line in intfile.readlines():
            columns = line.split('\t')
            lines += 1
            if lines % 50000 == 0:
                print("processed lines (HPRD): %d" % lines)


            # Creating the node dicts, if the node is already in the db assigning that to the node dict
            source_dict = insert_or_get_node_dict(columns[2], 'RefSeq', "taxid:9606", node_names_to_id, db_api)
            target_dict = insert_or_get_node_dict(columns[5], 'RefSeq', "taxid:9606", node_names_to_id, db_api)

            # Interaction detection method
            det_methods = map(lambda x: x.strip().lower(), columns[6].split(';'))
            det_dict = {
                'yeast 2-hybrid': 'MI:2254(two hybrid)',
                'in vitro': 'MI:0045(experimental interaction detection)',
                'in vivo': 'MI:0045(experimental interaction detection)',
                'y2h': 'MI:2254(two hybrid)',
                'vv': 'MI:0045(experimental interaction detection)',
                'vt': 'MI:0045(experimental interaction detection)',
                }
            det_methods = map(lambda x: "" if x not in det_dict else det_dict[x], det_methods)
            det_methods = filter(lambda x: x, det_methods)
            det_methods = set(det_methods)

            interaction_types = "is_directed:false|is_direct:true"

            pubmed_ids = map(lambda x: x.strip(), columns[7].split(','))
            pubmed_ids = filter(lambda x: re.search("^\\d+$", x), pubmed_ids)
            pubmed_ids = set(pubmed_ids)
            pubmed_ids.add("14525934")  # original HPRD publication
            pubmed_ids = map(lambda x: "pubmed:" + x, pubmed_ids)

            edge_dict = {
                'publication_ids': '|'.join(pubmed_ids),
                'layer': '3',
                'source_db': 'HPRD',
                'interaction_identifiers': None,
                'confidence_scores': None,
                'interaction_detection_method': '|'.join(det_methods),
                'interaction_types': interaction_types,
                'first_author': None
            }

            db_api.insert_edge(source_dict, target_dict, edge_dict)

        print("processed lines: %d" % lines)

    # Saving the to a DB_TYPE.db files
    db_api.save_db_to_file(EXPORT_DB_LOCATION)


if __name__ == '__main__':
    print("Parsing database...")
    main(logger = None)
    print("Parsing database is completed. SQLite database is saved to: " + EXPORT_DB_LOCATION)







