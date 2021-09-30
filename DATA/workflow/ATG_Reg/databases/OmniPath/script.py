"""
 Parse Omnipath pathway with tissue patterns as additional node information (www.omnipathdb.org/data)
    :argument: DATA_FILE: http://omnipathdb.org/interactions/?fields=sources&fields=references
    :argument: EXPORT_DB_LOCATION: saving location
    :argument: DB_TYPE: name of the database
    
"""

# Imports
from ARNlib.SQLiteDBApi.sqlite_db_api import PsimiSQL
import re

# Defining constants
SQL_SEED = '../../../../../ARNlib/SQLiteDBApi/network-db-seed.sql'
DATA_FILE = '../../../ATG_Reg/databases/OmniPath/files/interactions.txt'
EXPORT_DB_LOCATION = '../output/OmniPath'
DB_TYPE = 'omnipath'

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

    with open(DATA_FILE) as data:
        # Skipping the header
        data.readline()
        node_names_to_id = {}
        lines = 0

        for line in data:
            columns = line.strip().split('\t')
            lines += 1
            if lines % 50000 == 0:
                print("processed lines (OmniPath): %d" % lines)

            # Creating the node dicts, if the node is already in the db assigning that to the node dict
            source_dict = insert_or_get_node_dict(columns[0], "Uniprot", 'taxid:9606', node_names_to_id, db_api)
            target_dict = insert_or_get_node_dict(columns[1], "Uniprot", 'taxid:9606', node_names_to_id, db_api)

            # the link is indirect, unless it is directed or if it has a role as inhibitor or stimulator
            direct = 'false'

            if columns[2] == '1':
                directed = 'true'
                direct = 'true'
            elif columns[2] == '0':
                directed = 'false'
            else:
                print("WARNING: unknown direction flag in line: " + line)

            interaction_type_terms = []

            if columns[3] == '1':
                interaction_type_terms.append('MI:0624(stimulant)')
                direct = 'true'

            if columns[4] == '1':
                interaction_type_terms.append('MI:0623(inhibition)')
                direct = 'true'

            interaction_types = "is_directed:%s|is_direct:%s" % (directed, direct)
            if len(interaction_type_terms) > 0:
                interaction_types += "|" + "|".join(interaction_type_terms)


            pubmed_ids = map(lambda x: x.strip(), columns[7].split(';'))
            pubmed_ids = filter(lambda x: re.search("^\\d+$", x), pubmed_ids)
            pubmed_ids = set(pubmed_ids)
            pubmed_ids.add("27898060")  # OmniPath publication
            pubmed_ids = map(lambda x: "pubmed:" + x, pubmed_ids)

            edge_dict = {
                    'publication_ids': "|".join(pubmed_ids),
                    'layer': '3',
                    'source_db': 'OmniPath',
                    'interaction_identifiers': None,
                    'confidence_scores': None,
                    'interaction_detection_method': None,
                    'interaction_types': interaction_types,
                    'first_author': None
                }

            db_api.insert_edge(source_dict, target_dict, edge_dict)

        print("processed lines (OmniPath): %d" % lines)

    # Saving the to a DB_TYPE.db files
    db_api.save_db_to_file(EXPORT_DB_LOCATION)


if __name__ == '__main__':
    print("Parsing database...")
    main(logger = None)
    print("Parsing database is completed. SQLite database is saved to: " + EXPORT_DB_LOCATION)
