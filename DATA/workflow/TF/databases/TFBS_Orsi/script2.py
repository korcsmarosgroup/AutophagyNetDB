"""
  Parses TF data from Orsi
"""

from ARNlib.SQLiteDBApi.sqlite_db_api import PsimiSQL

# Defining constants
SQL_SEED = '../../../../../ARNlib/SQLiteDBApi/network-db-seed.sql'
DB_TYPE = 'TFlink'
EXPORT_DB_LOCATION = '../../output/TFBS_Orsi'
DATA_FILE = 'files/allInteractionData_v3_fixed.tsv'


def insert_or_get_node_dict(id, taxid, node_names_to_id, db_api):

    node_dict = {
        "name": f'Uniprot:{id.strip()}',
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

        node_names_to_id = {}

        for line in data:
            line = line.strip().split('\t')
            taxid_a = line[8].split('(')[0]
            taxid_b = line[9].split('(')[0]

            accepted_tax_ids = ['9606', '7227', '6239', '7955']

            if (taxid_a.split(":")[1] in accepted_tax_ids) and (taxid_b.split(":")[1] in accepted_tax_ids):

                id_a = line[0].split(":")[1]
                id_b = line[1].split(":")[1]

                # Creating the node dicts, if the node is already in the db assigning that to the node dict
                source_dict = insert_or_get_node_dict(id_a, taxid_a, node_names_to_id, db_api)
                target_dict = insert_or_get_node_dict(id_b, taxid_b, node_names_to_id, db_api)

                interaction_types = "%s|is_directed:%s|is_direct:%s" \
                                    % ('MI:2247(trascriptional regulation)', 'true', 'false')

                if 'pubmed' in line[7]:
                    pub = []
                    pubmed_ids = line[7].split("|")
                    for pubmed_id in pubmed_ids:
                        if pubmed_id != '':
                            id = pubmed_id.split(":")[1]
                            if "e" in id:
                                id = id.replace("e", "")
                            new = f'pubmed:{id}'
                            pub.append(new)
                else:
                    pub = ""

                interaction_detection_methods = line[6].strip().split("|")
                new_interaction_detection_methods = []
                for method in interaction_detection_methods:
                    if "MI:" in method:
                        m = method.split(":")[2].replace('"', '')
                        new_method = f'MI:{m}'
                        new_interaction_detection_methods.append(new_method)

                edge_dict = {
                    'publication_ids': "|".join(pub),
                    'layer': '6',
                    'source_db': DB_TYPE,
                    'interaction_identifiers': None,
                    'confidence_scores': None,
                    'interaction_detection_method': "|".join(new_interaction_detection_methods),
                    'interaction_types': interaction_types,
                    'first_author': None
                }

                db_api.insert_edge(source_dict, target_dict, edge_dict)

    # Saving the to a DB_TYPE.db files
    db_api.save_db_to_file(EXPORT_DB_LOCATION)


if __name__ == '__main__':
    print("Parsing database...")
    main(logger=None)
    print("Parsing database is completed. SQLite database is saved to: " + EXPORT_DB_LOCATION)
