"""
Parsing miR2Disease data
 :argument: DATA_FILE: miRNA-target data files http://watson.compbio.iupui.edu:8080/miR2Disease/download/miRtar.txt
 :argument: DB_DESTINATION: saving location of the database
"""

# Imports
from SLKlib.SQLiteDBApi.sqlite_db_api import PsimiSQL
import requests, json, os, re
import urllib.parse

# Defining constants
SQL_SEED = '../../../../../SLKlib/SQLiteDBApi/network-db-seed.sql'
DATA_FILE = 'files/miRtar.txt'
DB_DESTINATION = '../../output/miR2Disease'


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


def is_well_formed_id(id):
    return re.search("^[/\\.\\w-]+$", id)


def main(logger):

    # Initiating the parser
    db_api = PsimiSQL(SQL_SEED)
    all_data = []
    pubmed_id_map = {}

    if os.path.exists('./pubmed_id_map_cache_for_miR2Disease.json'):
        with open("./pubmed_id_map_cache_for_miR2Disease.json") as cache:
            pubmed_id_map = json.load(cache)

    with open(DATA_FILE, encoding='ISO-8859-1') as data:
        # Skipping the header
        data.readline()
        data.readline()
        data.readline()
        node_names_to_id = {}
        lines = 0

        for line in data:
            columns = line.split('\t')
            if len(columns) > 1:
                lines += 1
                if lines % 50 == 0:
                    print("processed lines (miR2Disease): %d" % lines)

                columns[3] = columns[3].strip()

                if not is_well_formed_id(columns[0].strip()) or not is_well_formed_id(columns[1].strip()):
                    print("Warning: malformed ID, link skipped")
                    continue

                all_data.append(columns)

                if columns[3] not in pubmed_id_map:
                    search_term = columns[3].replace(".", ' ').replace(' and ', ' ').replace(' or ', ' ').replace("'", '').strip()
                    search_term = "%s[pdat] AND %s" % (columns[2].strip(), search_term)
                    search_term = urllib.parse.quote(search_term, safe='')
                    URL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&report=uilist&retmode=json&term=' + search_term
                    resp = requests.get(URL)
                    pubmed_id_list = json.loads(resp.text)['esearchresult']['idlist']

                    if pubmed_id_list and len(pubmed_id_list)==1:
                        pubmed_id_map[columns[3]] = pubmed_id_list[0]
                    else:
                        print("WARNING: pmid not found")
                        # print("         pubmed ID list: %s" % str(pubmed_id_list))
                        # print("         %s %s" % (columns[2], columns[3]))
                        # print("         " + URL)
                        pubmed_id_map[columns[3]] = None

        print("processed lines (miR2Disease): %d" % lines)

    print("saving output db")
    for columns in all_data:

        source_dict = insert_or_get_node_dict(columns[0], 'miRBase', 'taxid:9606', node_names_to_id, db_api)
        target_dict = insert_or_get_node_dict(columns[1], 'GeneCards', 'taxid:9606', node_names_to_id, db_api)

        # Getting files from the web with a custom URL
        pubmed_ids = ['18927107'] # mir2Disease publication
        if columns[3] in pubmed_id_map and pubmed_id_map[columns[3]]:
            pubmed_ids.append(str(pubmed_id_map[columns[3]]).strip())
        pubmed_ids = set(map(lambda x: ("pubmed:"+x).strip(), pubmed_ids))

        interaction_types = "is_directed:true|is_direct:true|MI:0571(mrna cleavage)"

        # Inserting edges
        edge_dict = {
                'publication_ids': "|".join(pubmed_ids),
                'layer': '5',
                'source_db': 'miR2Disease',
                'interaction_identifiers': None,
                'confidence_scores': None,
                'interaction_detection_method': None,
                'interaction_types': interaction_types,
                'first_author': None
            }

        db_api.insert_edge(source_dict, target_dict, edge_dict)

    # Saving the to a DB_TYPE.db files
    db_api.save_db_to_file(DB_DESTINATION)

    # save pubmed_id_map, so that we can re-use next time
    with open("./pubmed_id_map_cache_for_miR2Disease.json", 'w') as cache:
        json.dump(pubmed_id_map, cache, indent=4, sort_keys=True)


if __name__ == '__main__':
    print("Parsing database...")
    main(logger=None)
    print("Parsing database is completed. SQLite database is saved to: " + DB_DESTINATION)


