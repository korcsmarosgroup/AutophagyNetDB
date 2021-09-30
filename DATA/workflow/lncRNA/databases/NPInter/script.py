"""
Parses NPInter database
 :argument: DATA_FILE: data files (linux version) http://www.bioinfo.org/NPInter/datadownload/interaction_NPInter[v3.0].txt.tar.gz
 :argument: DB_DESTINATION: saving location of database
 :argument: SPECIES_DICT: dictionary containing species' name and taxonomy ids

"""

# Imports
from ARNlib.SQLiteDBApi.sqlite_db_api import PsimiSQL
import re, sys

# Defining constants
SQL_SEED = '../../../../../ARNlib/SQLiteDBApi/network-db-seed.sql'
DATA_FILE = 'files/interaction_NPInter[v3.0].txt'
DB_DESTINATION = '../../output/NPInter'

SPECIES_DICT = {
    'Homo sapiens' : { 'tax_id': '9606', 'id_prefix':'hsa'},
    'Drosophila melanogaster' : { 'tax_id': '7227', 'id_prefix':'dme'},
    'Caenorhabditis elegans' : { 'tax_id': '6239', 'id_prefix':'cel'},
    'Danio rerio' : { 'tax_id': '7955', 'id_prefix':'dre'},
}




def insert_or_get_node_dict(id, taxid, node_names_to_id, db_api):
    node_dict = {
        "name": id.strip(),
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


def fix_id(id_type, id, name, tax_dict, id_type_map):
    new_id = id

    if id_type not in id_type_map:
        print("ERROR! unknown ID type: " + id_type)
        sys.exit(1)

    id_type_info = id_type_map[id_type]

    if id_type_info['use_name'] or id.lower() == "null":
        new_id = name
    
    new_id = new_id.replace(" ", "")
    new_id = new_id.replace("*", "")

    if id_type_info['id_type'] == 'miRBase':
        if new_id.lower().startswith('mir-'):
            new_id = "%s-%s" % (tax_dict['id_prefix'], new_id)
        if new_id.endswith("-3p") or new_id.endswith("-5p"):
            new_id = new_id[:-3]

    if re.search("^[/\\.\\w-]+$", new_id):

        return '%s:%s' % (id_type_info['id_type'], new_id)
    else:
        print("WARNING: malformed ID value: " + new_id)
        return None


def main(logger):
    db_api = PsimiSQL(SQL_SEED)

    node_names_to_id = {}
    with open(DATA_FILE, encoding='ISO-8859-1') as data:

        lines = 0

        for line in data:
            lines += 1
            if lines % 50000 == 0:
                print("processed lines: %d" % lines)

            columns = line.strip().split('\t')
            if columns[-1] == 'RNA-RNA':
                if columns[12] in SPECIES_DICT:

                    tax_id = 'taxid:' + SPECIES_DICT[columns[12]]['tax_id']

                    id_type_map = {
                        'NONCODE': { 'id_type':'NONCODE', 'use_name':False},
                        'miRBase': { 'id_type':'miRBase', 'use_name':True},
                        'UniProt': { 'id_type':'Uniprot', 'use_name':False},
                        'UniGene': { 'id_type':'GeneCards', 'use_name':True},
                        'RefSeq': { 'id_type':'RefSeq', 'use_name':False},
                    }

                    source_id = fix_id(columns[2].strip(), columns[3].strip(), columns[4].strip(), SPECIES_DICT[columns[12]], id_type_map)
                    target_id = fix_id(columns[6].strip(), columns[7].strip(), columns[4].strip(), SPECIES_DICT[columns[12]], id_type_map)

                    if not source_id or not target_id:
                        continue

                    source_dict = insert_or_get_node_dict(source_id, tax_id, node_names_to_id, db_api)
                    target_dict = insert_or_get_node_dict(target_id, tax_id, node_names_to_id, db_api)

                    interaction_types = "MI:0407(direct interaction)|is_directed:true|is_direct:true"

                    pubmed_ids = ['27087310']  # NPInter publication
                    pubmed_id = columns[11].strip()
                    if len(pubmed_id) > 0 and re.search("^\\d+$", pubmed_id):
                        pubmed_ids.append(pubmed_id)
                    pubmed_ids = set(map(lambda x: 'pubmed:' + x, pubmed_ids))

                    edge_dict = {
                        'publication_ids': "|".join(pubmed_ids),
                        'layer': '7',
                        'source_db': 'NPInter',
                        'interaction_identifiers': None,
                        'confidence_scores': None,
                        'interaction_detection_method': None,
                        'interaction_types': interaction_types,
                        'first_author': None
                    }

                    db_api.insert_edge(source_dict, target_dict, edge_dict)
        print("processed lines: %d" % lines)

    # Saving the to a DB_TYPE.db files
    db_api.save_db_to_file(DB_DESTINATION)
    print("NPInter finished")


if __name__ == '__main__':
    print("Parsing database...")
    main(logger = None)
    print("Parsing database is completed. SQLite database is saved to: " + DB_DESTINATION)

