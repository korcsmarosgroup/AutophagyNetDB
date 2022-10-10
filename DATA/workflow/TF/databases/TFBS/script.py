
from ARNlib.SQLiteDBApi.sqlite_db_api import PsimiSQL

# Defining constants
SQL_SEED = '../../../../../ARNlib/SQLiteDBApi/network-db-seed.sql'
DB_TYPE = 'PSSMprediction'
EXPORT_DB_LOCATION = '../../output/TFBS'
DATA_FILE = 'files/HITSallbodyresults.txt'


def insert_or_get_node_dict(name, taxid, node_names_to_id, db_api):

    node_name = name.split(";")
    alt_node_names = []
    new_node_name = []
    for n in node_name:
        if n.startswith("ENST"):
            new_n = f'Ensembl:{n}'
            if len(new_node_name) == 0:
                new_node_name.append(new_n)
            else:
                alt_node_names.append(new_n)
        else:
            new_n = f'HGNC:{n}'
            if "::" in n:
                n = n.replace("::", "/")
                new_n = f'HGNC:{n}'
            elif "(" in n:
                n = n.split("(")[0]
                new_n = f'HGNC:{n}'
            new_node_name.append(new_n)

    node_dict = {
        "name": "|".join(new_node_name),
        "tax_id": taxid,
        "alt_accession": "|".join(alt_node_names),
        'pathways': None,
        "aliases": None,
        "topology": "Transcription factor"
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
            columns = line.split('\t')
            lnc_name = columns[0].split('|')[3]
            taxid = 'taxid:9606'

            # Creating the node dicts, if the node is already in the db assigning that to the node dict
            source_dict = insert_or_get_node_dict(lnc_name, taxid, node_names_to_id, db_api)
            target_dict = insert_or_get_node_dict(columns[2], taxid, node_names_to_id, db_api)

            interaction_types = "%s|is_directed:%s|is_direct:%s" \
                                % ('MI:2247(trascriptional regulation)', 'true', 'false')
            edge_dict = {
                'publication_ids': 'pubmed:28591841|pubmed:29140473',
                'layer': '6',
                'source_db': 'PSSMprediction',
                'interaction_identifiers': None,
                'confidence_scores': None,
                'interaction_detection_method': 'MI:1176(sequence based prediction of gene regulatory region binding sites)',
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



