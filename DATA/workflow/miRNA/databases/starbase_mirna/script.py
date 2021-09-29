'''
    miRNA-mRNA interactions
    :argument: FILE_TO_TAXID: dictionary translating file names to taxids
               file_to_detmet: dictionary translating filenames to interaction detection method
'''

# Imports
from SLKlib.SQLiteDBApi.sqlite_db_api import PsimiSQL

# Defining constants
SQL_SEED = '../../../../../SLKlib/SQLiteDBApi/network-db-seed.sql'
DATA_FILE_LIST = ['files/starbase_v3_miRNAmRNA.txt',
                  'files/starbase_v3_miRNA_degradome_human.txt',
                  'files/starbase_v3_miRNA_degradome_worm.txt',
                  'files/starbase_v3_miRNA_valid.txt']
DB_DESTINATION = '../../output/starbase'

FILE_DICT = {}


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
    db_api = PsimiSQL(SQL_SEED)

    node_names_to_id = {}
    for file in DATA_FILE_LIST:
        print("processing data file: " + file)
        with open(file) as data:

            # Skipping the header
            data.readline()
            data.readline()
            data.readline()
            data.readline()

            metainfo = FILE_DICT[file.split("/")[-1]]

            for line in data:
                columns = line.split('\t')
                if len(columns) < 2:
                    continue

                rna_id = columns[metainfo['rna_id_column']].strip()
                if metainfo['rna_id_type'] == 'miRBase' and (rna_id.endswith("-3p") or rna_id.endswith("-5p")):
                    # during the rna mapping, we dont care about the 3p/5p postfix
                    rna_id = rna_id[:-3]
                if metainfo['rna_id_type'] == 'HGNC':
                    rna_id = rna_id.lower()

                # The wormBase IDs in the mapping DB contains only uppercase IDs
                gene_id = columns[metainfo['gene_id_column']].strip()
                if metainfo['gene_id_type'] == 'WormBase':
                    gene_id = gene_id.upper()

                source_dict = insert_or_get_node_dict(rna_id, metainfo['rna_id_type'], metainfo['tax_id'], node_names_to_id, db_api)
                target_dict = insert_or_get_node_dict(gene_id, metainfo['gene_id_type'], metainfo['tax_id'],
                                                      node_names_to_id, db_api)

                interaction_types = "is_directed:true|is_direct:true|MI:0571(mrna cleavage)"

                scores = []
                for score_definition in metainfo['scores']:
                    value = columns[score_definition['column']].strip()
                    score_name = score_definition['score_name'].strip()
                    scores.append("%s:%s" % (score_name, value))

                # Inserting edges
                edge_dict = {
                    'publication_ids': 'pubmed:24297251',  # StarBase v2.0 publication
                    'layer': '5',
                    'source_db': 'StarBase',
                    'interaction_identifiers': None,
                    'confidence_scores': "|".join(scores),
                    'interaction_detection_method': metainfo['detection_method'],
                    'interaction_types': interaction_types,
                    'first_author': None
                }

                db_api.insert_edge(source_dict, target_dict, edge_dict)

    # Saving the to a DB_TYPE.db files
    db_api.save_db_to_file(DB_DESTINATION)
    print("StarBase finished " + DB_DESTINATION)


if __name__ == '__main__':
    print("Parsing database...")
    main(logger=None)
    print("Parsing database is completed. SQLite database is saved to: " + DB_DESTINATION)
