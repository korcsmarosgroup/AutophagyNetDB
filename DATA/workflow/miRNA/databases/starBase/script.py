'''
    miRNA-mRNA interactions
    :argument: FILE_TO_TAXID: dictionary translating file names to taxids
               file_to_detmet: dictionary translating filenames to interaction detection method
'''

# Imports
from ARNlib.SQLiteDBApi.sqlite_db_api import PsimiSQL

# Defining constants
SQL_SEED = '../../../../../ARNlib/SQLiteDBApi/network-db-seed.sql'
DATA_FILE_LIST = ['files/starbase_v3_miRNAmRNA.txt',
                  'files/starbase_v3_miRNA_degradome_human.txt',
                  'files/starbase_v3_miRNA_degradome_worm.txt',
                  'files/starbase_v3_miRNA_valid.txt']
DB_DESTINATION = '../../output/starbase'

FILE_TO_TAXID = {
    'files/starbase_v3_miRNAmRNA.txt': "taxid:9606",
    'files/starbase_v3_miRNA_degradome_human.txt': 'taxid:9606',
    'files/starbase_v3_miRNA_degradome_worm.txt': 'taxid:6239',
    'files/starbase_v3_miRNA_valid.txt': 'taxid:9606'
    }
file_to_detmet = {
    'files/starbase_v3_miRNAmRNA.txt': "MI:1110(predicted interaction)",
    'files/starbase_v3_miRNA_degradome_human.txt': 'MI:0045(experimental interaction detection)',
    'files/starbase_v3_miRNA_degradome_worm.txt': 'MI:0045(experimental interaction detection)',
    'files/starbase_v3_miRNA_valid.txt': 'MI:2321(high throughput sequencing)'
}

def get_node_mirna(mirna_name, taxid, psi_mi_to_sql_object):
    """
    This function sets up a node dict and returns it. If the node is already in the SQLite database it fetches that node from the db, so it won't be inserted multiple times.
    """

    # Testing if the node is already in the database
    node_dict = psi_mi_to_sql_object.get_node(mirna_name, node_tax_id=taxid)

    if not node_dict:
        node_dict = {
            "name" : 'mirbase:' + mirna_name,
            "tax_id": taxid,
            "alt_accession": None,
            'pathways': None,
            "aliases": None,
            "topology": None
        }

    return node_dict


def get_node_gene(gene_name, taxid, psi_mi_to_sql_object):
    """
    This function sets up a node dict and returns it. If the node is already in the SQLite database it fetches that node from the db, so it won't be inserted multiple times.
    """

    # Testing if the node is already in the database
    node_dict = psi_mi_to_sql_object.get_node(gene_name, node_tax_id=taxid)

    if not node_dict:
        node_dict = {
            "name": 'HGNC:' + gene_name,
            "tax_id": taxid,
            "alt_accession": None,
            'pathways': None,
            "aliases": None,
            "topology": None
        }

    return node_dict


def main(logger):
    # Declaring variables and constants
    inserted_nodes = {}

    # Initiating the parser
    db_api = PsimiSQL(SQL_SEED)
    for file in DATA_FILE_LIST:
        print(file)
        with open(file) as data:
            # Skipping the header
            data.readline()
            data.readline()
            data.readline()
            data.readline()
            for line in data:
                columns = line.split('\t')
                taxid = FILE_TO_TAXID[file]
                if len(columns) != 1:
                    if file == 'miRNA/databases/starBase/files/starbase_v3_miRNAmRNA.txt':
                        mirna_name = columns[1]
                        gene_name = columns[3]
                    elif file == 'miRNA/databases/starBase/files/starbase_v3_degradome_human.txt'\
                            or file == 'miRNA/databases/starBase/files/starbase_v3_degradome_worm.txt':
                        mirna_name = columns[1]
                        gene_name = columns[2]
                    elif file == 'miRNA/databases/starBase/files/starbase_v3_miRNA_valid.txt':
                        mirna_name = columns[1]
                        gene_name = columns[4]
                    else:
                        mirna_name = None
                        gene_name = None
                    # Creating the node dicts, if the node is already in the db assigning that to the node dict
                    source_dict = get_node_mirna(mirna_name, taxid, db_api)
                    target_dict = get_node_mirna(gene_name, taxid, db_api)

                    # Nodes are inserted to the db if they are not in it yet
                    if not 'id' in source_dict:
                        db_api.insert_node(source_dict)

                    if not 'id' in target_dict:
                        db_api.insert_node(target_dict)

                    interaction_types = "effect:%s|is_directed:%s|is_direct:%s" \
                                        % ('MI:0256(rna interference)', 'directed', 'unknown')

                    # Inserting edges
                    edge_dict = {
                        'publication_ids': 'pubmed:24297251',
                        'layer': '5',
                        'source_db': 'starbase',
                        'interaction_identifiers': None,
                        'confidence_scores': None,
                        'interaction_detection_method': file_to_detmet[file],
                        'interaction_types': interaction_types,
                        'first_author': None
                    }

                    db_api.insert_edge(source_dict, target_dict, edge_dict)

    # Saving the to a DB_TYPE.db files
    db_api.save_db_to_file(DB_DESTINATION)


if __name__ == '__main__':
    print("Parsing database...")
    main(logger = None)
    print("Parsing database is completed. SQLite database is saved to: " + DB_DESTINATION)


