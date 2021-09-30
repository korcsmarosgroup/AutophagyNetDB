"""
Parses MiRDeath database
 :argument: DATA_FILE: http://www.rna-world.org/mirdeathdb/data/miRDeathDB_all_data.txt
 :argument: DB_DESTINATION: saving location of database
"""


# Imports
from ARNlib.SQLiteDBApi.sqlite_db_api import PsimiSQL

# Defining constants
SQL_SEED = '../../../../../ARNlib/SQLiteDBApi/network-db-seed.sql'
DATA_FILE = 'files/miRDeathDB_all_data.txt'
DB_DESTINATION = '../../output/miRDeathDB'


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


'''

COLUMN  HEADER                                   EXAMPLE
        
0       miRNA_symbol                             miR-106b
1       miRBase_mature_ID                       "hsa-miR-106b-5p,hsa-miR-106b-3p"
2       miRBase_ID                               "MIMAT0000680,MIMAT0004672"
3       Gene_Symbol                              CDKN1A
4       Pathway                                  apoptosis
5       Action_Mode                              down
6       Organism                                 human
7       Tissue                                   gastric cancer
8       PMID                                     18328430
9       tax_id                                   9606
10      geneid                                   1026
11      Synonyms                                 CAP20|CDKN1|CIP1|MDA-6|P21|SDI1|WAF1|p21CIP1
12      Links                                    HGNC:1784|MIM:116899|Ensembl:ENSG00000124762|HPRD:00298|Vega:OTTHUMG00000014603|M2D:hsa-miR-106b
13      chromosome                               6
14      map_location                             6p21.2
15      Description                              "cyclin-dependent kinase inhibitor 1A (p21, Cip1)"
16      type_of_gene                             protein-coding
17      Full_name_from_nomenclature_authority    "cyclin-dependent kinase inhibitor 1A (p21, Cip1)"
18      Other_designations                       CDK-interacting protein 1|CDK-interaction protein 1|DNA synthesis inhibitor|cyclin-dependent kinase inhibitor 1|melanoma differentiation associated protein 6|wild-type p53-activated fragment 1


'''


def main(logger):
    # Initiating the parser
    db_api = PsimiSQL(SQL_SEED)

    with open(DATA_FILE) as data:
        # Skipping the header
        data.readline()
        node_names_to_id = {}

        for line in data:
            columns = line.strip().split('\t')
            if len(columns) != 1:
                if columns[9] == '9606' or columns[9] == '7227':

                    mirbase_id = ("hsa-"+columns[0]) if columns[9] == '9606' else ("dme-"+columns[0])
                    # there are two malformed ID in the database:
                    # hsa-miR-149*       -->  hsa-miR-149
                    # hsa-"miR-34a,b,c"  -->  hsa-"miR-34"
                    mirbase_id = mirbase_id.replace('*', '').replace('\"', '').replace('a,b,c', '')

                    source_dict = insert_or_get_node_dict(mirbase_id, 'miRBase', 'taxid:' + columns[9], node_names_to_id, db_api)
                    target_dict = insert_or_get_node_dict(columns[3], 'GeneCards', 'taxid:' + columns[9], node_names_to_id, db_api)

                    interaction_types = "is_directed:true|is_direct:true|MI:0571(mrna cleavage)"

                    pubmed_ids = ['22743998']  # miRDeathDB publication
                    if len(columns[8].strip()) > 0:
                        pubmed_ids.append(columns[8].strip())
                    pubmed_ids = set(map(lambda x: 'pubmed:' + x, pubmed_ids))

                    # Inserting edges
                    edge_dict = {
                            'publication_ids': "|".join(pubmed_ids),
                            'layer': '5',
                            'source_db': 'miRDeathDB',
                            'interaction_identifiers': None,
                            'confidence_scores': None,
                            'interaction_detection_method': None,
                            'interaction_types': interaction_types,
                            'first_author': None
                        }

                    db_api.insert_edge(source_dict, target_dict, edge_dict)

    # Saving the to a DB_TYPE.db files
    db_api.save_db_to_file(DB_DESTINATION)
    print("miRDeathDB finished")


if __name__ == '__main__':
    print("Parsing database...")
    main(logger = None)
    print("Parsing database is completed. SQLite database is saved to: " + DB_DESTINATION)


