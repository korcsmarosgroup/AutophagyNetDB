"""
Parses the miRecords database
 :argument: DATA_FILE: validated target data files: http://c1.accurascience.com/miRecords/download_data.php?v=4
 :argument: DB_DESTINATION: saving location of database
 :argument: SPECIES_DICT: dictionary containing species' name and taxonomy ids

"""

# Imports
from ARNlib.SQLiteDBApi.sqlite_db_api import PsimiSQL

# Defining constants
SQL_SEED = '../../../../../ARNlib/SQLiteDBApi/network-db-seed.sql'
DATA_FILE = 'files/miRecords_version4.txt'
DB_DESTINATION = '../../output/miRecords'
SPECIES_DICT = {
    'Homo sapiens' : { 'tax_id': '9606', 'id_prefix':'hsa'},
    'Drosophila melanogaster' : { 'tax_id': '7227', 'id_prefix':'dme'},
    'Caenorhabditis elegans' : { 'tax_id': '6239', 'id_prefix':'cel'},
    'Danio rerio' : { 'tax_id': '7955', 'id_prefix':'dre'},
}


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
COLUMN    HEADER                                    EXAMPLE VALUE

0         Pubmed_id                                 15105502.0
1         Target gene_species_scientific            Homo sapiens
2         Target gene_name                          HOXC8
3         Target gene_Refseq_acc                    NM_022658
4         Target site_number                        1.0
5         miRNA_species                             Homo sapiens
6         miRNA_mature_ID                           hsa-miR-196a
7         miRNA_regulation                          overexpression
8         Reporter_target gene/region               luciferase
9         Reporter link element                     {target site}
10        Test_method_inter                         
11        Target gene mRNA_level                    {unchanged}{unchanged}
12        Original description                      {Segments from all three UTRs ... }{Segments from all ...  }
13        Mutation_target region                    
14        Post mutation_method                      
15        Original description_mutation_region      
16        Target site_position                      2017.0
17        miRNA_regulation_site                     
18        Reporter_target site                      
19        Reporter link element                     {target site}
20        Test_method_inter_site                    {activity assay}
21        Original description_inter_site           
22        Mutation_target site                      
23        Post mutation_method_site                 
24        Original description_mutation_site        
25        Additional note                           
                                      
'''


def main(logger):
    # Initiating the parser
    db_api = PsimiSQL(SQL_SEED)

    with open(DATA_FILE) as data:
        # Skipping the header
        data.readline()
        node_names_to_id = {}

        for line in data:
            columns = line.split('\t')

            if columns[1] in SPECIES_DICT and columns[5] in SPECIES_DICT:


                # there are a few kinds of malformed miRbase IDs, like:
                #   [has-let-7a3b]    -->   hsa-let-7a3b
                #   hsa-miR-34b*      -->   hsa-miR-34b
                #   miR-143           -->   <tax_id>-miR-143
                #
                rna_id = columns[6]
                rna_id = rna_id.replace("[","").replace("]","").replace("has-", "hsa-")

                if not rna_id.startswith("hsa-") and not rna_id.startswith("dme-") and not rna_id.startswith("cel-") and not rna_id.startswith("dre-"):
                    if rna_id.startswith("miR-"):
                        rna_id = SPECIES_DICT[columns[5]]['id_prefix'] + "-" + rna_id
                    else:
                        print("WARNING: skipping interaction due to malformed miRBase ID: " + rna_id)
                        continue
                rna_id = rna_id.replace("*","")

                source_dict = insert_or_get_node_dict(rna_id, 'miRBase','taxid:' + SPECIES_DICT[columns[5]]['tax_id'], node_names_to_id, db_api)
                target_dict = insert_or_get_node_dict(columns[3], 'RefSeq','taxid:' + SPECIES_DICT[columns[1]]['tax_id'], node_names_to_id, db_api)

                interaction_types = "is_directed:true|is_direct:true|MI:0571(mrna cleavage)"

                # pubmed id example: 15105502.0
                pubmed_id = columns[0].split('.')[0].strip()
                pubmed_ids = ['18996891']  # miRecords publication
                if len(pubmed_id) > 0:
                    pubmed_ids.append(pubmed_id)
                pubmed_ids = set(map(lambda x: 'pubmed:' + x, pubmed_ids))

                edge_dict = {
                        'publication_ids': "|".join(pubmed_ids),
                        'layer': '5',
                        'source_db': 'miRecords',
                        'interaction_identifiers': None,
                        'confidence_scores': None,
                        'interaction_detection_method': None,
                        'interaction_types': interaction_types,
                        'first_author': None
                    }

                db_api.insert_edge(source_dict, target_dict, edge_dict)

    # Saving the to a DB_TYPE.db files
    db_api.save_db_to_file(DB_DESTINATION)
    print("miRecords finished")


if __name__ == '__main__':
    print("Parsing database...")
    main(logger = None)
    print("Parsing database is completed. SQLite database is saved to: " + DB_DESTINATION)


