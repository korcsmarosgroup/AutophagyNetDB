'''
    miRNA-lncRNA interactions
    :argument: DATA_FILE: //bioinfo.life.hust.edu.cn/lncRInter/browse?species=&class=&level=RNA-RNA
'''

# Imports
from ARNlib.SQLiteDBApi.sqlite_db_api import PsimiSQL

# Defining constants
SQL_SEED = '../../../../../ARNlib/SQLiteDBApi/network-db-seed.sql'
DATA_FILE = 'files/result_all.txt'
DB_DESTINATION = '../../output/lncRInter'

ORGANISM_TO_TAXID = {
    "homo sapiens": "taxid:9606",
    "drosophila melanogaster": "taxid:7227",
    "caenorhabditis elegans": "taxid:6239",
    "danio rerio" : "taxid:7955"
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

def main(logger):

    # Initiating the parser
    db_api = PsimiSQL(SQL_SEED)
    node_names_to_id = {}

    with open(DATA_FILE) as data:
        # Skipping the header
        data.readline()

        for line in data:
            columns = line.split('\t')
            if columns[4].strip().lower() in ORGANISM_TO_TAXID:
                source_dict = insert_or_get_node_dict(columns[0], "HGNC", ORGANISM_TO_TAXID[columns[4].strip().lower()], node_names_to_id, db_api)
                target_dict = insert_or_get_node_dict(columns[1], "HGNC", ORGANISM_TO_TAXID[columns[4].strip().lower()], node_names_to_id, db_api)

                interaction_types = "is_directed:true|is_direct:true|MI:0407(direct interaction)"

                detmap = {
                    'pull-down assay': 'MI:0096(pull down)',
                    'qPCR, Western blot, RIP.': 'MI:1195(quantitative pcr)|MI:0113(western blot)|MI:1017(rna immunoprecipitation)',
                    'qRT-PCR, RNAi': 'MI:1196(quantitative reverse transcription pcr)',
                    'in vitro': 'MI:0045(experimental interaction detection)',
                    'In vitro': 'MI:0045(experimental interaction detection)',
                    'Luciferase reporter assay, Pulldown assay ': 'MI:2285(miRNA interference luciferase reporter assay)|MI:0096(pull down)',
                    'luciferase reporter assays and pull-down assay': 'MI:2285(miRNA interference luciferase reporter assay)|MI:0096(pull down)',
                    'RIP': 'MI:1017(rna immunoprecipitation)',
                    'luciferase reporter constructs': 'MI:2285(miRNA interference luciferase reporter assay)',
                    'in vitro and vivo': 'MI:0045(experimental interaction detection)',
                    'dual luciferase reporter assay?': 'MI:2285(miRNA interference luciferase reporter assay)',
                    'dual luciferase reporter assays': 'MI:2285(miRNA interference luciferase reporter assay)',
                    'qPCR, RNAi etc.': 'MI:1195(quantitative pcr)',
                    'ISH and Luciferase Assay': 'MI:2285(miRNA interference luciferase reporter assay)',
                    'In vitro RNA/dsDNA binding assay utilizing biotin tagged RNA oligos as bait': 'MI:0045(experimental interaction detection)',
                    'In vitro RNA/dsDNA binding assay': 'MI:0045(experimental interaction detection)',
                    'biotin-avidin pull-down system': 'MI:0096(pull down)',
                    'microRNA crosslinking and immunoprecipitation (miR-CLIP)': 'MI:2191(clip)',
                    'Luciferase reporter assay and qPCR': 'MI:2285(miRNA interference luciferase reporter assay)',
                    'RNA immunoprecipitation;luciferase reporter assays': 'MI:2285(miRNA interference luciferase reporter assay)|MI:1017(rna immunoprecipitation)',
                    'in vivo': 'MI:0045(experimental interaction detection)',
                    'luciferase reporter assays': 'MI:2285(miRNA interference luciferase reporter assay)',
                    'RNA immunoprecipitation and luciferase reporter assays': 'MI:2285(miRNA interference luciferase reporter assay)|MI:1017(rna immunoprecipitation)',
                    'EMSA': 'MI:0413(electrophoretic mobility shift assay)',
                    'luciferase reporter assay': 'MI:2285(miRNA interference luciferase reporter assay)',
                    'Luciferase assays': 'MI:2285(miRNA interference luciferase reporter assay)',
                    '-': 'MI:0045(experimental interaction detection)',
                    'RNA immunoprecipitation': 'MI:1017(rna immunoprecipitation)',
                    'RIP, Biotin-RNA Pull-Down Assay,qRT-PCR,EMSA': 'MI:1017(rna immunoprecipitation)|MI:0096(pull down)|MI:0413(electrophoretic mobility shift assay)|MI:1196(quantitative reverse transcription pcr)',
                    'Luciferase reporter assay, RIP assay and RNA pull-down assay': 'MI:2285(miRNA interference luciferase reporter assay)|MI:1017(rna immunoprecipitation)|MI:0096(pull down)',
                    'qPCR, Western blot and RNAi': 'MI:1195(quantitative pcr)|MI:0113(western blot)',
                    'luciferase  reporter  assay': 'MI:2285(miRNA interference luciferase reporter assay)',
                    'CLIP': 'MI:2191(clip)',
                    'RIP and ChIP assay ': 'MI:1017(rna immunoprecipitation)',
                    'in vitro or vivo': 'MI:0045(experimental interaction detection)',
                    'RNA pull-down assay': 'MI:0096(pull down)',
                    'immunoprecipitation (RIP) assay and RNA pull-down assay': 'MI:1017(rna immunoprecipitation)|MI:0096(pull down)',
                    'luciferase reporter': 'MI:2285(miRNA interference luciferase reporter assay)',
                    'in vitro and in vivo': 'MI:0045(experimental interaction detection)',
                    'in viro': 'MI:0045(experimental interaction detection)',
                    'co-RNA-FISH assays': 'MI:0045(experimental interaction detection)',
                    'luciferase reporter ': 'MI:2285(miRNA interference luciferase reporter assay)',
                    'microarray, qPCR': 'MI:1195(quantitative pcr)',
                    'In vitro and in vivo': 'MI:0045(experimental interaction detection)',
                    'Luciferase reporter assays': 'MI:2285(miRNA interference luciferase reporter assay)',
                    'RIP and Luciferase assays': 'MI:2285(miRNA interference luciferase reporter assay)|MI:1017(rna immunoprecipitation)',
                    'RNA-FISH': 'MI:0045(experimental interaction detection)',
                    'RNA FISH': 'MI:0045(experimental interaction detection)',
                    'FISH, Allele-specific RT-PCR': 'MI:1196(quantitative reverse transcription pcr)',
                    'RIP and RNA pull-down': 'MI:1017(rna immunoprecipitation)',
                    'RIP and ChIP assay': 'MI:0019(coimmunoprecipitation)'
                }


                detmethod = None
                if columns[8].strip() in detmap:
                    detmethod = detmap[columns[8].strip()]
                else:
                    print("WARNING: unknown detection method: " + columns[8].strip())

                pubmed_ids = ['28529080']  # lncRInter publication
                pubmed_id = columns[9].strip()
                if len(pubmed_id) > 0:
                    pubmed_ids.append(pubmed_id)
                pubmed_ids = set(map(lambda x: 'pubmed:' + x, pubmed_ids))



                # Inserting edges
                edge_dict = {
                            'publication_ids': "|".join(pubmed_ids),
                            'layer': '7',
                            'source_db': 'lncRInter',
                            'interaction_identifiers': None,
                            'confidence_scores': None,
                            'interaction_detection_method': detmethod,
                            'interaction_types': interaction_types,
                            'first_author': None
                        }

                db_api.insert_edge(source_dict, target_dict, edge_dict)

    # Saving the to a DB_TYPE.db files
    db_api.save_db_to_file(DB_DESTINATION)
    print("lncRInter finished")


if __name__ == '__main__':
    print("Parsing database...")
    main(logger = None)
    print("Parsing database is completed. SQLite database is saved to: " + DB_DESTINATION)


