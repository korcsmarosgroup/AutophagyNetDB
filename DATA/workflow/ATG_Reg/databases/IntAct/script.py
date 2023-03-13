"""
Parses IntAct data files
    :argument: DATA_FILE: ftp://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact-micluster.txt
    :argument: DB_TYPE: name of the database
    :argument: EXPORT_DB_LOCATION: saving location
"""

# Imports
from ARNlib.SQLiteDBApi.sqlite_db_api import PsimiSQL
import re

# Defining constants
SQL_SEED = '../../../../../ARNlib/SQLiteDBApi/network-db-seed.sql'
DB_TYPE = 'intact'
EXPORT_DB_LOCATION = '../../output/IntAct'
DATA_FILE = 'files/intact-micluster.txt'


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


def is_well_formed_uniprot_id(alt_id):
    #uniprotkb:P53999

    if not alt_id.strip().lower().startswith('uniprotkb'):
        return False

    id = alt_id[10:].strip()
    for bad_character in (":", "(", ")", " ", "/"):
        if bad_character in id:
            return False

    return True


'''
COLUMN       HEADER                                  EXAMPLE
0            ID(s) interactor A                      intact:EBI-998260
1            ID(s) interactor B                      intact:EBI-7333987
2            Alt. ID(s) interactor A                 uniprotkb:P53999
3            Alt. ID(s) interactor B                 uniprotkb:P04326
4            Alias(es) interactor A                  psi-mi:tcp4_human|psi-mi:SUB1|uniprotkb:Q96L29|uniprotkb:SUB1 homolog|uniprotkb:Positive cofactor 4|uniprotkb:p14|uniprotkb:SUB1|uniprotkb:PC4|uniprotkb:RPO2TC1
5            Alias(es) interactor B                  psi-mi:tat_hv112|psi-mi:tat|uniprotkb:tat|uniprotkb:Transactivating regulatory protein|intact:MINT-138666
6            Interaction detection method(s)         psi-mi:"MI:0096"(pull down)|psi-mi:"MI:0018"(two hybrid)|psi-mi:"MI:0019"(coimmunoprecipitation)
7            Publication 1st author(s)               Holloway et al. (2000)
8            Publication Identifier(s)               pubmed:10887206|mint:MINT-5212759
9            Taxid interactor A                      taxid:9606(human)
10           Taxid interactor B                      taxid:11679(hv112)
11           Interaction type(s)                     psi-mi:"MI:0407"(direct interaction)|psi-mi:"MI:0915"(physical association)|psi-mi:"MI:0915"(physical association)
12           Source database(s)                      psi-mi:"MI:0471"(MINT)
13           Interaction identifier(s)               intact:EBI-7334047|intact:EBI-7333982|intact:EBI-7334070
14           Confidence value(s)                     intact-miscore:0.59257585
'''


def main(logger):
    # Initiating the parser
    db_api = PsimiSQL(SQL_SEED)


    # Parsing files
    with open(DATA_FILE) as data:
        data.readline()
        missing_alt_id = 0
        edges_in_known_taxes = 0
        node_names_to_id = {}
        lines = 0

        for line in data:
            columns = line.split('\t')
            if len(columns) < 2:
                continue
            lines += 1
            if lines % 50000 == 0:
                print("processed lines (IntAct): %d" % lines)

            # tax id like: taxid:9606(human), taxid:83333(ecoli), taxid:-1(in vitro)
            tax_id_a = columns[9][:10]
            tax_id_b = columns[10][:10]

            if tax_id_a not in ('taxid:9606') or \
               tax_id_b not in ('taxid:9606'):
                continue

            edges_in_known_taxes += 1

            if is_well_formed_uniprot_id(columns[2]) and is_well_formed_uniprot_id(columns[3]):
                
                uniprot_id_a = columns[2][10:].strip()
                uniprot_id_b = columns[3][10:].strip()
                
                # Creating the node dicts, if the node is already in the db assigning that to the node dict
                source_dict = insert_or_get_node_dict(uniprot_id_a, 'Uniprot', tax_id_a, node_names_to_id, db_api)
                target_dict = insert_or_get_node_dict(uniprot_id_b, 'Uniprot', tax_id_b, node_names_to_id, db_api)

                # interaction detection methods: psi-mi:"MI:0096"(pull down)|psi-mi:"MI:0018"(two hybrid)
                detection_methods = columns[6].split("|")
                detection_methods = filter(lambda x: x.strip().lower().startswith('psi-mi:'), detection_methods)
                detection_methods = map(lambda x: x.strip()[7:].replace("\"","").strip(), detection_methods)
                detection_methods = set(detection_methods)

                # pubmed ids: pubmed:10887206|mint:MINT-5212759
                pubmed_ids = columns[8].split("|")
                pubmed_ids = filter(lambda x: x.strip().lower().startswith('pubmed:'), pubmed_ids)
                pubmed_ids = map(lambda x: x.strip()[7:], pubmed_ids)
                pubmed_ids = filter(lambda x: re.search("^\\d+$", x), pubmed_ids)
                pubmed_ids = set(pubmed_ids)
                pubmed_ids.add("24234451")  # intact publication
                pubmed_ids = map(lambda x: "pubmed:" + x, pubmed_ids)

                # interaction type: psi-mi:"MI:0407"(direct interaction)|psi-mi:"MI:0915"(physical association)
                interaction_type_terms = columns[11].split("|")
                interaction_type_terms = filter(lambda x: x.strip().lower().startswith('psi-mi:'), interaction_type_terms)
                interaction_type_terms = map(lambda x: x.strip()[7:].replace("\"","").strip(), interaction_type_terms)
                interaction_type_terms = set(interaction_type_terms)

                # we remove 'MI:0407(direct interaction)' term, as it is redundant with the is_direct attribute
                interaction_type_terms.discard("MI:0407(direct interaction)")

                interaction_type = "is_directed:false|is_direct:true"
                if len(interaction_type_terms)>0:
                    interaction_type += "|"+"|".join(interaction_type_terms)

                # interaction score examples in the IntAct input file:
                # - intact-miscore:0.558037
                # - author score:low
                # - author score:Retest score=6; Class=Core; confidence score set1/set2 =2
                # - author score:"Socio-affinity score: 6.11118"
                # - author-confidence:Z-score = 17.60
                # - replication-based confidence:4
                # we don't keep the author-type scores, as those are a mess and also contains several non-numeric scores
                confidence_scores = columns[14].split("|")
                confidence_scores = map(lambda x: x.strip(), confidence_scores)
                confidence_scores = filter(lambda x: not x.startswith("author score:") and not x.startswith("author-confidence:"), confidence_scores)
                confidence_scores = map(lambda x: x.replace("intact-miscore", "intact miscore"), confidence_scores)
                confidence_scores = map(lambda x: x if x.lower().startswith("intact") else "intact "+x, confidence_scores)
                confidence_scores = set(confidence_scores)

                edge_dict = {
                    'publication_ids': "|".join(pubmed_ids),
                    'layer': '2',
                    'source_db': "IntAct",
                    'interaction_identifiers': None,
                    'confidence_scores':  "|".join(confidence_scores),
                    'interaction_detection_method': "|".join(detection_methods),
                    'interaction_types': interaction_type,
                    'first_author': None
                }

                db_api.insert_edge(source_dict, target_dict, edge_dict)
            else:
                missing_alt_id += 1

        print("processed lines (IntAct): %d" % lines)
        print("number of links in the known species: %d" % edges_in_known_taxes)
        print("links with missing uniprot ID: %d" % missing_alt_id)

    # Saving the to a DB_TYPE.db files
    db_api.save_db_to_file(EXPORT_DB_LOCATION)


if __name__ == '__main__':
    print("Parsing database...")
    main(logger=None)
    print("Parsing database is completed. SQLite database is saved to: " + EXPORT_DB_LOCATION)



