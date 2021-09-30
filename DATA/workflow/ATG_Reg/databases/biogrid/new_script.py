"""
 New parsing
"""

# Imports
import re
from ARNlib.SQLiteDBApi.sqlite_db_api import PsimiSQL

# Defining constants
SQL_SEED = '../../../../../../../ARNlib/SQLiteDBApi/network-db-seed.sql'
path = '../../../ATG_Reg/databases/biogrid/files/BIOGRID-SYSTEM-3.4.145.mitab/'
RAW_FILE_LIST = [path + 'BIOGRID-SYSTEM-Positive_Genetic-3.4.145.mitab.txt']
DB_TYPE = 'biogrid'
EXPORT_DB_LOCATION = '../../output/'


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

# COLUMN                  HEADER                    EXAMPLE
0                     ID Interactor A	                entrez gene/locuslink:375
1                     ID Interactor B	                entrez gene/locuslink:23163
2                     Alt IDs Interactor A	            biogrid:106870|entrez gene/locuslink:ARF1
3                     Alt IDs Interactor B	            biogrid:116775|entrez gene/locuslink:GGA3
4                     Aliases Interactor A	            -
5                     Aliases Interactor B	            -
6                     Interaction Detection Method	    psi-mi:"MI:0018"(two hybrid)
7                     Publication 1st Author	        "DellAngelica EC (2000)"
8                     Publication Identifiers	        pubmed:10747089
9                     Taxid Interactor A	            taxid:9606
10                    Taxid Interactor B	            taxid:9606
11                    Interaction Types	                psi-mi:"MI:0407"(direct interaction)
12                    Source Database	                psi-mi:"MI:0463"(biogrid)
13                    Interaction Identifiers	        biogrid:586
14                    Confidence Values                 -

'''


def main(logger):
    # Initiating the parser
    db_api = PsimiSQL(SQL_SEED)

    for file in RAW_FILE_LIST:
        with open(file) as data:
            # Skipping the header
            data.readline()
            node_names_to_id = {}
            lines = 0

            for line in data:
                if line.strip() != '' and not line.strip().startswith("#"):
                    lines += 1
                    if lines % 50000 == 0:
                        print("processed lines (biogrid): %d" % lines)

                    columns = line.split('\t')

                    # if a cell contains only '-', we replace it with empty string
                    columns = list(map(lambda x: "" if x.strip() == "-" else x, columns))

                    tax_id_a = columns[9].strip().lower()
                    tax_id_b = columns[10].strip().lower()

                    if tax_id_a in ('taxid:9606', 'taxid:7227', 'taxid:6239', 'taxid:7955') and \
                       tax_id_b in ('taxid:9606', 'taxid:7227', 'taxid:6239', 'taxid:7955'):

                        biogrid_ids_a = filter(lambda x: x.strip().lower().startswith("biogrid:"), columns[2].split("|"))
                        biogrid_id_a = list(biogrid_ids_a)[0][8:]

                        biogrid_ids_b = filter(lambda x: x.strip().lower().startswith("biogrid:"), columns[3].split("|"))
                        biogrid_id_b = list(biogrid_ids_b)[0][8:]

                        # Creating the node dicts, if the node is already in the db assigning that to the node dict
                        source_dict = insert_or_get_node_dict(biogrid_id_a, 'BioGrid', tax_id_a, node_names_to_id, db_api)
                        target_dict = insert_or_get_node_dict(biogrid_id_b, 'BioGrid', tax_id_b, node_names_to_id, db_api)

                        # interaction types in biogrid:
                        # direct:
                        #    - psi-mi:"MI:0407"(direct interaction)
                        #    - psi-mi:"MI:0915"(physical association)
                        #    - psi-mi:"MI:0914"(association)
                        # indirect:
                        #    - psi-mi:"MI:0799"(additive genetic interaction defined by inequality)
                        #    - psi-mi:"MI:0403"(colocalization)
                        #    - psi-mi:"MI:0796"(suppressive genetic interaction defined by inequality)
                        #    - psi-mi:"MI:0794"(synthetic genetic interaction defined by inequality)
                        mi_number = columns[11][11:15]
                        if mi_number not in ("0407", "0915", "0914", "0799", "0403", "0796", "0794"):
                            print("warning: unknown interaction type: " + columns[11])
                        is_direct = True
                        if mi_number in ("0799", "0403", "0796", "0794"):
                            is_direct = False

                        # we add the MI term to the interaction_types
                        # but we skip MI:0407(direct interaction) -> this info is already presented in the is_direct attribute
                        output_mi_string = "|" + columns[11].replace("psi-mi:", "").replace("\"", "")
                        if "MI:0407" in output_mi_string:
                            output_mi_string = ""

                        interaction_types = "is_directed:false|is_direct:%s%s" % (str(is_direct).lower(), output_mi_string)

                        # Interaction detection methods: psi-mi:"MI:0018"(two hybrid)
                        detection_methods = columns[6].split("|")
                        detection_methods = map(lambda x: x[7:] if x.lower().startswith('psi-mi') else x, detection_methods)
                        detection_methods = map(lambda x: x.replace("\"", ""), detection_methods)

                        # pubmed ids: pubmed:10747089
                        pubmed_ids = columns[8].split("|")
                        pubmed_ids = map(lambda x: x[7:] if x.lower().startswith('pubmed') else x, pubmed_ids)
                        pubmed_ids = filter(lambda x: re.search("^\\d+$", x), pubmed_ids)
                        pubmed_ids = set(pubmed_ids)
                        pubmed_ids.add("30476227")  # latest biogrid publication
                        pubmed_ids = map(lambda x: "pubmed:"+x, pubmed_ids)

                        edge_dict = {
                            'publication_ids': "|".join(pubmed_ids),
                            'layer': '3',
                            'source_db': 'TheBiogrid',
                            'interaction_identifiers': None,
                            'confidence_scores': None,
                            'interaction_detection_method': "|".join(detection_methods),
                            'interaction_types': interaction_types,
                            'first_author': None
                        }

                        db_api.insert_edge(source_dict, target_dict, edge_dict)

            print("processed lines: %d" % lines)

    # Saving the to a DB_TYPE.db files
    db_api.save_db_to_file(EXPORT_DB_LOCATION)


if __name__ == '__main__':
    print("Parsing database...")
    main(logger=None)
    print("Parsing database is completed. SQLite database is saved to: " + EXPORT_DB_LOCATION)
