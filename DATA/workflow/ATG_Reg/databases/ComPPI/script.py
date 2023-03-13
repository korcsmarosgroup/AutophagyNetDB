"""
 Parse comPPI database, this will be additional edge information that will be added later on to the layers
    :argument: DATA_FILE_LIST: list of the data files of all species
    :argument: EXPORT_DB_LOCATION: saving location
    :argument: DB_TYPE: name of the database
"""

# Imports
from ARNlib.SQLiteDBApi.sqlite_db_api import PsimiSQL
import re

# Defining constants
SQL_SEED = '../../../../../ARNlib/SQLiteDBApi/network-db-seed.sql'
DATA_FILE_LIST = ['files/comppi--interactions--tax_hsapiens_loc_all.txt']
EXPORT_DB_LOCATION = '../../output/ComPPI'
DB_TYPE = 'comppi'
MI_TREM_ONTOLOGY_FILE = '../mi.owl'


def insert_or_get_node_dict(id, idtype, alt_id, taxid, node_names_to_id, db_api):
    idtype = idtype.split("/")[0]
    if idtype == "UniProtKB":
        idtype = "Uniprot"
    else:
        print("warning: id type=" + idtype)

    node_dict = {
        "name": idtype.strip() + ':' + id.strip(),
        "tax_id": 'taxid:' + taxid,
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


def build_mi_term_dict():
    mi_term_dict = {}
    with open(MI_TREM_ONTOLOGY_FILE) as f:
        mi_id = ""
        for line in f:
            if line.startswith("id: "):
                mi_id = line[4:].strip()
            if line.startswith("name: "):
                name = line[6:].strip()
                mi_term_dict[name.lower()] = "%s(%s)" % (mi_id, name)

    # a few extra ComPII specific items:
    mi_term_dict['surface plasmon resonance array (spr)'] = "MI:0921(surface plasmon resonance array)"
    mi_term_dict['proximity enzyme-linked immunosorbent assay'] = "MI:0411(enzyme linked immunosorbent assay)"
    mi_term_dict['co-immunoprecipitation'] = "MI:0019(coimmunoprecipitation)"
    mi_term_dict['lambda repressor two-hybrid'] = "MI:0655(lambda repressor two hybrid)"
    mi_term_dict['reverse two-hybrid'] = "MI:0726(reverse two hybrid)"
    mi_term_dict['mammalian protein-protein interaction trap'] = "MI:0231(mammalian protein protein interaction trap)"
    mi_term_dict['atomic force microscopy (afm)'] = "MI:0872(atomic force microscopy)"
    mi_term_dict['beta-lactamase complementation'] = "MI:0011(beta lactamase complementation)"
    mi_term_dict['far-western blotting'] = "MI:0047(far western blotting)"
    mi_term_dict['pull-down'] = "MI:0096(pull down)"
    mi_term_dict['anti-bait co-immunoprecipitation'] = "MI:0006(anti bait coimmunoprecipitation)"
    mi_term_dict['acetylation assay'] = "MI:0889(acetylase assay)"
    mi_term_dict['fluorescence resonance energy transfer (fret)'] = "MI:0055(fluorescent resonance energy transfer)"
    mi_term_dict['two-hybrid fragment pooling approach'] = "MI:0399(two hybrid fragment pooling approach)"
    mi_term_dict['gal4-vp16 complementation'] = "MI:0728(gal4 vp16 complementation)"
    mi_term_dict['experimental'] = "MI:0045(experimental interaction detection)"
    mi_term_dict['three-hybrid screening'] = "MI:0588(three hybrid)"
    mi_term_dict['x-ray diffraction'] = "MI:0114(x-ray crystallography)"
    mi_term_dict['bioluminescence resonance energy transfer (bret)'] = "MI:0012(bioluminescence resonance energy transfer)"
    mi_term_dict['anti-tag co-immunoprecipitation'] = "MI:0007(anti tag coimmunoprecipitation)"
    mi_term_dict['lexa dimerization assay'] = "MI:0369(lex-a dimerization assay)"
    mi_term_dict['two-hybrid screening'] = "MI:0018(two hybrid)"
    mi_term_dict['immunodepleted co-immunoprecipitation'] = "MI:0858(immunodepleted coimmunoprecipitation)"
    mi_term_dict['nuclear magnetic resonance (nmr)'] = "MI:0077(nuclear magnetic resonance)"
    mi_term_dict['toxr dimerization assay'] = "MI:0370(tox-r dimerization assay)"
    mi_term_dict['enzyme-linked immunosorbent assay (elisa)'] = "MI:0411(enzyme linked immunosorbent assay)"
    mi_term_dict['oxidoreductase activity electron transfer assay'] = "MI:0979(oxidoreductase assay)"
    mi_term_dict['kinase homogeneous time-resolved fluorescence'] = "MI:0420(kinase homogeneous time resolved fluorescence)"
    mi_term_dict['mst'] = "MI:1247(microscale thermophoresis)"
    mi_term_dict['homogeneous time-resolved fluorescence'] = "MI:0510(homogeneous time resolved fluorescence)"
    mi_term_dict['beta-galactosidase complementation'] = "MI:0010(beta galactosidase complementation)"
    mi_term_dict['comigration in sds-page'] = "MI:0808(comigration in sds page)"
    mi_term_dict['solid-phase assay'] = "MI:0892(solid phase assay)"

    # mi_term_dict['vt'] = ???
    # mi_term_dict['yeast t'] = ??? "MI:0018(two hybrid)"

    return mi_term_dict


def main(logger):
    # Initiating the parser
    db_api = PsimiSQL(SQL_SEED)
    unknown_det_methods = set()
    unknown_det_method_count = 0
    mi_term_dict = build_mi_term_dict()

    node_names_to_id = {}
    for file_loc in DATA_FILE_LIST:
        print("parsing ComPPI file: " + file_loc)
        with open(file_loc) as data:
            # Skipping the header
            data.readline()
            lines = 0

            for line in data:
                lines += 1
                if lines % 50000 == 0:
                    print("processed lines: %d" % lines)
                columns = line.split('\t')

                # Creating the node dicts, inserted to the db if they are not in it yet
                source_dict = insert_or_get_node_dict(columns[0], columns[1], columns[2], columns[3], node_names_to_id,
                                                      db_api)
                target_dict = insert_or_get_node_dict(columns[4], columns[5], columns[6], columns[7], node_names_to_id,
                                                      db_api)

                # Interaction detection method
                detmethods = []
                for detection_method in columns[9].split('|'):
                    detection_method = detection_method.replace("(Unknown)", "")
                    detection_method = detection_method.replace("(Experimental)", "")
                    detection_method = detection_method.replace("\"", "")
                    detection_method = detection_method.strip().lower()
                    if detection_method in mi_term_dict:
                        detmethods.append(mi_term_dict[detection_method])
                    else:
                        unknown_det_methods.add(detection_method)
                        unknown_det_method_count += 1

                interaction_types = "is_directed:false|is_direct:true"

                pubmed_ids = map(lambda x: "pubmed:" + x.strip(), columns[11].split('|'))
                pubmed_ids = filter(lambda x: re.search("^\\d+$", x), pubmed_ids)
                pubmed_ids = set(pubmed_ids)
                pubmed_ids.add("pubmed:25348397") # ComPPI publication

                edge_dict = {
                    'publication_ids': '|'.join(pubmed_ids),
                    'layer': '2',
                    'source_db': 'ComPPI',
                    'interaction_identifiers': None,
                    'confidence_scores': None,
                    'interaction_detection_method': "|".join(detmethods),
                    'interaction_types': interaction_types,
                    'first_author': None
                }

                db_api.insert_edge(source_dict, target_dict, edge_dict)
            print("processed lines: %d" % lines)

    print("not identified detection methods: \n%s\n\n" % "\n".join(unknown_det_methods))
    print("number of not identified detection methods: %d" % unknown_det_method_count)

    # Saving the to a DB_TYPE.db files
    print("exporting SQLite database to: " + EXPORT_DB_LOCATION)
    db_api.save_db_to_file(EXPORT_DB_LOCATION)


if __name__ == '__main__':
    print("Parsing database...")
    main(logger=None)
    print("Parsing database is completed. SQLite database is saved to: " + EXPORT_DB_LOCATION)
