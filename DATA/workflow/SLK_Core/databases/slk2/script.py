"""
This script parses a .csv files with the curated data
    :argument: FILE_LOCATION: location of the .csv files
    :argument: DB_TYPE: name of the source database
    :argument: DB_DESTINATION: saving destination
    :argument: ORGANISM_NAME_TO_MITAB_ID_MAP: organism name to taxid
    :argument: IS_DIRECT_MAP: dictionary of directness to MI ids
    :argument: IS_DIRECTED_MAP: dictionary to directedness to MI ids
    :argument: EFFECT_MAP: dictionary of effect to MI ids
    :argument: MOLECULAR_MAP: dictionary of molecular background to MI ids

"""

import csv
from ARNlib.SQLiteDBApi.sqlite_db_api import PsimiSQL

#Defining constants
#SQL_SEED = '../../../../../ARNlib/SQLiteDBApi/network-db-seed.sql'
DB_TYPE = 'slk2'
FILE_LOCATION = 'files/04122017-signalink-1uBbAI.csv'
DB_DESTINATION = '../../../layer0/output/slk2'
ORGANISM_NAME_TO_MITAB_ID_MAP = {
    "H. sapiens": "taxid:9606",
    "D. melanogaster": "taxid:7227",
    "C. elegans": "taxid:6239"
}
IS_DIRECT_MAP = {
    "direct": "MI:0407(direct interaction)",
    "predicted to be direct" : "MI:0407(direct interaction)",
    "indirect": "indirect",
    "UNK": "unknown"
}

IS_DIRECTED_MAP = {
    "ppi directed": 'directed'
}

EFFECT_MAP = {
    'Unknown': 'MI:0190(interaction type)',
    'down-regulates': 'MI:2240(down-regulates)',
    "down-regulates activity": 'MI:2241(down-regulates activity)',
    "down-regulates quantity": 'MI:2242(down-regulates quantity)',
    "down-regulates quantity by destabilization": 'MI:2244(down-regulates quantity by destabilization)',
    "down-regulates quantity by repression": 'MI:2243(down-regulates quantity by repression)',
    'unknown': 'MI:0190(interaction type)',
    'up-regulates': 'MI:2235(up-regulates)',
    "up-regulates activity": 'MI:2236(up-regulates activity)',
    "up-regulates quantity by expression": 'MI:2238(up-regulates quantity by expression)',
    "up-regulates quantity by stabilization": 'MI:2239(up-regulates quantity by stabilization)',
    "stimulation" : 'MI:0624(stimulation)',                 # MI id from stimulant tag
    "inhibition" : 'MI:0623(inhibition)'
}

accepted_pathways = {
    "TCR": 'T-cell receptor',
    "BCR": 'B-cell receptor',
    "TLR": 'Toll-like receptor',
    "IIP": 'Innate immune pathways',
    "JAK/STAT": 'JAK/STAT',
    "NHR": 'Nuclear hormone receptor',
    "RTK": 'Receptor tyrosine kinase',
    "Rho/Cytoskeleton": 'Rho pathway',
    "TGF": 'TGF',
    "Notch": 'Notch',
    "GPCR": 'G-protein coupled receptor',
    "WNT/Wingless": 'WNT/Wingless',
    "HIPPO": 'Hippo',
    "HH": 'Hedgehog',
    "TNF/Apoptosis": 'TNF pathway',
}


def get_mitab_pathways_list_string(pathways):
    pathways_list = pathways.split(',')
    new_pathway_list = []
    for p in pathways_list:
        new_pathway = p.split("(")[0]
        if new_pathway == "Hedgehog":
            new_pathway_short = "HH"
            new_pathway_list.append(new_pathway_short)
        else:
            new_pathway_list.append(new_pathway)

    return "|".join(new_pathway_list)


def get_mitab_publication_list_string(publications):
    publications_list = publications.split('|')
    publications_list.append("23331499")
    return "pubmed:"+'|pubmed:'.join(publications_list)


def insert_or_get_node_dict(node_name, alias, pathways, topology, taxname, node_names_to_id, db_api):

    mitab_name = "Uniprot:" + node_name

    topologies = set(map(lambda x: x.strip(), topology.split(",")))

    pathway = pathways.split("|")
    new_pathways = []
    for p in pathway:
        new_p = accepted_pathways[p]
        new_pathways.append(new_p)

    node_dict = {
        "name": mitab_name,
        "tax_id": ORGANISM_NAME_TO_MITAB_ID_MAP[taxname],
        "alt_accession": None,
        "aliases": alias,
        "pathways" : "|".join(new_pathways),
        "topology" : "|".join(topologies)
    }

    if node_dict['name'] in node_names_to_id:
        node_dict['id'] = node_names_to_id[node_dict['name']]
    else:
        db_api.insert_unique_node(node_dict)
        node_names_to_id[node_dict['name']] = node_dict['id']

    return node_dict


# Declaring variables and constants
inserted_nodes = {}


def main(logger):
    # Initiating the parser
    db_api = PsimiSQL(SQL_SEED)

    FILE = csv.reader(open(FILE_LOCATION), delimiter=';')

    # Skipping the header
    next(FILE)

    node_names_to_id = {}

    for row in FILE:
        mitab_source_pathways = get_mitab_pathways_list_string(row[5])
        mitab_target_pathways = get_mitab_pathways_list_string(row[11])

        # Creating the node dicts, if the node is already in the db assigning that to the node dict
        source_dict = insert_or_get_node_dict(row[1],row[0],mitab_source_pathways,row[4].strip(), row[3], node_names_to_id, db_api)
        target_dict = insert_or_get_node_dict(row[7],row[6],mitab_target_pathways,row[10].strip(), row[9], node_names_to_id, db_api)

        effect = EFFECT_MAP[row[15]]

        is_direct = IS_DIRECT_MAP[row[14].lower()]
        if "MI:0407(direct interaction)" in is_direct:
            is_direct = "true"
        else:
            is_direct = "false"

        is_directed = IS_DIRECTED_MAP[row[13].lower()]
        if is_directed == "directed":
            is_directed = "true"
        else:
            is_directed = "false"

        # Setting up the interaction type
        interaction_types = "%s|is_directed:%s|is_direct:%s" \
                            % (effect, is_directed, is_direct)

        new_scores = []

        if row[18] != '':
            scores = row[18].split(",")
            for s in scores:
                confidence_score_name = s.split(":")[0]
                if " " in confidence_score_name:
                    confidence_score_name = confidence_score_name.replace(" ", "")
                confidence_score_value = s.split(":")[1]
                if " " in confidence_score_value:
                    confidence_score_value = confidence_score_value.replace(" ", "")
                score = f'SLK2 {confidence_score_name}:{confidence_score_value}'
                if score not in new_scores:
                    new_scores.append(score)

        edge_dict = {
            'interaction_detection_method' : None,
            'first_author': None,
            'publication_ids': get_mitab_publication_list_string(row[16]),
            'interaction_types': interaction_types,
            'source_db' : "SLKv2.0",
            'interaction_identifiers': None,
            'confidence_scores': "|".join(new_scores),
            'layer' : "8"
        }

        db_api.insert_edge(source_dict, target_dict, edge_dict)

    # Saving the to a DB_TYPE.db files
    db_api.save_db_to_file(DB_DESTINATION)


if __name__ == '__main__':
    main(logger=None)
