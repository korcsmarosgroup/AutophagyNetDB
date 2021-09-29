"""
This script parses the signalink 2.1 csv files to a Psi-Mi formatted SQLite database files
    :argument: DB_TYPE: name of the database
    :argument: SLK_21_FILE_LOCATION: location of thedata files
    :argument: TAX_ID_MAP_FILE_LOCATION: files with two columns, uniprot id and tax id
    :argument: DB_DESTINATION: saving destination
    :argument: IS_DIRECT_MAP: dictionary of directness to MI ids
    :argument: IS_DIRECTED_MAP : dictionary of directedness to MI ids
    :argument: EFFECT_MAP: dictionary of effect to MI ids
    :argument: MOLECULAR_MAP: dictionary of molecular background to MI ids 
"""

# Imports
import csv
from SLKlib.SQLiteDBApi.sqlite_db_api import PsimiSQL

# Defining constants
SQL_SEED = '../../../../../SLKlib/SQLiteDBApi/network-db-seed.sql'
DB_TYPE = 'slk21'
SLK_21_FILE_LOCATION = 'files/SLK21.txt'
TAX_ID_MAP_FILE_LOCATION = 'files/entrytotaxid.txt'
DB_DESTINATION = '../../../layer0/output/slk21'

IS_DIRECT_MAP = {
    "direct": "MI:0407(direct interaction)",
    "dierct": "MI:0407(direct interaction)",
    "indirect": "indirect",
    "direct/indirect": "unknown",
    '2': 'MI:0407(direct interaction)'
}

IS_DIRECTED_MAP = {
    'directed': 'directed',
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
    "inhibition" : 'MI:0623(inhibition)',
    "inhibiton" : 'MI:0623(inhibition)',
    "stimulation" : 'MI:0624(stimulation)'
}

MOLECULAR_MAP = {
    'binding' : 'MI:0462(bind)',
    'transcriptional regulation' : 'MI:2247(transcriptional regulation)',
    'phosphorylation' : 'MI:0217(phosphorylation reaction)',
    'phosphorilation' : 'MI:0217(phosphorylation reaction)',
    'phosphoryltion' : 'MI:0217(phosphorylation reaction)',
    'phosphotilation' : 'MI:0217(phosphorylation reaction)',
    '' : 'MI:0190(interaction type)',
    'ubiquitination' : 'MI:0220(ubiquitination reaction)',
    'ubiquitinisation' : 'MI:0220(ubiquitination reaction)',
    'relocalization' : 'MI:2256(relocalization)',
    'dephosphorylation' : 'MI:0203(dephosphorylation  reaction)',
    'dephosphotilation' : 'MI:0203(dephosphorylation reaction)',
    'dephosphorilation' : 'MI:0203(dephosphorylation reaction)',
    'cleavage' : 'MI:0194(cleavage reaction)',
    'deubiquitination' : 'MI:0204(deubiquitination reaction)',
    'deubiquitinisation' : 'MI:0204(deubiquitination reaction)',
    'guanine nucleotide exchange factor' : 'MI:2252(guanine nucleotide exchange factor)',
    'sumoylation' : 'MI:0566(sumoylation reaction)',
    'dimethylation' : 'MI:0871(dimethylation reaction)',
    'acetylation' : 'MI:0192(acetylation reaction)',
    'deacetylation' : 'MI:0197(deacetylation reaction)',
    'deacetyltion' : 'MI:0197(deacetylation reaction)',
    'autophosphorylation' : 'MI:0217(phosphorylation reaction)',
    'autophosphorilation': 'MI:0217(phosphorylation reaction)',
    'methylation' : 'MI:0213(methylation reaction)',
    'monoubuquitinisation' : 'MI:0220(ubiquitination reaction)',
    'polyubiquitinisation' : 'MI:0220(ubiquitination reaction)',
    'autophosphorilation of mapk14' : 'MI:0217(phosphorylation reaction)'

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
        new_pathway = p.split("-")[0]
        if new_pathway == "Hedgehog":
            new_pathway_short = "HH"
            new_pathway_list.append(new_pathway_short)
        else:
            new_pathway_list.append(new_pathway)

    return "|".join(new_pathway_list)


def get_mitab_publication_list_string(publications):
    publications_list = publications.split(',')
    if "23331499" not in publications_list:
        publications_list.append("23331499")
    return "pubmed:" + '|pubmed:'.join(publications_list)


def get_taxid_map_dict(tax_id_map_file_location):
    """
    This function returns a dict where the keys are uniprot ids and the values are organism ids
    :param tax_id_map_file_location: The location of the files that will be the base of the uniprot to organism id dict
    :type tax_id_map_file_location: str
    :return: A dict k: uniprot id v: organism id
    """

    # Declaring the dict that will be returned
    tax_id_dict = {}

    # Opening the files and looping through it
    with open(tax_id_map_file_location) as tax_id_map_file:

        # Skipping header
        tax_id_map_file.readline()

        # Looping throgh files and adding entries to the tax_id_dict
        for line in tax_id_map_file:

            # Assigning the line to variables
            uniprot_id, tax_id = line.strip().split('\t')

            # Adding the entries to the dict
            tax_id_dict[uniprot_id] = tax_id

    return tax_id_dict


def insert_or_get_node_dict(node_map, node_name, pathways, topology, node_names_to_id, db_api):
    try:
        tax_id = node_map[node_name]
    except KeyError:
        tax_id = 'None'
        pass
    mitab_tax_id = "taxid:" + tax_id
    mitab_name = "Uniprot:" + node_name

    topologies = set(map(lambda x: x.strip(), topology.split(",")))
    if "#HIANYZIK" in topologies:
        topologies.remove("#HIANYZIK")

    pathway = pathways.split("|")
    new_pathways = []
    for p in pathway:
        new_p = accepted_pathways[p]
        new_pathways.append(new_p)

    node_dict = {
        "name": mitab_name,
        "tax_id": mitab_tax_id,
        "alt_accession": None,
        "aliases": None,
        "pathways" : "|".join(new_pathways),
        "topology" : "|".join(topologies)
    }

    if node_dict['name'] in node_names_to_id:
        node_dict['id'] = node_names_to_id[node_dict['name']]
    else:
        db_api.insert_unique_node(node_dict)
        node_names_to_id[node_dict['name']] = node_dict['id']

    return node_dict


def main(logger):
    # Declaring variables and constants
    inserted_nodes = {}
    UNIPROT_TO_TAX_ID_MAP = get_taxid_map_dict(TAX_ID_MAP_FILE_LOCATION)

    # Initiating the parser
    db_api = PsimiSQL(SQL_SEED)

    SLK_21_FILE = csv.reader(open(SLK_21_FILE_LOCATION, encoding="ISO-8859-1"), delimiter='\t', quotechar='"')
    # Skipping the header
    next(SLK_21_FILE)

    node_names_to_id = {}

    for row in SLK_21_FILE:

        mitab_source_pathways = get_mitab_pathways_list_string(row[1])
        mitab_target_pathways = get_mitab_pathways_list_string(row[4])

        if (row[0] not in UNIPROT_TO_TAX_ID_MAP) or (row[3] not in UNIPROT_TO_TAX_ID_MAP):
            continue

        # Creating the node dicts, if the node is already in the db assigning that to the node dict
        source_dict = insert_or_get_node_dict(UNIPROT_TO_TAX_ID_MAP, row[0], mitab_source_pathways, row[2], node_names_to_id, db_api)
        target_dict = insert_or_get_node_dict(UNIPROT_TO_TAX_ID_MAP, row[3], mitab_target_pathways, row[5], node_names_to_id, db_api)


        effect = EFFECT_MAP[row[8]]

        is_direct = IS_DIRECT_MAP[row[6].lower()]
        if "MI:0407(direct interaction)" in is_direct:
            is_direct = "true"
        else:
            is_direct = "false"

        is_directed = IS_DIRECTED_MAP[row[7].lower()]
        if is_directed == "directed":
            is_directed = "true"
        else:
            is_directed = "false"

        edge_dict = {
            'interaction_detection_method' : None,
            'first_author' : None,
            'publication_ids' : get_mitab_publication_list_string(row[9]),
            'interaction_types' : "%s|is_directed:%s|is_direct:%s" % (effect, is_directed, is_direct),
            'source_db' : 'SLKv2.1',
            'interaction_identifiers' : None,
            'confidence_scores' : None,
            'layer' : "0"
        }

        db_api.insert_edge(source_dict, target_dict, edge_dict)

    # Saving the to a DB_TYPE.db files
    db_api.save_db_to_file(DB_DESTINATION)


if __name__ == '__main__':
    main(logger=None)
