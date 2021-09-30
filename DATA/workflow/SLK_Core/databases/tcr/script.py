"""
This script parses the TCR signalling database .tsv files and inserts the extracted data to a PsiMI SQLite .db files
    :argument: TCR_DATA_LOC: data files
    :argument: DESTINATION: saving location
    :argument: DB_TYPE: name of the database
    :argument: IS_DIRECT_MAP: dictionary of directness to MI ids
    :argument: IS_DIRECTED_MAP: dictionary of directedness to MI ids
    :argument: EFFECT_MAP: dictionary of effect to MI ids
"""

from ARNlib.SQLiteDBApi.sqlite_db_api import PsimiSQL

#Defining constants
SQL_SEED = '../../../../../ARNlib/SQLiteDBApi/network-db-seed.sql'
TCR_DATA_LOC = 'files/T_cell_activation_signaling_pathway.txt'
DESTINATION = '../../../layer0/output/tcr'
DB_TYPE = 'tcr'

IS_DIRECT_MAP = {
    "direct": "true",
    "direct " : "true",
    "indirect": "false",
    "direct/indirect": "false",
    '2' : 'true'
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
    "inhibition": 'MI:0623(inhibition)',
    "inhibiton": 'MI:0623(inhibition)',
    "stimulation": 'MI:0624(stimulation)',
    "activation": 'MI:2254(chemical activation reaction)',

}


def get_mitab_publication_list_string(publications):
    publications_list = publications.split(';')
    if "23331499" not in publications_list:
        publications_list.append("23109003")
    return "pubmed:"+'|pubmed:'.join(publications_list)


def main(logger):
    TCR_DATA_FILE = open(TCR_DATA_LOC, encoding="ISO-8859-1")
    # Skipping the header line, and assigning the files's content to a list
    lines = TCR_DATA_FILE.readline()

    # Initiating a PsimiSQL object
    parser = PsimiSQL(SQL_SEED)

    for line in TCR_DATA_FILE:
        cells = line.split('\t')

        # Storing the needed properties in variables
        name_a = cells[1].strip()
        name_b = cells[3].strip()

        alt_accession_a = cells[0]
        alt_accession_b = cells[2]

        if name_a == '':
            continue

        # Building the node dictionaries, and inserting them to the db with the parser
        node_a_dict = {
            'name' : "Uniprot:"+name_a,
            'alt_accession' : "entrez gene/locuslink:"+alt_accession_a,
            'tax_id' : "taxid:9606",
            'pathways' : "T-cell receptor",
            'aliases' : None
        }

        parser.insert_node(node_a_dict)

        if name_b == '':
            continue

        node_b_dict = {
            'name' : "Uniprot:"+name_b,
            'alt_accession' : "entrez gene/locuslink:"+alt_accession_b,
            'tax_id' : "taxid:9606",
            'pathways' : "T-cell receptor",
            'aliases' : None
        }

        parser.insert_node(node_b_dict)

        # Gathering the edge's properies, and inserting the edge to the db

        interaction_direction = IS_DIRECT_MAP[cells[5].lower()]
        interaction_effect = EFFECT_MAP[cells[6].lower().strip()]

        pubmed_ids = cells[8]

        interaction_types = "%s|is_directed:%s|is_direct:%s" % (interaction_effect, "true", interaction_direction)

        edge_dict = {
            'interaction_detection_method' : None,
            'first_author' : None,
            'publication_ids' : get_mitab_publication_list_string(pubmed_ids),
            'interaction_types' : interaction_types,
            'source_db' : "TCRcuration",
            'interaction_identifiers' : None,
            'confidence_scores' : None,
            'layer' : '0'
        }

        parser.insert_edge(node_a_dict,node_b_dict,edge_dict)

    # Saving the db to a files
    parser.save_db_to_file(DESTINATION)


if __name__ == '__main__':
    main()
