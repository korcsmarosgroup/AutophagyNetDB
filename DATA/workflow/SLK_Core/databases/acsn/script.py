# -*- coding: utf-8 -*-

"""
    This script parses data from the ACSN database.
    :argument: SQLITE_DB_API: location of the sqlite api
    :argument: SQL_SEED: sql seed files location
    :argument: PATHWAY_FILE: PPI interactions in .gmt format https://acsn.curie.fr/files/acsn_master_curated.gmt
    :argument: ALL_EDGE_FILE_LOCATION: Correspondance between ACSN entities and HUGO names in .gmt format: https://acsn.curie.fr/files/acsn_names.gmt
    :argument: EXPORT_DB_LOCATION: Saving location of the created database
    :argument: CURATED_PROTEIN_LIST_FILE_LOCATION: PPI interaction in .sif format https://acsn.curie.fr/files/acsn_ppi.sif
    :argument: DB_TYPE: name of the source db
"""

from ARNlib.SQLiteDBApi.sqlite_db_api import PsimiSQL
import requests

#Defining constants
SQL_SEED = '../../../../../ARNlib/SQLiteDBApi/network-db-seed.sql'
PATHWAY_FILE_LOCATION = 'files/pathways.tsv'
EXPORT_DB_LOCATION = '../../../layer0/output/acsn'
CURATED_PROTEIN_LIST_FILE_LOCATION = 'files/acsn_ppi_ver2.txt'
ALL_EDGE_FILE_LOCATION = 'files/acsn_names.gmt'
DB_TYPE = 'ACSN'

PATHWAY_MAP = {
    "AKT_MTOR": "Receptor tyrosine kinase",
    "APOPTOSIS": "TNF pathway",
    "CASPASES": "TNF pathway",
    "TNF_RESPONSE": "TNF pathway",
    "HEDGEHOG": "Hedgehog",
    "MAPK": "Receptor tyrosine kinase",
    "PI3K_AKT_MTOR": "Receptor tyrosine kinase",
    "WNT_CANONICAL": "WNT/Wingless",
    "WNT_NON_CANONICAL": "WNT/Wingless",
    "CYTOSKELETON_POLARITY": "Rho pathway"
}

EFFECT_MAP = {
    "CATALYSIS": "stimulation",
    "INHIBITION": "inhibition",
    "MODULATION": "none",
    "PHYSICAL_STIMULATION": "stimulation",
    "TRIGGER": "none",
    "UNKNOWN_CATALYSIS": "stimulation",
    "UNKNOWN_INHIBITION": "inhibition",
    "UNKNOWN_POSITIVE_INFLUENCE": "stimulation",
    "activates": "stimulation",
    "inhibits": "inhibition",
    "null": "unknown",
    "catalysis-precedes": "none",
    "controls-phosphorylation-of|controls-state-change-of": "none",
    "controls-state-change-of": "none",
    "controls-expression-of": "none"
}


# Getting pathways
def get_pathways(pathway_file):
    """
    This function parses the pathway files and returns each pathway as a dictionary.
    :param pathway_file:
    :return: A list of dictionaries. Each dictionary has a pathway name (string) property and a genes (list of strings) property
    """
    pathways_list = []

    for line in pathway_file:
        cells = line.strip().split('\t')

        if cells[0] in PATHWAY_MAP.keys():
            pathway = PATHWAY_MAP[cells[0]]
        else:
            pathway = ""

        pathway_dict = {'pathway_name': pathway, 'genes': cells[2:len(cells)]}
        pathways_list.append(pathway_dict)

    return pathways_list


def get_pathway_list(gene, pathways_list):
    """
    This function searches the given pathways for the given gene and returns the name of the pathways where the given gene name is listed
    :param gene: A gene name in HUGO (HGNC) format
    :type gene: string
    :param pathways_list: List of pathways and genes that play role in these pathways
    :type pathways_list: List of dictionaries
    :return: A list of pathway names where the given gene plays a role
    """

    contributed_pathways = []

    for pathway in pathways_list:
        if gene in pathway['genes']:
            if pathway['pathway_name'] not in contributed_pathways:
                contributed_pathways.append(pathway['pathway_name'])

    return contributed_pathways


def main(logger):
    # Making a set from the curated files, the set will contain the proteins that does map to a unique value
    # (has 1 ->* mapping)
    not_valid_node_set = set()
    valid_node_set = set()

    # Getting all nodes that contain a * character, that means that the node represents more than one molecule
    # Making a list from these not unique proteins
    with open(CURATED_PROTEIN_LIST_FILE_LOCATION) as curated_protein_list_file:
        for line in curated_protein_list_file:
            line = line.strip()
            cells = line.split('\t')
            if len(cells) > 4:
                not_valid_node_set.add(cells[0])
            else:
                # Collecting protein nodes
                valid_node_set.add(cells[0])

    # Collecting pathways from the pathway files
    PATHWAY_FILE = PATHWAY_FILE_LOCATION
    pathways = get_pathways(open(PATHWAY_FILE))

    # Initialising a PsimiTOSQL object
    parser = PsimiSQL(SQL_SEED)

    # Generating a dictionary that holds unique node objects, in the same time the node sql table is filled up
    nodes = {}
    edges = {}

    with open(CURATED_PROTEIN_LIST_FILE_LOCATION) as PPI_FILE:

        PPI_FILE.readline()

        # Although this is a SIF formatted files, it only contains two interactors in a line
        # (a SIF files can contain more than 2 interactors in a line)
        for line in PPI_FILE:
            # getting the names of interacting genes in HUGO format
            cells = line.strip().split('\t')

            try:

                inttype = cells[1]

                gene_a = cells[0]
                gene_b = cells[2]

                pubmed_ids = cells[3]

            except IndexError:
                continue

            if (gene_a not in valid_node_set) or (gene_b not in valid_node_set):
                continue

            if pubmed_ids:
                pubmed_list = pubmed_ids.split(';')
                pubmed_list.append("26192618")
                if 'N/A' in pubmed_list:
                    pubmed_list.remove('N/A')
                pubmed_ids = 'pubmed:'+'|pubmed:'.join(pubmed_list)

                edge_id = gene_a + '@' + gene_b

                for type in inttype.lower().split(';'):
                    final_inttype = []
                    if 'association' in type:
                        selected_type = 'MI:0914(association)'
                        final_inttype.append(selected_type)
                    else:
                        selected_type = 'MI:0190(interaction type)'
                        final_inttype.append(selected_type)

                    if edge_id not in edges:
                        edges[edge_id] = {
                            'inserted': False,
                            'is_complex': None,
                            'pubmed': pubmed_ids,
                            'effect': '|'.join(final_inttype)
                        }
            else:
                continue

    with open(CURATED_PROTEIN_LIST_FILE_LOCATION) as PPI_FILE:

        PPI_FILE.readline()

        for line in PPI_FILE:
            # Resetting variables
            edge_id = None
            gene_a = None
            gene_b = None
            effect = None
            edge_id = None

            try:
                cells = line.split('\t')

                gene_a = cells[0]
                gene_b = cells[2]

            except IndexError:
                continue

            not_accepted_characters = [" ", "?", "~", ","]
            characters_in_gene_a = [e for e in not_accepted_characters if e in gene_a]
            if len(characters_in_gene_a) > 0:
                continue
            characters_in_gene_b = [e for e in not_accepted_characters if e in gene_b]
            if len(characters_in_gene_b) > 0:
                continue

            if (gene_a not in valid_node_set) or (gene_b not in valid_node_set):
                continue

            edge_id = gene_a + '@' + gene_b

            if edge_id in edges:
                if edges[edge_id]['is_complex'] is True or edges[edge_id]['inserted'] is True or "Reference" in edges[edge_id]['effect'] or "neighbor-of" in edges[edge_id]['effect']:
                    continue
                else:
                    pubmed_ids = edges[edge_id]['pubmed']
                    effect = edges[edge_id]['effect']
            else:
                continue

            """ creating and inserting edges to the db """

            gene_a_pathway_list = get_pathway_list(gene_a.replace('*', ''), pathways)

            gene_b_pathway_list = get_pathway_list(gene_b.replace('*', ''), pathways)

            # If the node is in the not_valid_node set, it is not inserted
            if gene_a not in not_valid_node_set:
                gene_a = gene_a.replace('*', '')
                if gene_a in nodes:
                    interactor_a = nodes[gene_a]
                else:
                    interactor_a = nodes[gene_a] = {
                        'name': 'HGNC:' + gene_a,
                        'alt_accession' : 'HGNC:' + gene_a,
                        'tax_id' : 'taxid:9606',
                        'pathways': '|'.join(gene_a_pathway_list),
                        'aliases': None
                    }
                    parser.insert_node(interactor_a)
            else:
                continue

            if gene_b not in not_valid_node_set:
                gene_b = gene_b.replace('*','')
                if gene_b in nodes:
                    interactor_b = nodes[gene_b]
                else:
                    interactor_b = nodes[gene_b] = {
                        'name': 'HGNC:' + gene_b,
                        'alt_accession' : 'HGNC:' + gene_b,
                        'tax_id': 'taxid:9606',
                        'pathways': '|'.join(gene_b_pathway_list),
                        'aliases': None
                    }
                    parser.insert_node(interactor_b)
            else:
                continue

            interaction_types = "%s|is_directed:%s|is_direct:%s" \
                                % (effect, "true", "false")

            edge_dict = {
                'interaction_detection_method': None,
                'first_author': None,
                'publication_ids': pubmed_ids,
                'interaction_types': interaction_types,
                'source_db': DB_TYPE,
                'interaction_identifiers': None,
                'confidence_scores': None,
                'layer': "0"
            }

            parser.insert_edge(interactor_a, interactor_b, edge_dict)

            edges[edge_id]['inserted'] = True

    parser.save_db_to_file(EXPORT_DB_LOCATION)


if __name__ == '__main__':
    print("Parsing database...")
    main()
    print("Parsing database is completed. SQLite database is saved to: " + EXPORT_DB_LOCATION)
