"""
    This script parses the reactome human database Psi-Mi formatted .tsv files and copies the acquired data to a SQLite .db files
    :argument: DATA_FILE: human PPI pairs in psi-mitab format: http://www.reactome.org/download/current/homo_sapiens.mitab.interactions.txt.gz
    :argument: DB_TYPE: name of the source database
    :argument: PATHWAY_FILE_LOCATION: files with 2 columns, source pathway id and signalink id
    :argument: PUBMED_INTERACTION_FILE:human PPI pairs in tab deliminated format: http://www.reactome.org/download/current/homo_sapiens.interactions.txt.gz
    :argument: EXPORT_DB_LOCATION: saving location of the created database
    :argument: UNI_TO_PATHWAY: files with uniprot ids and their pathways from: https://reactome.org/download/current/UniProt2Reactome_All_Levels.txt
"""

from ARNlib.SQLiteDBApi.sqlite_db_api import PsimiSQL
import sys
import sqlite3
#Defining constants
SQL_SEED = '../../../../../ARNlib/SQLiteDBApi/network-db-seed.sql'
DATA_FILE = '../../../layer0/databases/reactome/files/reactome.homo_sapiens.interactions.psi-mitab.txt'
DB_TYPE = 'Reactome'
PATHWAY_FILE_LOCATION = '../../../layer0/databases/reactome/files/pathway_map.txt'
EXPORT_DB_LOCATION = '../../../layer0/output/reactome'
UNI_TO_PATHWAY = 'files/UniProt2Reactome_All_Levels.txt'


def insert_or_get_node_dict(id, idtype, alt_id, alias, taxid, pw, node_names_to_id, db_api):
    if idtype == "uniprotkb":
        idtype = "Uniprot"

    if pw == "TCR":
        pw = "T-cell receptor"
    elif pw == "BCR":
        pw = "B-cell receptor"
    elif pw == "TLR":
        pw = "Toll-like receptor"
    elif pw == "IIP":
        pw = "Innate immune pathways"
    elif pw == "NHR":
        pw = "Nuclear hormone receptor"
    elif pw == "RTK":
        pw = "Receptor tyrosine kinase"
    elif pw == "Rho/Cytoskeleton":
        pw = "Rho pathway"
    elif pw == "GPCR":
        pw = "G-protein coupled receptor"
    elif pw == "WNT":
        pw = "WNT/Wingless"
    elif pw == "HIPPO":
        pw = "Hippo"
    elif pw == "HH":
        pw = "Hedgehog"

    node_dict = {
        "name": idtype.strip() + ':' + id.strip(),
        "tax_id": taxid,
        "alt_accession": alt_id,
        'pathways': pw,
        "aliases": alias,
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
    path_dict = {}

    # Uncomment if using run_all_auto
    # data = DATA_FILE.split('\n')

    # Getting only human data
    # path_map = UNI_TO_PATHWAY.split('\n')

    with open(DATA_FILE) as data, open(UNI_TO_PATHWAY) as path_map:

        for line in path_map:

            line = line.strip().split('\t')

            if len(line) > 4:
                if line[5] == 'Homo sapiens':
                    path_dict[line[0]] = line[1]

        data.readline()

        reactome_to_signalink_pathway_map = {}

        pathway_file = open(PATHWAY_FILE_LOCATION)
        next(pathway_file)

        for line in pathway_file:
            reactome_pathway_id, signalink_pathway = line.strip().split('\t')
            reactome_to_signalink_pathway_map[reactome_pathway_id] = signalink_pathway

        node_names_to_id = {}
        for line in data:

            columns = line.strip().split('\t')

            if len(columns) > 1:

                id_a = columns[0].strip().split(":")[1]
                id_type_a = columns[0].strip().split(":")[0]
                id_b = columns[1].strip().split(":")[1]
                id_type_b = columns[1].strip().split(":")[0]
                # Building the pathway dict for SLK3 pathways

                if not id_a in path_dict.keys() or not id_b in path_dict.keys():
                    continue

                if not path_dict[id_a] in reactome_to_signalink_pathway_map \
                        or not path_dict[id_b] in reactome_to_signalink_pathway_map:
                    continue

                interactor_a_tax_id = columns[9].split("(")[0]
                interactor_b_tax_id = columns[10].split("(")[0]
                if (interactor_a_tax_id != "taxid:9606") or (interactor_b_tax_id != "taxid:9606"):
                    continue

                # Creating the node dicts, if the node is already in the db assigning that to the node dict
                source_dict = insert_or_get_node_dict(id_a, id_type_a, columns[2], columns[4], interactor_a_tax_id, reactome_to_signalink_pathway_map[path_dict[id_a]], node_names_to_id, db_api)
                target_dict = insert_or_get_node_dict(id_b, id_type_b, columns[3], columns[5], interactor_b_tax_id, reactome_to_signalink_pathway_map[path_dict[id_b]], node_names_to_id, db_api)

                # Setting up the interaction type
                effect = columns[11].replace('psi-mi:', '').replace('"', '')
                interaction_types = "%s|is_directed:%s|is_direct:%s" \
                                    % (effect, 'true', 'false')

                if columns[8] != '-':
                    pubmed = columns[8].split("|")
                    pubmed.append("pubmed:29145629")
                    pubmed_ids = "|".join(pubmed)
                else:
                    pubmed_ids = "pubmed:29145629"

                edge_dict = {
                    'publication_ids': pubmed_ids,
                    'layer': '0',
                    'source_db': 'Reactome',
                    'interaction_identifiers': None,
                    'confidence_scores': columns[14].split("(")[0],
                    'interaction_detection_method': columns[6].replace('psi-mi:', '').replace('"', ''),
                    'interaction_types': interaction_types,
                    'first_author': columns[7]
                }

                db_api.insert_edge(source_dict, target_dict, edge_dict)

    # Saving the to a DB_TYPE.db files
    db_api.save_db_to_file(EXPORT_DB_LOCATION)


if __name__ == '__main__':
    print("Parsing database...")
    main(logger=None)
    print("Parsing database is completed")
