'''
 :argument: DATA_FILE: from SLK 2.0 endocytosis data
'''

from ARNlib.SQLiteDBApi.sqlite_db_api import PsimiSQL

# Defining constants
SQL_SEED = '../../../../../ARNlib/SQLiteDBApi/network-db-seed.sql'
DB_TYPE = 'SignaLink3'
EXPORT_DB_LOCATION = '../../output/signalink3'
DATA_FILE = 'files/SLK3_human_core.csv'


def get_node_a(name, taxid, alt_acc, pathway, topology, psi_mi_to_sql_object):
    """
    This function sets up a node dict and returns it. If the node is already in the SQLite database it fetches that node from the db, so it won't be inserted multiple times
    """

    # Testing if the node is already in the database
    node_dict = psi_mi_to_sql_object.get_node(name, node_tax_id=taxid)

    if not node_dict:
        node_dict = {
            "name" : 'Uniprot:' + name,
            "tax_id": taxid,
            "alt_accession": alt_acc,
            'pathways': pathway,
            "aliases": None,
            "topology": topology
        }

    return node_dict


def get_node_b(name, taxid, alt_acc, pathway, topology, psi_mi_to_sql_object):
    """
    This function sets up a node dict and returns it. If the node is already in the SQLite database it fetches that node from the db, so it won't be inserted multiple times.

    """

    # Testing if the node is already in the database
    node_dict = psi_mi_to_sql_object.get_node(name, node_tax_id=taxid)

    if not node_dict:
        node_dict = {
            "name": 'Uniprot:' + name,
            "tax_id": taxid,
            "alt_accession": alt_acc,
            'pathways': pathway,
            "aliases": None,
            "topology": topology
        }

    return node_dict

def get_mitab_pathways_list(pathways):
    pathways_list = pathways.split('|')
    new_pathway_list = []
    for p in pathways_list:
        new_pathway = p.split("(")[0]
        if new_pathway == "Hedgehog":
            new_pathway_short = "HH"
            new_pathway_list.append(new_pathway_short)
        elif new_pathway == "NOTCH":
            new_pathway_short = "Notch"
            new_pathway_list.append(new_pathway_short)
        else:
            new_pathway_list.append(new_pathway)

    return new_pathway_list

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
    "WNT": 'WNT/Wingless',
    "HIPPO": 'Hippo',
    "HH": 'Hedgehog',
    "TNF/Apoptosis": 'TNF pathway',
    "TNF": 'TNF pathway',
}

def main(logger):
    # Initiating the parser
    db_api = PsimiSQL(SQL_SEED)

    with open(DATA_FILE) as data:

        # Skipping the header
        data.readline()

        for line in data:
            columns = line.strip().split(',')

            pathways_a = get_mitab_pathways_list(columns[5])
            new_pathways_a = []
            for p in pathways_a:
                pathway_a = p
                if " " in p:
                    pathway_a = p.replace(" ", "")
                elif '"' in p:
                    pathway_a = p.replace('"', "")
                new_p = accepted_pathways[pathway_a]
                new_pathways_a.append(new_p)
            pathways_b = get_mitab_pathways_list(columns[11])
            new_pathways_b = []
            for p in pathways_b:
                pathway_b = p
                if " " in p:
                    pathway_b = p.replace(" ", "")
                elif '"' in p:
                    pathway_b = p.replace('"', "")
                new_p_b = accepted_pathways[pathway_b]
                new_pathways_b.append(new_p_b)

            if columns[4] == 'None':
                topo_a = None
            else:
                topo_a = columns[4]

            if columns[10] == 'None':
                topo_b = None
            else:
                topo_b = columns[10]

            # Creating the node dicts, if the node is already in the db assigning that to the node dict
            source_dict = get_node_a(columns[1], 'taxid:9606', columns[0], "|".join(new_pathways_a), topo_a, db_api)
            target_dict = get_node_b(columns[7], 'taxid:9606', columns[6], "|".join(new_pathways_b), topo_b, db_api)

            # Nodes are inserted to the db if they are not in it yet
            if not 'id' in source_dict:
                db_api.insert_node(source_dict)

            if not 'id' in target_dict:
                db_api.insert_node(target_dict)

            # Pubmed references
            pub_id = '|pubmed:'.join(columns[16].split('|'))

            # Directedness
            effect = columns[13]

            interaction_types = "is_direct:true|is_directed:true|%s" \
                                % effect

            edge_dict = {
                'publication_ids': columns[16],
                'layer': '3',
                'source_db': 'SignaLink3',
                'interaction_identifiers': None,
                'confidence_scores': None,
                'interaction_detection_method': None,
                'interaction_types': interaction_types,
                'first_author': None
            }

            db_api.insert_edge(source_dict, target_dict, edge_dict)

    # Saving the to a DB_TYPE.db files
    db_api.save_db_to_file(EXPORT_DB_LOCATION)


if __name__ == '__main__':
    print("Parsing database...")
    main(logger = None)
    print("Parsing database is completed. SQLite database is saved to: " + EXPORT_DB_LOCATION)



