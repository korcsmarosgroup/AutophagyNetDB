# setting up imports
import sqlite3
import sys, logging

# Adding the psimi_to_sql module to sys
from ARNlib.SQLiteDBApi.sqlite_db_api import PsimiSQL

# Declaring constants
SQL_SEED_LOCATION = '../SQLiteDBApi/network-db-seed.sql'
DESTINATION = '../merger'


def merge_strings(string_1, string_2, separator="|"):
    """
    This function merges two stings, that elements separated by a string or character.
    :param string_1:
    :type string_1: str
    :param string_2:
    :param separator:
    :type string_2: str
    :return:
    """

    if string_1 is None:
        string_1 = ""

    if string_2 is None:
        string_2 = ""

    # converting the strings to lists
    list_1 = string_1.split(separator)
    list_2 = string_2.split(separator)

    list_1 = list(filter(lambda item: item != '-', list_1))
    list_2 = list(filter(lambda item: item != '-', list_2))

    # adding the lists
    sum_lists = list_1 + list_2

    # 1) converting the lists to sets 2) removing empty elements
    # 3) getting the union of these sets 4) converting back to list
    merged_lists = list(set(filter(None, sum_lists)))

    # joining the list to a string
    result_string = separator.join(merged_lists)

    return result_string


def get_union_of_nodes(node_dict_1, node_dict_2):
    """
    This function returns the union of two node dicts properties, except the alt_accession and the tax_id. If the alt accessions are not the same the first node's alt accession will be the new alt accession, and the second node's alt accession becomes an alias.
    :param node_dict_1:
    :param node_dict_2:
    :return:
    """

    # merging the two nodes properties that needs to be merged

    merged_pathways = merge_strings(node_dict_1['pathways'], node_dict_2['pathways'])

    merged_aliases = merge_strings(node_dict_1["aliases"], node_dict_2["aliases"])
    merged_topology = merge_strings(node_dict_1["topology"], node_dict_2["topology"])

    alt_acccession = node_dict_1['alt_accession']

    # dealing with the alt accession (gene symbol)

    if node_dict_1['alt_accession'] != node_dict_2['alt_accession']:
        if merged_aliases:
            merged_aliases += "|" + node_dict_2['alt_accession']
        else:
            merged_aliases = node_dict_2['alt_accession']

    # initiating the new node
    new_node = {
        "name": node_dict_1["name"],
        "alt_accession": alt_acccession,
        "tax_id": node_dict_1["tax_id"],
        "pathways": merged_pathways,
        "aliases": merged_aliases,
        "topology": merged_topology
    }

    return new_node


def main(log):
    print("The following dbs will be parsed: " + ", ".join(map(lambda line: line.split("/")[-1], SOURCE_DB_FILE_LIST)))

    # Declaring the dicts that will hold the data
    nodes = {}
    collected_edges_undirected = {}
    collected_edges_directed = {}

    merged_edge_counter = 0
    all_edge_counter = 0
    directed_counter = 0
    undirected_counter = 0

    # the number of pieces (.db files)
    sum_files = len(SOURCE_DB_FILE_LIST)

    # filling up the nodes dictionary with the data contained in db piece files
    for db_file in SOURCE_DB_FILE_LIST:
        print("parsing DB file: " + db_file)
        # executing a query that selects everything (but the node id) from the current SQLite .db files
        db = sqlite3.connect(db_file)
        cursor = db.cursor()
        cursor.execute("SELECT * FROM node WHERE tax_id = 'taxid:9606' OR "
                       "tax_id = 'taxid:6239' OR tax_id = 'taxid:7227' OR tax_id = 'taxid:7955'")

        # iterating trough the db row by row
        while True:
            row = cursor.fetchone()
            # until the last row
            if row is None:
                break
            # if unique, inserting the node (row) to the nodes dictionary
            id, name, alt_accession, tax_id, pathways, aliases, topology = row
            node = {
                "name" : name,
                'alt_accession' : alt_accession,
                'tax_id' : tax_id,
                'pathways' : pathways,
                'aliases' : aliases,
                'topology' : topology
            }

            if name not in nodes:
                nodes[name] = node
            else:
                nodes[name] = get_union_of_nodes(nodes[name], node)

        # closing the current db
        db.close()
        # logging out some info
        current_file = SOURCE_DB_FILE_LIST.index(db_file)
        print("Building the node dictionary: Processing %d files out of %d\r" % (current_file, sum_files))
        print("total number of nodes: %d" % len(nodes))

    print("all input db processed successfully\n\n")
    # making a memory database and inserting the unique nodes from the nodes dictionary
    #log.debug('Inserting nodes to database')
    parser = PsimiSQL(SQL_SEED_LOCATION)
    for node in nodes:
        parser.insert_unique_node(nodes[node])
        nodes[node]['id'] = parser.cursor.lastrowid

    # looping through the files again to make an edge list
    #log.debug("Started building edge dict")
    file_counter = 1
    for db_file in SOURCE_DB_FILE_LIST:
        print("\nprocessing edges from file: " + db_file)
        # executing a query that selects everything (but the node id) from the current SQLite .db files
        db = sqlite3.connect(db_file)
        cursor = db.cursor()
        cursor.execute("SELECT * FROM edge")

        while True:
            row = cursor.fetchone()
            # if there aren't any more nodes break out of the loop
            if not row:
                print("finished... total number of edges so far: %d" % all_edge_counter)
                break
            else:
                all_edge_counter += 1
                # deconstructing the row (list)
                edge_row_id, old_interactor_a_node_id, old_interactor_b_node_id, interactor_a_node_name, interactor_b_node_name, interaction_detection_method , first_author, publication_ids, interaction_types, source_db, interaction_identifiers, confidence_scores, layer = row

                # because in the nodes dict building process the query only asks for human||drosophila||celegans nodes
                # we have to make sure that we don't try to insert edges whose nodes are in the nodes dict
                # (=does not contain any other organisms node)
                if interactor_a_node_name in nodes and interactor_b_node_name in nodes:
                    # generating an edge id that will be the key in the edge dict
                    edge_id = interactor_a_node_name + "@" + interactor_b_node_name + "@" + str(layer)
                    edge_id_reverse = interactor_b_node_name + "@" + interactor_a_node_name + "@" + str(layer)

                    # generating an edge dict, that will be a value for the key
                    current_edge = {
                        'interaction_detection_method' : interaction_detection_method,
                        'first_author': first_author,
                        'publication_ids': publication_ids,
                        'interaction_types': interaction_types,
                        'source_db': source_db,
                        'interaction_identifiers': interaction_identifiers,
                        'confidence_scores': confidence_scores,
                        'layer': layer
                    }

                    # Checking if the edge is directed
                    if current_edge['interaction_types'] is not None:
                        directed = 'is_directed:true' in current_edge['interaction_types'].split('|')
                    else:
                        directed = False

                    if directed:
                        directed_counter += 1
                        # if the edge dict does not have this edge_id. the edge is stored in the edge dict with it's id
                        if not edge_id in collected_edges_directed:
                            collected_edges_directed[edge_id] = current_edge
                        else:
                            # if an edge is already in the dict it will be merged with the current edge
                            collected_edge = collected_edges_directed[edge_id]
                            collected_edge['interaction_types'] = merge_strings(collected_edge['interaction_types'],current_edge['interaction_types'])
                            collected_edge['first_author'] = merge_strings(collected_edge['first_author'],current_edge['first_author'])
                            collected_edge['publication_ids'] = merge_strings(collected_edge['publication_ids'],current_edge['publication_ids'])
                            collected_edge['source_db'] = merge_strings(collected_edge['source_db'], current_edge['source_db'])
                            collected_edge['interaction_identifiers'] = merge_strings(collected_edge['interaction_identifiers'], current_edge['interaction_identifiers'])
                            collected_edge['interaction_detection_method'] = merge_strings(collected_edge['interaction_detection_method'], current_edge['interaction_detection_method'])
                            collected_edge['confidence_scores'] = merge_strings(collected_edge['confidence_scores'], current_edge['confidence_scores'])
                            merged_edge_counter += 1

                    else:
                        # not directed:
                        # if the edge dict does not have this edge_id. the edge is stored in the edge dict with it's id
                        if not edge_id in collected_edges_undirected and not edge_id_reverse in collected_edges_undirected:
                            undirected_counter += 1
                            collected_edges_undirected[edge_id] = current_edge
                        else:
                            # if we are here, then we know that the undirected edge is present already
                            # first we find if it is present normally or reverse
                            if edge_id in collected_edges_undirected:
                                collected_edge = collected_edges_undirected[edge_id]
                            else:
                                collected_edge = collected_edges_undirected[edge_id_reverse]

                            collected_edge['interaction_types'] = merge_strings(collected_edge['interaction_types'],current_edge['interaction_types'])
                            collected_edge['first_author'] = merge_strings(collected_edge['first_author'],current_edge['first_author'])
                            collected_edge['publication_ids'] = merge_strings(collected_edge['publication_ids'], current_edge['publication_ids'])
                            collected_edge['source_db'] = merge_strings(collected_edge['source_db'], current_edge['source_db'])
                            collected_edge['interaction_identifiers'] = merge_strings(collected_edge['interaction_identifiers'], current_edge['interaction_identifiers'])
                            collected_edge['interaction_detection_method'] = merge_strings(collected_edge['interaction_detection_method'], current_edge['interaction_detection_method'])
                            collected_edge['confidence_scores'] = merge_strings(collected_edge['confidence_scores'], current_edge['confidence_scores'])
                            merged_edge_counter += 1

    print("\n\nTotal number of edges read from the input DBs: %d" % all_edge_counter)
    print("Total number of edges to be written to the output DB: %d" % (len(collected_edges_directed) + len(collected_edges_undirected)))
    print("Number of directed edges: %d" % directed_counter)
    print("Number of undirected edges: %d" % undirected_counter)
    print("Number of merged edges: %d" % merged_edge_counter)
    print("\n\nBuilding edge dict done!")
    print("Started inserting edges to the db")

    # iterating through edges dictionary and inserting nodes to the SQLite db
    for edge_id, edge_to_insert in collected_edges_undirected.items():
        node_a, node_b, layer = edge_id.split('@')
        parser.insert_edge(nodes[node_a], nodes[node_b], edge_to_insert)

    for edge_id, edge_to_insert in collected_edges_directed.items():
        node_a, node_b, layer = edge_id.split('@')
        parser.insert_edge(nodes[node_a], nodes[node_b], edge_to_insert)

    parser.save_db_to_file(DESTINATION)

    print("merging of the given layer is finished successfully :)")


if __name__ == '__main__':
    main(log=None)
