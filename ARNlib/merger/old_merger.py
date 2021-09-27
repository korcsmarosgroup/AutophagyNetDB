#setting up imports
import sqlite3
import sys
import os
import re

#adding the psimi_to_sql module to sys
from SLKlib.SQLiteDBApi.sqlite_db_api import PsimiSQL

#declaring constants
SQL_SEED_LOCATION = '../SQLiteDBApi/network-db-seed.sql'
DESTINATION = '../merger'
SOURCE_DB_FILE_LIST = ['../mapper/protein/output/miRDeathDB_mapped.db',]

print("The following dbs will be parsed: " + ", ".join(map(lambda line: line.split("/")[-1], SOURCE_DB_FILE_LIST)))

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

    #converting the strings to lists
    list_1 = string_1.split(separator)
    list_2 = string_2.split(separator)

    list_1 = list(filter(lambda item: item != '-', list_1))
    list_2 = list(filter(lambda item: item != '-', list_2))

    #adding the lists
    sum_lists = list_1 + list_2

    #1) converting the lists to sets 2) removing empty elements 3) getting the union of these sets 4) converting back to list
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

    #merging the two nodes properties that needs to be merged

    merged_pathways = merge_strings(node_dict_1['pathways'], node_dict_2['pathways'])

    merged_aliases = merge_strings(node_dict_1["aliases"], node_dict_2["aliases"])

    alt_acccession = node_dict_1['alt_accession']

    #dealing with the alt accession (gene symbol)

    if node_dict_1['alt_accession'] != node_dict_2['alt_accession']:
        if merged_aliases:
            merged_aliases += "|" + node_dict_2['alt_accession']
        else:
            merged_aliases = node_dict_2['alt_accession']

    #initiating the new node
    new_node = {
        "name" : node_dict_1["name"],
        "alt_accession" : alt_acccession,
        "tax_id" : node_dict_1["tax_id"],
        "pathways" : merged_pathways,
        "aliases" : merged_aliases
    }

    return new_node


def main(log):
    # Declaring the dicts that will hold the data
    nodes = {}
    collected_edges = {}

    merged_edge_counter = 0
    not_merged_edge = 0

    all_counter = 0

    # the number of pieces (.db files)
    sum_files = len(SOURCE_DB_FILE_LIST)

    # filling up the nodes dictionary with the data contained in db piece files
    for db_file in SOURCE_DB_FILE_LIST:

        #executing a query that selects everything (but the node id) from the current SQLite .db files
        db = sqlite3.connect(db_file)
        cursor = db.cursor()
        cursor.execute("SELECT * FROM node WHERE tax_id = 'taxid:9606' OR tax_id = 'taxid:6239' OR tax_id = 'taxid:7227' OR 'taxid:7955'")

        # iterating trough the db row by row
        while True:
            row = cursor.fetchone()
            # until the last row
            if row == None:
                print(all_counter)
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
            if not name in nodes:
                nodes[name] = node
            else:
                nodes[name] = get_union_of_nodes(nodes[name], node)
        # closing the current db
        db.close()
        # logging out some info
        current_file = SOURCE_DB_FILE_LIST.index(db_file)
        sys.stdout.write("Building the node dictionary: Processing %d files out of %d\r" % (current_file, sum_files))

    # making a memory database and inserting the unique nodes from the nodes dictionary
    print('Inserting nodes to database')
    parser = PsimiSQL(SQL_SEED_LOCATION)
    for node in nodes:
        parser.insert_unique_node(nodes[node])
        nodes[node]['id'] = parser.cursor.lastrowid

    # looping through the files again to make an edge list
    print("Started building edge dict")
    file_counter = 1
    for db_file in SOURCE_DB_FILE_LIST:

        sys.stdout.write("Inserting edges to edge dict from '%s' (%d/%d)\r" % (db_file, file_counter, sum_files))

        #executing a query that selects everything (but the node id) from the current SQLite .db files
        db = sqlite3.connect(db_file)
        cursor = db.cursor()
        cursor.execute("SELECT * FROM edge")

        while True:
            row = cursor.fetchone()
            # if there aren't any more nodes break out of the loop
            if not row:
                break
            else:
                # deconstructing the row (list)
                edge_row_id, old_interactor_a_node_id, old_interactor_b_node_id, interactor_a_node_name, interactor_b_node_name, interaction_detection_method , first_author, publication_ids, interaction_types, source_db, interaction_identifiers, confidence_scores, layer = row

                # because in the nodes dict building process the query only asks for human||drosophila||celegans nodes
                # we have to make sure that we don't try to insert edges whose nodes are in the nodes dict (=does not contain any other organisms node)
                if interactor_a_node_name in nodes and interactor_b_node_name in nodes:
                    # generating an edge id that will be the key in the edge dict
                    edge_id = interactor_a_node_name + "@" + interactor_b_node_name

                    # generating an edge dict, that will be a value for the key
                    current_edge = {
                        'interaction_detection_method' : interaction_detection_method,
                        'first_author' : first_author,
                        'publication_ids' : publication_ids,
                        'interaction_types' : interaction_types,
                        'source_db' : source_db,
                        'interaction_identifiers' : interaction_identifiers,
                        'confidence_scores' : confidence_scores,
                        'layer' : layer
                    }

                    # if the edge dict does not have this edge_id. the edge is stored in the edge dict with it's id
                    if not edge_id in collected_edges:
                        collected_edges[edge_id] = []
                        collected_edges[edge_id].append(current_edge)
                    else:
                        # if the edge has this id, the array of edges is searched until the edge with the same interaction type found
                        merged = False
                        for collected_edge in collected_edges[edge_id]:
                            # if an edge is already in the dict with the same interaction type, it will be merged with the current edge
                            # 'already in' means that the two edges have the same effect in their interaction types
                            # TODO: FIX!! only works with layer0, other layers don't have interaction type data
                                collected_edge['interaction_types'] = merge_strings(collected_edge['interaction_types'].split('|')[0], current_edge['interaction_types'].split('|')[0])
                                collected_edge['first_author'] = merge_strings(collected_edge['first_author'],current_edge['first_author'])
                                collected_edge['source_db'] = merge_strings(collected_edge['source_db'], current_edge['source_db'])
                                collected_edge['interaction_identifiers'] = merge_strings(collected_edge['interaction_identifiers'], current_edge['interaction_identifiers'])
                                collected_edge['interaction_detection_method'] = merge_strings(collected_edge['interaction_detection_method'], current_edge['interaction_detection_method'])
                                collected_edge['confidence_scores'] = merge_strings(collected_edge['confidence_scores'], current_edge['confidence_scores'])
                                merged = True
                                merged_edge_counter += 1
                                break
                        if merged == False:
                            # if current edge cannot be merged with an edge that is already in the edges dict, it will be appended to a new edge
                            # TODO: log this!!!
                            collected_edges[edge_id].append(current_edge)
                            not_merged_edge += 1
                            if collected_edge['source_db'] != current_edge['source_db']:
                                pass

    print("Building edge dict done!")
    print("Started inserting edges to the db")

    # iterating through edges dictionary and inserting nodes to the SQLite db
    for k, v in collected_edges.items():
        for edge_to_insert in v:
            # getting the nodes
            node_a, node_b = k.split('@')

            node_a_dict = nodes[node_a]
            node_b_dict = nodes[node_b]

            parser.insert_edge(node_a_dict, node_b_dict, edge_to_insert)

    print("Saving db")
    parser.save_db_to_file(DESTINATION)


if __name__ == '__main__':
    main(log = None)