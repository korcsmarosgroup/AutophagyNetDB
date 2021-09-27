"""
 Checks for duplicates in final json files
"""

# Imports
import json


def is_duplicate_edge(filename):
    edge_id_list = []
    duplicate_list = []
    jsonfile = json.load(open(filename))
    for dict in jsonfile:
        edge_id = dict['source'] + '@' + dict['target']
        if edge_id not in edge_id_list:
            edge_id_list.append(edge_id)
        else:
            duplicate_list.append(edge_id)

    return duplicate_list


def is_duplicate_node(filename):
    node_id_list = []
    duplicate_list = []
    jsonfile = json.load(open(filename))
    for dict in jsonfile:
        node_id = dict['name']
        if node_id not in node_id_list:
            node_id_list.append(node_id)
        else:
            duplicate_list.append(node_id)

    print(duplicate_list)


is_duplicate_node('nodes.json')
# is_duplicate_edge('edges.json')


