"""
This script counts how mane edges are in a certain pathway
"""
__author__ = 'blaise'

import sqlite3
import json
import pprint

#queryiing the database
query = """
SELECT edge.interactor_a_node_name, edge.interactor_b_node_name, A.pathways as a_pathways, B.pathways as b_pathways
FROM edge
JOIN node AS A ON edge.interactor_a_node_id = A.id
JOIN node AS B ON edge.interactor_b_node_id = B.id
"""

db = sqlite3.connect("reactome.db")

cursor = db.cursor()

cursor.execute(query)

pathway_counters = {}

while True:
    row = cursor.fetchone()
    if not row:
        break
    else:
        source_name, target_name, source_pathways, target_pathways = row
        source_pathways_list = source_pathways.split('|')
        target_pathways_list = target_pathways.split('|')
        edge_pathways = source_pathways_list + target_pathways_list
        edge_pathway_set = set(edge_pathways)

        for pathway in edge_pathway_set:
            if "signalink" in pathway:
                if pathway_counters.has_key(pathway):
                    pathway_counters[pathway] += 1
                else:
                    pathway_counters[pathway] = 1

file_content = json.dumps(pathway_counters)
file_content = pprint.pformat(file_content)

with open("result", "w") as _file:
    _file.write(file_content)