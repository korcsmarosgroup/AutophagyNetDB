"""
    Collapse layer3 into 14 pathways - each as a node - node size: how many interactors are in the pathway
        If a node is annotated to multiple pathways - count it in each pathway
    Count number of edges for each core node - where core node is the target and pathway node is the source - sum edge number for each pathway
        Edge width will correspond to number of connections
"""

import pandas as pd
import csv

# Xenophagy

# Read in node file
node_annot = pd.read_csv('xeno_ptwconn_L3_node.csv')

accepted_pathways = [
    'T-cell receptor',
    'B-cell receptor',
    'Toll-like receptor',
    'Innate immune pathways',
    'JAK/STAT',
    'Nuclear hormone receptor',
    'Receptor tyrosine kinase',
    'Rho pathway',
    'TGF',
    'Notch',
    'G-protein coupled receptor',
    'WNT/Wingless',
    'Hippo',
    'Hedgehog',
    'TNF pathway',
	'RTK',
	'TCR',
	'NHR'
]

path_node_dict = {}
node_size_dict = {}
for path in accepted_pathways:
    path_node_dict[path] = node_annot[node_annot['pathways'].astype("string").str.contains(path)].name.tolist()
    ptw_count = len(node_annot[node_annot['pathways'].astype("string").str.contains(path)])
    node_size_dict[path] = ptw_count

# Read in interaction file
int_df = pd.read_csv('ptw_conn_int.csv')

int_df = int_df.drop('interaction', axis=1)
int_df = int_df.drop('selected', axis=1)
int_df = int_df.drop('shared interaction', axis=1)
int_df = int_df.drop('shared name', axis=1)

int_df['source'] = int_df['name'].astype("string").str.split(' ', expand=True)[0]
int_df['target'] = int_df['name'].astype("string").str.split(' ', expand=True)[3]

int_df.drop('name', inplace = True, axis=1)


# Count number of interactions for each pathway between core and l3 (if node is source)
# TODO: count number of intreactions of a pathway with each core node (replace source node with pathway name - import this into cytoscape)
edge_count_dict = {}
for ptw, nodes in path_node_dict.items():
    int_num = 0
    for prot in nodes:
        int_num += len(int_df[int_df['source'].astype("string").str.contains(prot)])
    edge_count_dict[ptw] = int_num

# Write to file
with open('xeno_edge_width.csv', 'w') as f:  # You will need 'wb' mode in Python 2.x
    for key in edge_count_dict.keys():
        f.write("%s, %s\n" % (key, edge_count_dict[key]))

with open('xeno_node_size.csv', 'w') as f:  # You will need 'wb' mode in Python 2.x
    for key in node_size_dict.keys():
        f.write("%s, %s\n" % (key, node_size_dict[key]))
