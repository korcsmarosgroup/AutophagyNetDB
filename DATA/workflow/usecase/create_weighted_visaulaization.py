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
node_annot = node_annot[node_annot['selected'] == True]

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
node_to_ptw_dict = node_annot[['name', 'pathways']].set_index('name').T.to_dict('list')

path_node_dict = {}
node_size_dict = {}
for path in accepted_pathways:
    path_node_dict[path] = node_annot[node_annot['pathways'].astype("string").str.contains(path)].name.tolist()
    ptw_count = len(node_annot[node_annot['pathways'].astype("string").str.contains(path)])
    node_size_dict[path] = ptw_count

# Get l1 nodes
l1_nodes = pd.read_csv('xeno_L1_nodes.csv')
l1_nodes = l1_nodes[l1_nodes['selected'] == True]
l1_nodes['new_name'] = 'L1'
L1_dict = dict(zip(l1_nodes.name, l1_nodes.new_name))

# Read in interaction file of L1 and L3 interactions
int_df = pd.read_csv('ptw_conn_int.csv')
int_df = int_df.drop('interaction', axis=1)
int_df = int_df.drop('selected', axis=1)
int_df = int_df.drop('shared interaction', axis=1)
int_df = int_df.drop('shared name', axis=1)

int_df['source'] = int_df['name'].astype("string").str.split(' ', expand=True)[0]
int_df['target'] = int_df['name'].astype("string").str.split(' ', expand=True)[3]

int_df.drop('name', inplace = True, axis=1)


# Count number of interactions for each pathway between core and l3 (if node is source)
# TODO: count number of interactions of a pathway with each core node (replace source node with pathway name - import this into cytoscape)

# Rename node df so column names match int_df
node_annot.rename({'name': 'source'}, axis=1, inplace=True)
# If source node names match add pathway column from node df to int df
merged_df = int_df.merge(node_annot[['source', 'pathways']], how='left')
# Change the source name to pathways
merged_df.drop('source', axis=1, inplace=True)
merged_df.rename({'pathways': 'source'}, axis=1, inplace=True)
# If multiple pathways for a node, enter it as multiple separate rows for each pathway
merged_df['source'] = merged_df['source'].str.split('|')
merged_df = merged_df.explode('source')

# get dict and df of which ptw each L1 node is connecting with
ptw_to_l1_dict = merged_df.groupby('source')['target'].apply(list).to_dict()
ptw_to_l1 = pd.DataFrame(merged_df.groupby(['source', 'target']).size().reset_index(name = "Group_Count"))

# Group L1 nodes
df2=merged_df.replace({"target": L1_dict})
merged_df['target'] = merged_df['target'].map(L1_dict).fillna(merged_df['target'])

# Group by target node and count number of interactions per node for each pathway
ptw_count = df2.groupby(['source', 'target']).size()

# print(ptw_count)

core_nodes = ['Uniprot:Q9BXW4', 'Uniprot:Q8TDY2', 'Uniprot:Q9H0R8', 'Uniprot:O95166', 'Uniprot:Q7Z6L1', 'Uniprot:Q9H1Y0', 'Uniprot:Q8IYT8', 'Uniprot:Q9NT62', 'Uniprot:Q14457',
              'Uniprot:Q9GZQ8', 'Uniprot:O94817', 'Uniprot:Q9Y4P8', 'Uniprot:Q676U5', 'Uniprot:P60520']
# Get L1-core interaction counts
# Rename node df so column names match int_df
l1_nodes.rename({'name': 'target'}, axis=1, inplace=True)

l1_nodes = l1_nodes.merge(ptw_to_l1[['target', 'source']], how='left')
l1_nodes.rename({'source': 'pathway'}, axis=1, inplace=True)
l1_nodes.rename({'target': 'source'}, axis=1, inplace=True)

# If source node names match add pathway column from node df to int df
l1_merged_df = int_df.merge(l1_nodes[['source', 'pathway']], how='left')
l1_merged_df = l1_merged_df[l1_merged_df['target'].isin(core_nodes)]

# Change the source name to pathways
l1_merged_df.drop('source', axis=1, inplace=True)
l1_merged_df.rename({'pathway': 'source'}, axis=1, inplace=True)

# Group by target node and count number of interactions per node for each pathway
l1_ptw_count = l1_merged_df.groupby(['source', 'target']).size()
print(l1_ptw_count)

# Merge two dataframes adding edge weights of Ptw-core and L1-core together
whole_df = pd.concat([ptw_count, l1_ptw_count]).groupby(['source', 'target']).sum().reset_index()
print(whole_df)

# Write to file
ptw_count.reset_index().to_csv('xeno_ptw_l1_edge_width.csv', sep=',', header=True, index=False)
l1_ptw_count.reset_index().to_csv('xenol1_core_edge_width.csv', sep=',', header=True, index=False)
whole_df.reset_index().to_csv('xeno_final_edge_width.csv', sep=',', header=True, index=False)

# Add node counts
for pathway in ptw_to_l1_dict.keys():
    node_size_dict[pathway] += len(ptw_to_l1_dict[pathway])

print(node_size_dict)
# Write to file
with open('xeno_node_size.csv', 'w') as f:  # You will need 'wb' mode in Python 2.x
    for key in node_size_dict.keys():
        f.write("%s, %s\n" % (key, node_size_dict[key]))
