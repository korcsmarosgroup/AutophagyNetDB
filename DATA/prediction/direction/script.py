"""
Direction prediction based on learning dataset from reactome
PPI direction calculated from domain interaction directions
"""

import sqlite3, csv

REACTOME_DB = '../../workflow/layer0/output/reactome.db'
PFAM_FILE_LIST = ['files/uniprot-pfam_human.tab']
pfam_dict = {}

conn = sqlite3.connect(REACTOME_DB)
with conn:
    c = conn.cursor()
    c.execute("SELECT interactor_a_node_name, interactor_b_node_name FROM edge")
    while True:
        row = c.fetchone()
        if row is None:
            break
        else:
            A_node = row[0].split(':')[1]
            B_node = row[1].split(':')[1]
            for file in PFAM_FILE_LIST:
                with open(file) as infile:
                    infile.readline()
                    for line in infile:
                        line = line.strip().split('\t')
                        if (line[0] == A_node or line[0] == B_node) and len(line) == 4:
                                pfam_dict[line[0]] = line[3]

print(pfam_dict)






