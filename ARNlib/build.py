"""
    This script builds the multi-structured database. Each nodeis present in every layer.
    The higher layer has the edges from the previous layer and their first neighbours.
"""

import sqlite3
from collections import OrderedDict as od

# Building mock data set
DB = od([('signor', 'SLK_Core'),
         ('PSP', 'layer1'),
         ('PhosphoSite', 'PTM'),
         ('OmniPath', 'ATG_Reg'),
         ('miRDeathDB', 'miRNA'),
         ('starbase', 'lncRNA')])


def get_mockdata():
    with open('mockdata.tsv', 'w') as mockfile:
        # Adding header
        mockfile.write(
            'id\tinteractor_node_a_id\tinteractor_node_b_id\tinteractor_node_a_name\tinteractor_node_b_name\tinteraction_detection_method\tfirst_author\tpublication_ids\tinteraction_types\tsource_db\tinteraction_identifier\tconfidence_score\tlayer\n')
        path = 'merger.db'

        conn1 = sqlite3.connect(path)
        c1 = conn1.cursor()
        c1.execute('SELECT * FROM edge')
        counter = 0
        while True:
            edge = c1.fetchone()
            counter += 1
            # Until the last row
            # For just 100 edges from each database
            if edge is None:
                break
            mockfile.write("%s\n" % '\t'.join(map(str, edge)))


# Calling the function
get_mockdata()

# Creating output table
conn = sqlite3.connect('SLK3_layers.db')
with conn:
    c = conn.cursor()
    c.execute("DROP TABLE IF EXISTS layer0")
    c.execute("DROP TABLE IF EXISTS layer1")
    c.execute("DROP TABLE IF EXISTS layer2")
    c.execute("DROP TABLE IF EXISTS layer3")
    c.execute("DROP TABLE IF EXISTS layer5")
    c.execute("DROP TABLE IF EXISTS layer6")
    c.execute("DROP TABLE IF EXISTS layer7")
    c.execute('''CREATE TABLE layer0 (id INT,
                                    interactor_a_node_id INT,
                                    interactor_b_node_id INT,
                                    interactor_a_node_name TEXT,
                                    interactor_b_node_name TEXT,
                                    interaction_detection_method TEXT,
                                    first_author TEXT,
                                    publication_ids TEXT,
                                    interaction_identifiers TEXT,
                                    source_db TEXT,
                                    interaction_types TEXT,
                                    confidence_scores TEXT,
                                    layer TEXT)''')
    c.execute('''CREATE TABLE layer1 (id INT,
                                    interactor_a_node_id INT,
                                    interactor_b_node_id INT,
                                    interactor_a_node_name TEXT,
                                    interactor_b_node_name TEXT,
                                    interaction_detection_method TEXT,
                                    first_author TEXT,
                                    publication_ids TEXT,
                                    interaction_identifiers TEXT,
                                    source_db TEXT,
                                    interaction_types TEXT,
                                    confidence_scores TEXT,
                                    layer TEXT)''')
    c.execute('''CREATE TABLE layer2 (id INT ,
                                    interactor_a_node_id INT,
                                    interactor_b_node_id INT,
                                    interactor_a_node_name TEXT,
                                    interactor_b_node_name TEXT,
                                    interaction_detection_method TEXT,
                                    first_author TEXT,
                                    publication_ids TEXT,
                                    interaction_identifiers TEXT,
                                    source_db TEXT,
                                    interaction_types TEXT,
                                    confidence_scores TEXT,
                                    layer TEXT)''')
    c.execute('''CREATE TABLE layer3 (id INT ,
                                    interactor_a_node_id INT,
                                    interactor_b_node_id INT,
                                    interactor_a_node_name TEXT,
                                    interactor_b_node_name TEXT,
                                    interaction_detection_method TEXT,
                                    first_author TEXT,
                                    publication_ids TEXT,
                                    interaction_identifiers TEXT,
                                    source_db TEXT,
                                    interaction_types TEXT,
                                    confidence_scores TEXT,
                                    layer TEXT)''')
    c.execute('''CREATE TABLE layer5 (id INT ,
                                    interactor_a_node_id INT,
                                    interactor_b_node_id INT,
                                    interactor_a_node_name TEXT,
                                    interactor_b_node_name TEXT,
                                    interaction_detection_method TEXT,
                                    first_author TEXT,
                                    publication_ids TEXT,
                                    interaction_identifiers TEXT,
                                    source_db TEXT,
                                    interaction_types TEXT,
                                    confidence_scores TEXT,
                                    layer TEXT)''')
    c.execute('''CREATE TABLE layer6 (id INT ,
                                    interactor_a_node_id INT,
                                    interactor_b_node_id INT,
                                    interactor_a_node_name TEXT,
                                    interactor_b_node_name TEXT,
                                    interaction_detection_method TEXT,
                                    first_author TEXT,
                                    publication_ids TEXT,
                                    interaction_identifiers TEXT,
                                    source_db TEXT,
                                    interaction_types TEXT,
                                    confidence_scores TEXT,
                                    layer TEXT)''')
    c.execute('''CREATE TABLE layer7 (id INT ,
                                    interactor_a_node_id INT,
                                    interactor_b_node_id INT,
                                    interactor_a_node_name TEXT,
                                    interactor_b_node_name TEXT,
                                    interaction_detection_method TEXT,
                                    first_author TEXT,
                                    publication_ids TEXT,
                                    source_db TEXT,
                                    interaction_identifiers TEXT,
                                    interaction_types TEXT,
                                    confidence_scores TEXT,
                                    layer TEXT)''')

    # Extracting the data
    with open('mockdata.tsv') as mockfile:
        l0_nodes = []
        l1_nodes = []
        l2_nodes = []
        l3_nodes = []
        l5_nodes = []
        l6_nodes = []
        l7_nodes = []
        # Skipping the header
        mockfile.readline()
        for line in mockfile:
            line = line.strip().split('\t')
            c = conn.cursor()
            # Extracting L0 nodes
            if line[-1] == '0':
                l0_nodes.append(line[3])
                l0_nodes.append(line[4])
                c.execute("INSERT INTO layer0 VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)",
                          line)
            # Checking if there are edges with previous layer nodes in the layer above
            # For each layer
            if line[-1] == '1':
                if line[3] in l0_nodes or line[4] in l0_nodes:
                    l1_nodes.append(line[3])
                    l1_nodes.append(line[4])
                    c.execute("INSERT INTO layer1 VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)",
                              line)
            if line[-1] == '2':
                if line[3] in l1_nodes or line[4] in l1_nodes:
                    l2_nodes.append(line[3])
                    l2_nodes.append(line[4])
                    c.execute("INSERT INTO layer2 VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)",
                              line)
            if line[-1] == '3':
                if line[3] in l2_nodes or line[4] in l2_nodes:
                    l3_nodes.append(line[3])
                    l3_nodes.append(line[4])
                    c.execute("INSERT INTO layer3 VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)",
                              line)
            if line[-1] == '5':
                if line[4] in l3_nodes or line[3] in l3_nodes:
                    l5_nodes.append(line[3])
                    l5_nodes.append(line[4])
                    c.execute("INSERT INTO layer5 VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)",
                              line)
            if line[-1] == '6':
                if line[4] in l5_nodes or line[3] in l5_nodes:
                    l6_nodes.append(line[3])
                    l6_nodes.append(line[4])
                    c.execute("INSERT INTO layer6 VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)",
                              line)
            if line[-1] == '7':
                if line[4] in l6_nodes or line[3] in l6_nodes:
                    l7_nodes.append(line[3])
                    l7_nodes.append(line[4])
                    c.execute("INSERT INTO layer7 VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)",
                              line)
    # Adding additional edge data based on layer0
    #  Use another cursor if using UPDATE statement
    c2 = conn.cursor()
    num = 0
    while num < 7:
        # Getting data from lower layer
        c.execute('''SELECT interactor_a_node_name, interactor_b_node_name, interaction_identifiers
                          FROM layer%d''' % num)
        while True:
            row = c.fetchone()
            if row is None:
                if num == 3:
                    num = 4
                num += 1
                break
            else:
                # Only need directedness information from interaction_identifiers of layer0
                if num == 0:
                    direct = row[2].split('|')[1]
                else:
                    direct = row[2]
                if num+1 != 4:
                    # Updating interaction_identifiers columns, based on the previous layer, with adding already existing data
                    # SQLite string concatenation: ||
                    # TODO: fix multiple adding of interaction_identifiers
                    c2.execute('''UPDATE layer%d SET interaction_identifiers = '%s' || '|' || layer%d.interaction_identifiers
                                          WHERE interactor_a_node_name = '%s' OR interactor_b_node_name = '%s'
                                          OR interactor_a_node_name = '%s' OR interactor_b_node_name = "%s"'''
                               % (num+1, direct, num+1, row[0], row[0], row[1], row[1]))
