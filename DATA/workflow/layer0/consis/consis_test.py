'''
This script checks whether interaction data from multiple outside sources are correct. Not directed.
'''

import sqlite3, sys
from collections import Counter
import csv

# Defining constants
DB_LIST = ["../mapper/protein/output/acsn_mapped.db", "../mapper/protein/output/innatedb_mapped.db", "../mapper/protein/output/reactome_mapped.db", "../mapper/protein/output/signor_mapped.db"]
merger = '../merger.db'
'''
dbconn = sqlite3.connect('../mapper/protein/output/signor_mapped.db')
a_dict = {}
b_dict = {}
int_dict = {}
with dbconn:
    c = dbconn.cursor()
    c.execute("SELECT * FROM edge")
    while True:
        row = c.fetchone()
        if row is None:
            break
        # Inserting each A node and what nodes it interacts with into a dictionary
        else:
            if not row[3] in a_dict:
                a_dict[row[3]] = [row[4]]
            else:
                if row[4] in a_dict[row[3]]:
                    #print(row[3], row[4], row[8])
                    if not row[3] + '-' + row[4] in int_dict:
                        int_dict[row[3] + '-' + row[4]] = [row[8]]
                    else:
                        int_dict[row[3] + '-' + row[4]].append(row[8])
                a_dict[row[3]].append(row[4])

            # Inserting each B node and what it interacts with
            if not row[4] in b_dict:
                b_dict[row[4]] = [row[3]]
            else:
                b_dict[row[4]].append(row[3])

    #for key in int_dict:
     #   if all(val == int_dict[key][0] for val in int_dict[key]) is False:
      #          print(key.split('-'), Counter(int_dict[key]))


    #print(int_dict)

    for key in b_dict:
        for key2 in a_dict:
            if key == key2:
                for k, v in Counter(a_dict[key]).items():
                    if v > 1:
                        #print(key, k)
                        c.execute("SELECT interactor_a_node_name, interactor_b_node_name, interaction_types, source_db FROM EDGE WHERE interactor_a_node_name = '%s' AND interactor_b_node_name = '%s'" % (
                                key, k))
                        same_row = c.fetchone()
                        #if same_row is not None:
                            #print(same_row)

                for val in a_dict[key]:
                    if val in b_dict[key]:

                        c.execute(
                            "SELECT interactor_a_node_name, interactor_b_node_name, interaction_types, source_db FROM EDGE WHERE interactor_a_node_name = '%s' AND interactor_b_node_name = '%s'" % (
                                key2, val))
                        a_row = c.fetchone()
                        #print(a_row)
                        c.execute(
                            "SELECT interactor_a_node_name, interactor_b_node_name, interaction_types, source_db FROM EDGE WHERE interactor_a_node_name = '%s' AND interactor_b_node_name = '%s'" % (
                                val, key))
                        b_row = c.fetchone()
                        if a_row is not None and b_row is not None:
                            if a_row[2] != b_row[2]:
                                print(a_row, b_row)
'''
conn = sqlite3.connect(merger)

with open('non_supporting.tsv', 'w') as tsvfile1, open('supporting.tsv', 'w') as tsvfile2:
    writer = csv.writer(tsvfile1, delimiter='\t')
    writer2 = csv.writer(tsvfile2, delimiter='\t')
    writer.writerow(['A_node_name', 'B_node_name', 'Int_type', 'Source_db'])
    writer2.writerow(['A_node_name', 'B_node_name', 'Int_type', 'Source_db'])
    db_list = ['acsn', 'innatedb', 'signor', 'reactome']
    # Getting edge data from specific database
    for db in db_list:
        db_dict = {}
        db_b_dict = {}
        int_dict = {}
        with conn:
            c = conn.cursor()
            c.execute("SELECT * FROM EDGE WHERE source_db = 'source database:%s'" % db)
            while True:
                edge_row = c.fetchone()
                #print(edge_row)
                if edge_row is None:
                    break
                # Inserting each A node and what nodes it interacts with into a dictionary
                else:
                    if not edge_row[3] in db_dict:
                        db_dict[edge_row[3]] = [edge_row[4]]
                    else:
                        if not edge_row[3] + '-' + edge_row[4] in int_dict:
                            int_dict[edge_row[3] + '-' + edge_row[4]] = [edge_row[8]]
                        else:
                            int_dict[edge_row[3] + '-' + edge_row[4]].append(edge_row[8])
                        db_dict[edge_row[3]].append(edge_row[4])

                    # Inserting each B node and what it interacts with
                    #if not edge_row[4] in db_b_dict:
                        #db_b_dict[edge_row[4]] = [edge_row[3]]
                    #else:
                        #db_b_dict[edge_row[4]].append(edge_row[3])
            #print(int_dict)
            for key in int_dict:
                if all(val == int_dict[key][0] for val in int_dict[key]) is False:
                        for inter in int_dict[key]:
                            writer.writerow([key.split('-')[0], key.split('-')[1], inter, db])
                else:
                    for inter in int_dict[key]:
                        writer2.writerow([key.split('-')[0], key.split('-')[1], inter, db])

'''
    # Checking if there are edges represented multiple times in the database
    for key in db_b_dict:
        for key2 in db_dict:
            if key == key2:
                #print(key, db_dict[key], db_b_dict[key])
                for val in db_dict[key]:
                    if val in db_b_dict[key]:
                        c.execute(
                            "SELECT interactor_a_node_name, interactor_b_node_name, interaction_types, source_db FROM EDGE WHERE source_db = 'source database:%s' AND interactor_a_node_name = '%s' AND interactor_b_node_name = '%s'" % (
                            db, key2, val))
                        a_row = c.fetchone()
                        c.execute(
                            "SELECT interactor_a_node_name, interactor_b_node_name, interaction_types, source_db FROM EDGE WHERE source_db = 'source database:%s' AND interactor_a_node_name = '%s' AND interactor_b_node_name = '%s'" % (
                                db, val, key))
                        b_row = c.fetchone()
                        #if a_row is not None and b_row is not None:
                         #   if a_row[2] != b_row[2]:
                          #      print(a_row, b_row)

    for key in db_dict:
        seen = set()
        uniq = []
        same = []
        for val in db_dict[key]:
            if val not in seen:
               uniq.append(val)
               seen.add(val)
            else:
                same.append(val)
        if len(same) != 0:
            #print(key, same)
            for i in range(len(same)):
                c.execute(
                    "SELECT interactor_a_node_name, interactor_b_node_name, interaction_types, source_db FROM EDGE WHERE source_db = 'source database:%s' AND interactor_a_node_name = '%s' AND interactor_b_node_name = '%s'" % (
                        db, key, same[i]))
                #print(c.fetchone())
'''

with open('non_supp_pair.tsv', 'w') as tsvnon, open('supp_pair.tsv', 'w') as tsvsupp:
    writernon = csv.writer(tsvnon, delimiter='\t')
    writersupp = csv.writer(tsvsupp, delimiter='\t')
    writernon.writerow(['A_node_name', 'B_node_name', 'Source_db_1', 'Int_type_1', 'Source_db_2', 'Int_type_2'])
    writersupp.writerow(['A_node_name', 'B_node_name', 'Source_db_1', 'Int_type_1', 'Source_db_2', 'Int_type_2'])

    for db in db_list:
        db1 = db
        for other_db in db_list:
            db2 = other_db
            if db1 != db2:
                db_dict1 = {}
                with conn:
                    c = conn.cursor()
                    c.execute("SELECT * FROM EDGE WHERE source_db = 'source database:%s'" % db1)
                    while True:
                        edge_row = c.fetchone()
                        if edge_row is None:
                            break
                        else:
                            if not edge_row[3] in db_dict1:
                                db_dict1[edge_row[3]] = [edge_row[4]]
                            else:
                                db_dict1[edge_row[3]].append(edge_row[4])

                db_dict2 = {}
                with conn:
                    c = conn.cursor()
                    c.execute("SELECT * FROM EDGE WHERE source_db = 'source database:%s'" % db2)
                    while True:
                        edge_row = c.fetchone()
                        if edge_row is None:
                            break
                        else:
                            if not edge_row[3] in db_dict2:
                                db_dict2[edge_row[3]] = [edge_row[4]]
                            else:
                                db_dict2[edge_row[3]].append(edge_row[4])


            # Checking if node As from different databases interact with the same B nodes
                for key in list(db_dict1.keys()):
                    if key in list(db_dict2.keys()):
                        #print(key, db_dict1[key], db_dict2[key])

                        for val in db_dict1[key]:
                            if val in db_dict2[key]:
                                c.execute("SELECT interactor_a_node_name, interactor_b_node_name, interaction_types, source_db FROM EDGE WHERE interactor_a_node_name = '%s' AND interactor_b_node_name = '%s'" % (
                                    key, val))
                                edge_row = c.fetchall()
                                #print(edge_row)
                                if edge_row[0][2] != edge_row[1][2]:
                                    writernon.writerow([edge_row[0][0], edge_row[0][1], edge_row[0][3], edge_row[0][2], edge_row[1][3], edge_row[1][2]])
                                else:
                                    writersupp.writerow([edge_row[0][0], edge_row[0][1], edge_row[0][3], edge_row[0][2], edge_row[1][3], edge_row[1][2]])

                        for val2 in db_dict2[key]:
                            if val2 in db_dict1[key]:
                                c.execute(
                                    "SELECT interactor_a_node_name, interactor_b_node_name, interaction_types, source_db FROM EDGE WHERE interactor_a_node_name = '%s' AND interactor_b_node_name = '%s'" % (
                                    key, val2))
                                edge_row2 = c.fetchall()
                                #print(edge_row2)
                                if edge_row2 != edge_row:
                                    #print(edge_row2)
                                    if edge_row2[0][2] != edge_row2[1][2]:
                                        writernon.writerow([edge_row[0][0], edge_row[0][1], edge_row[0][3], edge_row[0][2], edge_row[1][3], edge_row[1][2]])
                                    else:
                                        writersupp.writerow([edge_row[0][0], edge_row[0][1], edge_row[0][3], edge_row[0][2], edge_row[1][3], edge_row[1][2]])





