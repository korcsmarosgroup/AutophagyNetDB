'''
This script translates all the node names or IDs to miRbase identifiers
'''

import csv
import sqlite3



def create_miRBase_map(all_id):
    # Creating map dataset
    conn = sqlite3.connect('miRNAmap.db')
    with conn:
        c = conn.cursor()
        c.execute('''DROP TABLE if EXISTS map''')
        c.execute('''CREATE TABLE map (id INTEGER PRIMARY KEY AUTOINCREMENT,
                                       mirName CHAR (100),
                                       mirAC CHAR (100))''')
    # Inserting data
    with open(all_id) as infile:
        infile.readline()
        with conn:
            c = conn.cursor()
            insert_name = []
            insert_id = []
            for line in infile:
                line = line.split(',')
                name = line[1]
                name2 = line[5]
                name3 = line[8]
                insert_name.append(str(name))
                insert_name.append(str(name2))
                insert_name.append(str(name3))
                id = line[0]
                id2 = line[4]
                id3 = line[7]
                insert_id.append(str(id))
                insert_id.append(str(id2))
                insert_id.append(str(id3))

            for i in range(len(insert_id)):
                c.execute('''INSERT INTO map (mirName, mirAC) VALUES (?,?)''',
                          (insert_name[i], insert_id[i]))


#create_miRBase_map(all_id='miRNA.csv')
