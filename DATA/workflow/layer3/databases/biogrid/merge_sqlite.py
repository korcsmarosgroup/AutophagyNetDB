"""
This script merges SQLite .db files with the same schema
:argument sys.argv[1]: Path to PsimiSQL class
"""

import os
import sys
import sqlite3

#adding the psimi_to_sql module to sys
sys.path.append(sys.argv[1])
from  sqlite_db_api import PsimiSQL

#declaring CONSTANTS
#getting the list of files from the current folder that will be merged
PIECE_LIST = filter(lambda filename: 'db_piece' in filename, os.listdir('./'))
SQL_SEED = 'C:/Users/Luca/PycharmProjects/SLK3/DATA/workflow/lib/SQLiteDBApi/network-db-seed.sql'

#declaring the dicts that will hold the data
nodes = {}

#the number of pieces (.db files)
sum_files = len(PIECE_LIST)

#filling up the nodes dictionary with the data contained in db piece files
for filename in PIECE_LIST:

    #executing a query that selects everything (but the node id) from the current SQLite .db files
    db = sqlite3.connect(filename)
    cursor = db.cursor()
    cursor.execute("SELECT name, alt_accession, tax_id from NODE")

    #iterating trough the db row by row
    while True:
        row = cursor.fetchone()
        #until the last row
        if row == None:
            break
        #if unique, inserting the node (row) to the nodes dictionary
        name, alt_accession, tax_id = row
        node = {
            "name" : name,
            'alt_accession' : alt_accession,
            'tax_id' : tax_id,
            'pathways' : None
        }
        if not nodes.has_key(name):
            nodes[name] = node
        else:
            nodes[name].update(node)
    #closing the current db
    db.close()
    #logging out some info
    current_file = PIECE_LIST.index(filename)
    sys.stdout.write("Building the node dictionary: Processing %d files out of %d\r" % (current_file, sum_files))

#making a memory database and inserting the unique nodes from the nodes dictionary
print('Inserting nodes to database')
parser = PsimiSQL(SQL_SEED)
for node in nodes:
    parser.insert_unique_node(nodes[node])

#now that we have the nodes in the final db, the edges can be inserted
#there is no need for a edges dictionary, because reading it from the files costs less memory
#iterating through the .db piece files again
print('Inserting edges to database')
for filename in PIECE_LIST:
    db = sqlite3.connect(filename)
    query = "SELECT * FROM edge"
    cursor = db.cursor()
    cursor.execute(query)

    #iterating trough the current piece .db files
    while True:
        edge_row = cursor.fetchone()
        if edge_row == None:
            break
        edge_dict = {
            'interaction_detection_method' : edge_row[5],
            'first_author' : edge_row[6],
            'publication_ids' : edge_row[7],
            'interaction_types' : edge_row[8],
            'source_db' : edge_row[9],
            'interaction_identifiers' : edge_row[10],
            'confidence_scores' : edge_row[11],
            'layer' : "3"
        }
        parser.insert_edge(nodes[edge_row[3]],nodes[edge_row[4]],edge_dict)


print('Saving database as biogrid_merged.db')
parser.save_db_to_file('biogrid_merged')
pass