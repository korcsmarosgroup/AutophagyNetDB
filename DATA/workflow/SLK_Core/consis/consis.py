"""
This script checks whether interaction data from multiple outside sources contain the same information
And if there are non supporting entries within a database
"""

# Imports
import sqlite3
from collections import Counter
import csv

# Defining constants
db_list = ['acsn',
           'innatedb',
           'signor',
           'reactome',
           'slk2',
           'slk21',
           'slk3',
           'tcr']

# Connecting to the database
merger = '../../merger.db'
conn = sqlite3.connect(merger)

def single_db_compare():
    """
    Checks if there are duplicate interactions with different information
    :return: 2 .tsv files
    """

    rank_dict = {
        'signor': 8,
        'slk3': 7,
        'slk21': 6,
        'slk2': 5,
        'innatedb': 4,
        'reactome': 3,
        'tcr': 2,
        'acsn': 1
    }

    # Creating .tsv files for supporting and non-supporting edges
    with open('non_supporting.tsv', 'w') as tsvfile1, open('supporting.tsv', 'w') as tsvfile2:
        writernon = csv.writer(tsvfile1, delimiter='\t')
        writersupp = csv.writer(tsvfile2, delimiter='\t')
        writernon.writerow(['A_node_name',
                            'B_node_name',
                            'Int_type_1',
                            'Source_db'])

        writersupp.writerow(['A_node_name',
                             'B_node_name',
                             'Int_type',
                             'Source_db'])

        # Getting edge data from specific database
        for db in db_list:
            db_dict = {}
            db_b_dict = {}
            int_dict = {}
            int_b_dict = {}
            with conn:
                c = conn.cursor()
                c.execute("SELECT * FROM EDGE")
                # Until last row
                while True:
                    edge_row = c.fetchone()
                    if edge_row is None:
                        break
                    # Inserting each A node and what nodes it interacts with into a dictionary
                    else:
                        for i in range(len(edge_row[9].split("|"))):
                            if edge_row[9].split("|")[i] == 'source database:' + db:

                                if not edge_row[3] in db_dict:
                                    db_dict[edge_row[3]] = [edge_row[4]]
                                else:
                                    # Adding each interaction to a dictionary with corresponding interaction type
                                    if not edge_row[3] + '-' + edge_row[4] in int_dict:
                                        if edge_row[8].split(':')[-1] != edge_row[8].split(':')[2]:
                                            int_dict[edge_row[3] + '-' + edge_row[4]] = [edge_row[8].split(':')[-1], edge_row[8].split(':')[2].split('|')[0]]
                                        else:
                                            int_dict[edge_row[3] + '-' + edge_row[4]] = [edge_row[8].split(':')[-1]]
                                    else:
                                        int_dict[edge_row[3] + '-' + edge_row[4]].append(edge_row[8].split(':')[-1])
                                        if edge_row[8].split(':')[-1] != edge_row[8].split(':')[2]:
                                                int_dict[edge_row[3] + '-' + edge_row[4]].append(edge_row[8].split(':')[2].split('|')[0])
                                    db_dict[edge_row[3]].append(edge_row[4])

                                # Inserting each B node and what it interacts with
                                if not edge_row[4] in db_b_dict:
                                    db_b_dict[edge_row[4]] = [edge_row[3]]
                                else:
                                    if not edge_row[4] + '-' + edge_row[3] in int_b_dict:
                                        if edge_row[8].split(':')[-1] != edge_row[8].split(':')[2]:
                                            int_b_dict[edge_row[4] + '-' + edge_row[3]] = [edge_row[8].split(':')[-1], edge_row[8].split(':')[2].split('|')[0]]
                                        else:
                                            int_dict[edge_row[3] + '-' + edge_row[4]] = [edge_row[8].split(':')[-1]]
                                    else:
                                        int_b_dict[edge_row[4] + '-' + edge_row[3]].append(edge_row[8].split(':')[-1])
                                        if edge_row[8].split(':')[-1] != edge_row[8].split(':')[2]:
                                                int_b_dict[edge_row[4] + '-' + edge_row[3]].append(edge_row[8].split(':')[2].split('|')[0])
                                    db_b_dict[edge_row[4]].append(edge_row[3])

                nonlist = []
                supplist = []

                for key in int_dict:
                    # If the same edge is present with different int type data, inserting the edge to .tsv
                    if all(val == int_dict[key][0] for val in int_dict[key]) is False:
                        for int_type in int_dict[key]:
                            nonlist.append([key.split('-')[0], key.split('-')[1], int_type, db])
                            writernon.writerow(
                                [key.split('-')[0],
                                 key.split('-')[1],
                                 int_type, db])
                            c.execute(
                            "DELETE FROM EDGE WHERE interactor_a_node_name = 'uniprot:%s' AND interactor_b_node_name = 'uniprot:%s' AND source_db = 'source database:%s'" % (
                                key.split('-')[0].split(':')[1],
                                key.split('-')[1].split(':')[1],
                                db))

                    # Otherwise, there are no duplicate interactions
                    else:
                        for inter in int_dict[key]:
                            supplist.append(
                                [key.split('-')[0],
                                 key.split('-')[1],
                                 inter,
                                 db])

                            writersupp.writerow(
                                [key.split('-')[0],
                                 key.split('-')[1],
                                 inter,
                                 db])

                for key_b in int_b_dict:
                    # If the same edge is present with different int type data, inserting the edge to .tsv
                    if all(val == int_b_dict[key_b][0] for val in int_b_dict[key_b]) is False:
                        for int_b_type in int_b_dict[key_b]:
                            if [key_b.split('-')[1],
                                key_b.split('-')[0],
                                int_b_type, db] not in nonlist:
                                writernon.writerow(
                                    [key_b.split('-')[1],
                                     key_b.split('-')[0],
                                     int_b_type,
                                     db,
                                     2])
                                c.execute(
                                "DELETE FROM EDGE WHERE interactor_a_node_name = 'uniprot:%s' AND interactor_b_node_name = 'uniprot:%s' AND source_db = 'source database:%s'" % (
                                    key_b.split('-')[1].split(':')[1],
                                    key_b.split('-')[0].split(':')[1],
                                    db))

                    # Otherwise, there are no duplicate interactions
                    else:
                        for inter_b in int_b_dict[key_b]:
                            if [key_b.split('-')[1], key_b.split('-')[0], inter_b, db] not in supplist:
                                writersupp.writerow(
                                    [key_b.split('-')[1],
                                     key_b.split('-')[0],
                                     inter_b,
                                     db,
                                     2])


def pair_db_compare():
    """
    Checks if an edge is present in another database, checks if interaction data is the same
    :return: 2 .tsv files
    """

    rank_dict = {
        'signor': 8,
        'slk3': 7,
        'slk21': 6,
        'slk2': 5,
        'innatedb': 4,
        'reactome': 3,
        'tcr': 2,
        'acsn': 1
    }

    # Creating .tsv files for supporting and non-supporting edges
    with open('non_supp_pair.tsv', 'w') as tsvnon, open('supp_pair.tsv', 'w') as tsvsupp:
        writernon = csv.writer(tsvnon, delimiter='\t')
        writersupp = csv.writer(tsvsupp, delimiter='\t')
        writernon.writerow(['A_node_name',
                            'B_node_name',
                            'Source_db_1',
                            'Int_type_1',
                            'Source_db_2',
                            'Int_type_2'])

        writersupp.writerow(['A_node_name',
                             'B_node_name',
                             'Source_db_1',
                             'Int_type_1',
                             'Source_db_2',
                             'Int_type_2'])

        # For each database
        for db in db_list:
            db1 = db
            for other_db in db_list:
                db2 = other_db
                if db1 != db2:
                    db_dict1 = {}
                    b_db_dict1 ={}
                    with conn:
                        c = conn.cursor()
                        c.execute("SELECT * FROM EDGE WHERE source_db = 'source database:%s'" % db1)
                        # Until last row
                        while True:
                            edge_row = c.fetchone()
                            if edge_row is None:
                                break
                            else:
                                if not edge_row[3] in db_dict1:
                                    db_dict1[edge_row[3]] = [edge_row[4]]
                                else:
                                    db_dict1[edge_row[3]].append(edge_row[4])

                                if not edge_row[4] in b_db_dict1:
                                    b_db_dict1[edge_row[4]] = [edge_row[3]]
                                else:
                                    b_db_dict1[edge_row[4]].append(edge_row[3])

                    db_dict2 = {}
                    b_db_dict2 = {}
                    with conn:
                        c = conn.cursor()
                        c.execute("SELECT * FROM EDGE WHERE source_db = 'source database:%s'" % db2)
                        # Until last row
                        while True:
                            edge_row = c.fetchone()
                            if edge_row is None:
                                break
                            else:
                                if not edge_row[3] in db_dict2:
                                    db_dict2[edge_row[3]] = [edge_row[4]]
                                else:
                                    db_dict2[edge_row[3]].append(edge_row[4])

                                if not edge_row[4] in b_db_dict2:
                                    b_db_dict2[edge_row[4]] = [edge_row[3]]
                                else:
                                    b_db_dict2[edge_row[4]].append(edge_row[3])

                    # Getting edges with multiple source_dbs
                    with conn:
                        c = conn.cursor()
                        c.execute("SELECT interactor_a_node_name, interactor_b_node_name, interaction_types, source_db  FROM EDGE WHERE source_db = 'source database:%s|source database:%s'" % (db1, db2))
                        # Until last row
                        while True:
                            edge_row = c.fetchone()
                            if edge_row is None:
                                break
                            else:
                                writersupp.writerow(
                                    [edge_row[0],
                                     edge_row[1],
                                     edge_row[2],
                                     edge_row[3].split('|')[0],
                                     edge_row[3].split('|')[1]])

                    # Checking if node As from different databases interact with the same B nodes
                    anon = []
                    asupp = []

                    for key in list(db_dict1.keys()):
                        if key in list(db_dict2.keys()):
                            for val in db_dict1[key]:
                                if val in db_dict2[key]:
                                    c.execute("SELECT interactor_a_node_name, interactor_b_node_name, interaction_types, source_db FROM EDGE WHERE interactor_a_node_name = '%s' AND interactor_b_node_name = '%s'" % (
                                        key, val))
                                    edge_row = c.fetchall()

                                    # If there is only one type of interaction
                                    if edge_row[0][2].split(':')[-1] == edge_row[0][2].split(':')[2]:
                                        if edge_row[0][2].split(':')[2].split('|')[0] != edge_row[1][2].split(':')[2].split('|')[0]:
                                            anon.append([edge_row[0][0],
                                                         edge_row[0][1],
                                                         edge_row[0][3],
                                                         edge_row[0][2],
                                                         edge_row[1][3],
                                                         edge_row[1][2]])
                                            writernon.writerow(
                                                [edge_row[0][0],
                                                 edge_row[0][1],
                                                 edge_row[0][3],
                                                 edge_row[0][2],
                                                 edge_row[1][3],
                                                 edge_row[1][2]])
                                            '''
                                            # Deleting lower rank database edge data from merger database
                                            if rank_dict[edge_row[0][3].split(':')[1]] > rank_dict[edge_row[1][3].split(':')[1]]:
                                                c.execute(
                                                    "DELETE FROM EDGE WHERE interactor_a_node_name = 'uniprot:%s' AND interactor_b_node_name = 'uniprot:%s' AND source_db = 'source database:%s'" % (
                                                        edge_row[1][0].split(':')[1],
                                                        edge_row[1][1].split(':')[1],
                                                        edge_row[1][3].split(':')[1]))
                                            else:
                                                c.execute(
                                                    "DELETE FROM EDGE WHERE interactor_a_node_name = 'uniprot:%s' AND interactor_b_node_name = 'uniprot:%s' AND source_db = 'source database:%s'" % (
                                                        edge_row[0][0].split(':')[1],
                                                        edge_row[0][1].split(':')[1],
                                                        edge_row[0][3].split(':')[1]))
                                            '''
                                        else:
                                            asupp.append(
                                                [edge_row[0][0],
                                                 edge_row[0][1],
                                                 edge_row[0][3],
                                                 edge_row[0][2],
                                                 edge_row[1][3],
                                                 edge_row[1][2]])

                                            writersupp.writerow(
                                                [edge_row[0][0],
                                                 edge_row[0][1],
                                                 edge_row[0][3],
                                                 edge_row[0][2],
                                                 edge_row[1][3],
                                                 edge_row[1][2]])

                                    # Otherwise
                                    else:
                                        # If interaction types do not match
                                        # (comparison between effect-effect and molecular background-molecular background, is_direct-is_direct)
                                        if edge_row[0][2].split(':')[2].split('|')[0] != edge_row[1][2].split(':')[2].split('|')[0] or \
                                            edge_row[0][2].split(':')[-1] != edge_row[1][2].split(':')[-1] or \
                                            edge_row[0][2].split(':')[5].split('|')[0] != edge_row[1][2].split(':')[5].split('|')[0]:

                                            anon.append(
                                                [edge_row[0][0],
                                                 edge_row[0][1],
                                                 edge_row[0][3],
                                                 edge_row[0][2],
                                                 edge_row[1][3],
                                                 edge_row[1][2]])

                                            writernon.writerow(
                                                [edge_row[0][0],
                                                 edge_row[0][1],
                                                 edge_row[0][3],
                                                 edge_row[0][2],
                                                 edge_row[1][3],
                                                 edge_row[1][2]])
                                            '''
                                            # Deleting lower rank database edge data from merger database
                                            if rank_dict[edge_row[0][3].split(':')[1]] > rank_dict[
                                                edge_row[1][3].split(':')[1]]:
                                                c.execute(
                                                    "DELETE FROM EDGE WHERE interactor_a_node_name = 'uniprot:%s' AND interactor_b_node_name = 'uniprot:%s' AND source_db = 'source database:%s'" % (
                                                        edge_row[1][0].split(':')[1],
                                                        edge_row[1][1].split(':')[1],
                                                        edge_row[1][3].split(':')[1]))
                                            else:
                                                c.execute(
                                                    "DELETE FROM EDGE WHERE interactor_a_node_name = 'uniprot:%s' AND interactor_b_node_name = 'uniprot:%s' AND source_db = 'source database:%s'" % (
                                                        edge_row[0][0].split(':')[1],
                                                        edge_row[0][1].split(':')[1],
                                                        edge_row[0][3].split(':')[1]))
                                            '''
                                        else:
                                            asupp.append(
                                                [edge_row[0][0],
                                                 edge_row[0][1],
                                                 edge_row[0][3],
                                                 edge_row[0][2],
                                                 edge_row[1][3],
                                                 edge_row[1][2]])

                                            writersupp.writerow(
                                                [edge_row[0][0],
                                                 edge_row[0][1],
                                                 edge_row[0][3],
                                                 edge_row[0][2],
                                                 edge_row[1][3],
                                                 edge_row[1][2]])
                            '''
                            # Same for the other database
                            for val2 in db_dict2[key]:
                                if val2 in db_dict1[key]:
                                    c.execute(
                                        "SELECT interactor_a_node_name, interactor_b_node_name, interaction_types, source_db FROM EDGE WHERE interactor_a_node_name = '%s' AND interactor_b_node_name = '%s'" % (
                                        key, val2))
                                    edge_row2 = c.fetchall()
                                    if edge_row2 != edge_row:
                                        print(edge_row2)
                                        print(edge_row)
                                        if edge_row2[0][2].split(':')[-1] == edge_row2[0][2].split(':')[2]:
                                            if edge_row2[0][2].split(':')[2].split('|')[0] != \
                                                    edge_row2[1][2].split(':')[2].split('|')[0] or \
                                                            edge_row2[0][2].split(':')[-1] != \
                                                            edge_row2[1][2].split(':')[-1]:
                                                if [edge_row2[0][0],
                                                    edge_row2[0][1],
                                                    edge_row2[0][3],
                                                    edge_row2[0][2],
                                                    edge_row2[1][3],
                                                    edge_row2[1][2]] not in anon and bnon:
                                                    writernon.writerow(
                                                        [edge_row2[0][0],
                                                         edge_row2[0][1],
                                                         edge_row2[0][3],
                                                         edge_row2[0][2],
                                                         edge_row2[1][3],
                                                         edge_row2[1][2],
                                                         2])
                                                    
                                                    # Deleting lower rank database edge data from merger database
                                                    if rank_dict[edge_row2[0][3].split(':')[1]] > rank_dict[
                                                        edge_row2[1][3].split(':')[1]]:
                                                        c.execute(
                                                            "DELETE FROM EDGE WHERE interactor_a_node_name = 'uniprot:%s' AND interactor_b_node_name = 'uniprot:%s' AND source_db = 'source database:%s'" % (
                                                                edge_row2[1][0].split(':')[1],
                                                                edge_row2[1][1].split(':')[1],
                                                                edge_row2[1][3].split(':')[1]))
                                                    else:
                                                        c.execute(
                                                            "DELETE FROM EDGE WHERE interactor_a_node_name = 'uniprot:%s' AND interactor_b_node_name = 'uniprot:%s' AND source_db = 'source database:%s'" % (
                                                                edge_row2[0][0].split(':')[1],
                                                                edge_row2[0][1].split(':')[1],
                                                                edge_row2[0][3].split(':')[1]))
                                                    
                                            else:
                                                if [edge_row2[0][0],
                                                    edge_row2[0][1],
                                                    edge_row2[0][3],
                                                    edge_row2[0][2],
                                                    edge_row2[1][3],
                                                    edge_row2[1][2]] not in asupp:
                                                    writersupp.writerow(
                                                        [edge_row2[0][0],
                                                         edge_row2[0][1],
                                                         edge_row2[0][3],
                                                         edge_row2[0][2],
                                                         edge_row2[1][3],
                                                         edge_row2[1][2],
                                                         2])
                                        else:
                                            print(edge_row2)
                                            if edge_row2[0][2].split(':')[2].split('|')[0] != edge_row2[1][2].split(':')[2].split('|')[0] or \
                                                edge_row2[0][2].split(':')[-1] != edge_row2[1][2].split(':')[-1] or \
                                                edge_row2[0][2].split(':')[5].split('|')[0] != edge_row2[1][2].split(':')[5].split('|')[0]:
                                                if [edge_row2[0][0],
                                                    edge_row2[0][1],
                                                    edge_row2[0][3],
                                                    edge_row2[0][2],
                                                    edge_row2[1][3],
                                                    edge_row2[1][2]] not in anon and bnon:
                                                    writernon.writerow(
                                                        [edge_row2[0][0],
                                                         edge_row2[0][1],
                                                         edge_row2[0][3],
                                                         edge_row2[0][2],
                                                         edge_row2[1][3],
                                                         edge_row2[1][2],
                                                         2])
                                                    
                                                    # Deleting lower rank database edge data from merger database
                                                    if rank_dict[edge_row2[0][3].split(':')[1]] > rank_dict[
                                                        edge_row2[1][3].split(':')[1]]:
                                                        c.execute(
                                                            "DELETE FROM EDGE WHERE interactor_a_node_name = 'uniprot:%s' AND interactor_b_node_name = 'uniprot:%s' AND source_db = 'source database:%s'" % (
                                                                edge_row2[1][0].split(':')[1],
                                                                edge_row2[1][1].split(':')[1],
                                                                edge_row2[1][3].split(':')[1]))
                                                    else:
                                                        c.execute(
                                                            "DELETE FROM EDGE WHERE interactor_a_node_name = 'uniprot:%s' AND interactor_b_node_name = 'uniprot:%s' AND source_db = 'source database:%s'" % (
                                                                edge_row2[0][0].split(':')[1],
                                                                edge_row2[0][1].split(':')[1],
                                                                edge_row2[0][3].split(':')[1]))
                                                    
                                            else:
                                                if [edge_row2[0][0],
                                                    edge_row2[0][1],
                                                    edge_row2[0][3],
                                                    edge_row2[0][2],
                                                    edge_row2[1][3],
                                                    edge_row2[1][2]] not in asupp:
                                                    writersupp.writerow(
                                                        [edge_row2[0][0],
                                                         edge_row2[0][1],
                                                         edge_row2[0][3],
                                                         edge_row2[0][2],
                                                         edge_row2[1][3],
                                                         edge_row2[1][2],
                                                         2])
                            '''
                    # Checking if node Bs from different databases interact with the same A nodes
                    bnon = []
                    bsupp = []

                    for key in list(b_db_dict1.keys()):
                        if key in list(b_db_dict2.keys()):
                            for val in b_db_dict1[key]:
                                if val in b_db_dict2[key]:
                                    c.execute(
                                        "SELECT interactor_a_node_name, interactor_b_node_name, interaction_types, source_db FROM EDGE WHERE interactor_a_node_name = '%s' AND interactor_b_node_name = '%s'" % (
                                             val, key))
                                    edge_row = c.fetchall()
                                    # print(edge_row)
                                    if edge_row[0][2].split(':')[-1] == edge_row[0][2].split(':')[2]:
                                        if edge_row[0][2].split(':')[2].split('|')[0] != edge_row[1][2].split(':')[2].split('|')[0]:
                                            bnon.append(
                                                [edge_row[0][0],
                                                 edge_row[0][1],
                                                 edge_row[0][3],
                                                 edge_row[0][2],
                                                 edge_row[1][3],
                                                 edge_row[1][2]])

                                            writernon.writerow(
                                                [edge_row[0][0],
                                                 edge_row[0][1],
                                                 edge_row[0][3],
                                                 edge_row[0][2],
                                                 edge_row[1][3],
                                                 edge_row[1][2],
                                                 'b'])
                                            '''
                                            # Deleting lower rank database edge data from merger database
                                            if rank_dict[edge_row[0][3].split(':')[1]] > rank_dict[
                                                edge_row[1][3].split(':')[1]]:
                                                c.execute(
                                                    "DELETE FROM EDGE WHERE interactor_a_node_name = 'uniprot:%s' AND interactor_b_node_name = 'uniprot:%s' AND source_db = 'source database:%s'" % (
                                                        edge_row[1][0].split(':')[1],
                                                        edge_row[1][1].split(':')[1],
                                                        edge_row[1][3].split(':')[1]))
                                            else:
                                                c.execute(
                                                    "DELETE FROM EDGE WHERE interactor_a_node_name = 'uniprot:%s' AND interactor_b_node_name = 'uniprot:%s' AND source_db = 'source database:%s'" % (
                                                        edge_row[0][0].split(':')[1],
                                                        edge_row[0][1].split(':')[1],
                                                        edge_row[0][3].split(':')[1]))
                                            '''
                                        else:
                                            bsupp.append(
                                                [edge_row[0][0],
                                                 edge_row[0][1],
                                                 edge_row[0][3],
                                                 edge_row[0][2],
                                                 edge_row[1][3],
                                                 edge_row[1][2]])

                                            writersupp.writerow(
                                                [edge_row[0][0],
                                                 edge_row[0][1],
                                                 edge_row[0][3],
                                                 edge_row[0][2],
                                                 edge_row[1][3],
                                                 edge_row[1][2],
                                                 'b'])
                                    else:
                                        if edge_row[0][2].split(':')[2].split('|')[0] != edge_row[1][2].split(':')[2].split('|')[0] or \
                                            edge_row[0][2].split(':')[-1] != edge_row[1][2].split(':')[-1] or \
                                            edge_row[0][2].split(':')[5].split('|')[0] != edge_row[1][2].split(':')[5].split('|')[0]:
                                            bnon.append(
                                                [edge_row[0][0],
                                                 edge_row[0][1],
                                                 edge_row[0][3],
                                                 edge_row[0][2],
                                                 edge_row[1][3],
                                                 edge_row[1][2]])

                                            writernon.writerow(
                                                [edge_row[0][0],
                                                 edge_row[0][1],
                                                 edge_row[0][3],
                                                 edge_row[0][2],
                                                 edge_row[1][3],
                                                 edge_row[1][2],
                                                 'b'])
                                            '''
                                            # Deleting lower rank database edge data from merger database
                                            if rank_dict[edge_row[0][3].split(':')[1]] > rank_dict[
                                                edge_row[1][3].split(':')[1]]:
                                                c.execute(
                                                    "DELETE FROM EDGE WHERE interactor_a_node_name = 'uniprot:%s' AND interactor_b_node_name = 'uniprot:%s' AND source_db = 'source database:%s'" % (
                                                        edge_row[1][0].split(':')[1],
                                                        edge_row[1][1].split(':')[1],
                                                        edge_row[1][3].split(':')[1]))
                                            else:
                                                c.execute(
                                                    "DELETE FROM EDGE WHERE interactor_a_node_name = 'uniprot:%s' AND interactor_b_node_name = 'uniprot:%s' AND source_db = 'source database:%s'" % (
                                                        edge_row[0][0].split(':')[1],
                                                        edge_row[0][1].split(':')[1],
                                                        edge_row[0][3].split(':')[1]))
                                            '''
                                        else:
                                            bsupp.append(
                                                [edge_row[0][0],
                                                 edge_row[0][1],
                                                 edge_row[0][3],
                                                 edge_row[0][2],
                                                 edge_row[1][3],
                                                 edge_row[1][2]])

                                            writersupp.writerow(
                                                [edge_row[0][0],
                                                 edge_row[0][1],
                                                 edge_row[0][3],
                                                 edge_row[0][2],
                                                 edge_row[1][3],
                                                 edge_row[1][2],
                                                 'b'])
                            '''
                            # Same for the other database
                            for val2 in b_db_dict2[key]:
                                if val2 in b_db_dict1[key]:
                                    c.execute(
                                        "SELECT interactor_a_node_name, interactor_b_node_name, interaction_types, source_db FROM EDGE WHERE interactor_a_node_name = '%s' AND interactor_b_node_name = '%s'" % (
                                             val2, key))
                                    edge_row2 = c.fetchall()
                                    if edge_row2 != edge_row:
                                        # print(edge_row2)
                                        if edge_row2[0][2].split(':')[-1] == edge_row2[0][2].split(':')[2]:
                                            if edge_row2[0][2].split(':')[2].split('|')[0] != \
                                                    edge_row2[1][2].split(':')[2].split('|')[0] or \
                                                            edge_row2[0][2].split(':')[-1] != \
                                                            edge_row2[1][2].split(':')[-1]:
                                                if [edge_row2[0][0],
                                                    edge_row2[0][1],
                                                    edge_row2[0][3],
                                                    edge_row2[0][2],
                                                    edge_row2[1][3],
                                                    edge_row2[1][2]] not in anon and bnon:
                                                    writernon.writerow(
                                                        [edge_row2[0][0],
                                                         edge_row2[0][1],
                                                         edge_row2[0][3],
                                                         edge_row2[0][2],
                                                         edge_row2[1][3],
                                                         edge_row2[1][2],
                                                         'b',
                                                         2])
                                                                                                        
                                                    # Deleting lower rank database edge data from merger database
                                                    if rank_dict[edge_row2[0][3].split(':')[1]] > rank_dict[
                                                        edge_row2[1][3].split(':')[1]]:
                                                        c.execute(
                                                            "DELETE FROM EDGE WHERE interactor_a_node_name = 'uniprot:%s' AND interactor_b_node_name = 'uniprot:%s' AND source_db = 'source database:%s'" % (
                                                                edge_row2[1][0].split(':')[1],
                                                                edge_row2[1][1].split(':')[1],
                                                                edge_row2[1][3].split(':')[1]))
                                                    else:
                                                        c.execute(
                                                            "DELETE FROM EDGE WHERE interactor_a_node_name = 'uniprot:%s' AND interactor_b_node_name = 'uniprot:%s' AND source_db = 'source database:%s'" % (
                                                                edge_row2[0][0].split(':')[1],
                                                                edge_row2[0][1].split(':')[1],
                                                                edge_row2[0][3].split(':')[1]))
                                                    
                                            else:
                                                if [edge_row2[0][0],
                                                    edge_row2[0][1],
                                                    edge_row2[0][3],
                                                    edge_row2[0][2],
                                                    edge_row2[1][3],
                                                    edge_row2[1][2]] not in asupp and bsupp:
                                                    writersupp.writerow(
                                                        [edge_row2[0][0],
                                                         edge_row2[0][1],
                                                         edge_row2[0][3],
                                                         edge_row2[0][2],
                                                         edge_row2[1][3],
                                                         edge_row2[1][2],
                                                         'b',
                                                         2])
                                        else:
                                            if edge_row2[0][2].split(':')[2].split('|')[0] != edge_row2[1][2].split(':')[2].split('|')[0] or \
                                            edge_row2[0][2].split(':')[-1] != edge_row2[1][2].split(':')[-1] or \
                                            edge_row2[0][2].split(':')[5].split('|')[0] != edge_row2[1][2].split(':')[5].split('|')[0]:
                                                if [edge_row2[0][0],
                                                    edge_row2[0][1],
                                                    edge_row2[0][3],
                                                    edge_row2[0][2],
                                                    edge_row2[1][3],
                                                    edge_row2[1][2]] not in bnon and anon:
                                                    writernon.writerow(
                                                        [edge_row2[0][0],
                                                         edge_row2[0][1],
                                                         edge_row2[0][3],
                                                         edge_row2[0][2],
                                                         edge_row2[1][3],
                                                         edge_row2[1][2],
                                                         'b',
                                                         2])
                                                    
                                                    # Deleting lower rank database edge data from merger database
                                                    if rank_dict[edge_row2[0][3].split(':')[1]] > rank_dict[
                                                        edge_row2[1][3].split(':')[1]]:
                                                        c.execute(
                                                            "DELETE FROM EDGE WHERE interactor_a_node_name = 'uniprot:%s' AND interactor_b_node_name = 'uniprot:%s' AND source_db = 'source database:%s'" % (
                                                                edge_row2[1][0].split(':')[1],
                                                                edge_row2[1][1].split(':')[1],
                                                                edge_row2[1][3].split(':')[1]))
                                                    else:
                                                        c.execute(
                                                            "DELETE FROM EDGE WHERE interactor_a_node_name = 'uniprot:%s' AND interactor_b_node_name = 'uniprot:%s' AND source_db = 'source database:%s'" % (
                                                                edge_row2[0][0].split(':')[1],
                                                                edge_row2[0][1].split(':')[1],
                                                                edge_row2[0][3].split(':')[1]))
                                                    
                                            else:
                                                if [edge_row2[0][0],
                                                    edge_row2[0][1],
                                                    edge_row2[0][3],
                                                    edge_row2[0][2],
                                                    edge_row2[1][3],
                                                    edge_row2[1][2]] not in bsupp and asupp:
                                                    writersupp.writerow(
                                                        [edge_row2[0][0],
                                                         edge_row2[0][1],
                                                         edge_row2[0][3],
                                                         edge_row2[0][2],
                                                         edge_row2[1][3],
                                                         edge_row2[1][2],
                                                         'b',
                                                         2])
                            '''


# Running each function
pair_db_compare()
single_db_compare()