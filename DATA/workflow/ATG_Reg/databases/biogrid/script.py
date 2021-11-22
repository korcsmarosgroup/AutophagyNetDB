# -*- coding: utf-8 -*-

"""
This script parses the automatically downloaded Biogrid (psimi-tab) files to a SQLite db
    :argument: RAW_FILE: psimi-tab files
    :argument: DB_TYPE: name of the database
    :argument: DB_DESTINATION: saving location
"""

# Imports
import sys, os
from ARNlib.SQLiteDBApi.sqlite_db_api import PsimiSQL

# Defining constants
SQL_SEED = '../../../../../../../ARNlib/SQLiteDBApi/network-db-seed.sql'
RAW_FILE = '../../../ATG_Reg/databases/biogrid/files/BIOGRID-SYSTEM-3.4.145.mitab/BIOGRID-SYSTEM-Affinity_Capture-Luminescence-3.4.145.mitab.txt'  # BioGRID (psimi-tab)
DB_TYPE = 'biogrid'
DB_DESTINATION = '../../output/'


def main(logger):
    file_ = open(RAW_FILE)
    file_.seek(0, os.SEEK_END)

    filesize = file_.tell()
    filesize_mb = filesize / (1024 * 1024)

    # reseting the iterator to the begining of the files
    file_.seek(0)
    file_.readline()

    # Creating a psimi to sql db to every 15Mb of the raw Biogrid files

    # Setting the size of the pice
    mb = 1024*1024
    piece_size = 10 * mb

    # The number of the little files
    file_counter = 0

    while file_.tell() < filesize:
        starting_position = file_.tell()
        parser = PsimiSQL(SQL_SEED)

        while file_.tell() < starting_position + piece_size:
            sys.stdout.write("Parsing piece: %d Mb / %d Mb  Total: %d Mb / %d Mb \r" % ((file_.tell()-starting_position)/(1024*1024) , piece_size/(1024 * 1024),file_.tell()/(1024*1024),filesize_mb))

            # Dealing with the data
            line = file_.readline()
            cells = line.split("\t")

            try:
                # Extracting node a's properties
                node_a_dict = {
                    'name' : extract_property("biogrid", cells[2]),
                    'alt_accession' : extract_property("locuslink", cells[2]),
                    'tax_id' : cells[9],
                    'pathways' : None,
                    'aliases' : None
                }

                # Extracting node b's properties
                node_b_dict = {
                    'name' : extract_property("biogrid", cells[3]),
                    'alt_accession' : extract_property("locuslink", cells[3]),
                    'tax_id' : cells[10],
                    'pathways' : None,
                    'aliases' : None
                }

                # Interaction types
                inttype = cells[11].replace('psi-mi:', '').replace('"', '')

                if inttype == 'MI:0407(direct interaction)':
                    is_direct = inttype
                    effect = 'MI:0190(interaction type)'
                else:
                    is_direct = 'unknown'
                    effect = inttype

                interaction_types = "effect:%s|is_directed:%s|is_direct:%s" \
                                    % (effect, "directed", is_direct)

                # Extracting the edge's properties
                edge_dict = {
                    'interaction_detection_method' : cells[6].replace('psi-mi:', '').replace('"', ''),
                    'first_author' : cells[7],
                    'publication_ids' : cells[8],
                    'interaction_types' : interaction_types,
                    'source_db' : 'biogrid',
                    'interaction_identifiers': None,
                    'confidence_scores' : cells[14],
                    'layer' : "1"
                }

                # Inserting interactor a to the node table
                parser.insert_node(node_a_dict)

                # Inserting interactor b to the node table
                parser.insert_node(node_b_dict)

                # After insertion the node dictionaries will contain a lastrowid property

                # Inserting edge
                parser.insert_edge(node_a_dict,node_b_dict,edge_dict)

                # Inserting aliases

                #aliases_a = cells[4]
                #aliases_b = cells[5]

                #parser.insert_aliases(node_a_dict,aliases_a)
                #parser.insert_aliases(node_b_dict,aliases_b)
            except IndexError:
                break

        parser.save_db_to_file(DB_DESTINATION + "db_piece_%d" % file_counter)
        parser.db.close()
        sum_files = filesize/piece_size
        #sys.stdout.write('%d / %d SQLite db pieces created\r' % (file_counter, sum_files))
        file_counter += 1

    print("Data insertion is completed")


def extract_property(search_str, str):
    """
    This function extracts the first string that matches a given search string from Psi-Mi formatted table's cell
    :param search_str: A Psi-Mi property that needs to be extracted from a Psi-Mi cell
    :type search_str: string
    :param str: A cell of the Psi-Mi formatted table
    :type str: string
    :return:
    """
    if "|" in str:
        str_arr = str.split("|")
        for word in str_arr:
            if search_str in word:
                return word
    elif search_str in str:
        return str

    return ""

if __name__ == '__main__':
    main(logger = None)
