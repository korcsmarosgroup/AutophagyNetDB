"""
  Parses the mapped Human autophagy database list and finds any edges from third party databases
  if they contain the listed genes
  :param merger = merger database file THIS USES UNIPROT IDS
"""

import sqlite3

HA_DATA_FILE = 'HAdb_mapped.txt'
ATG_DATA_FILE = 'atg_mapped.txt'
mergerdb = '../../merger.db'

# Setting up connection to merger database
conn = sqlite3.connect(mergerdb)


def main(log):
    with conn:
        c = conn.cursor()
        with open(HA_DATA_FILE) as ha_infile:
            ha_infile.readline()
            ATG_PROT_LIST = []
            for line in ha_infile:
                line = line.strip().split('\t')
                uniprot = line[1]
                if uniprot not in ATG_PROT_LIST:
                    ATG_PROT_LIST.append(uniprot)
        with open(ATG_DATA_FILE) as atg_infile:
            atg_infile.readline()
            for line in atg_infile:
                line = line.strip().split('\t')
                uniprot = line[1]
                if uniprot not in ATG_PROT_LIST:
                    ATG_PROT_LIST.append(uniprot)

        c.execute("SELECT * FROM edge")
        while True:
            row = c.fetchone()
            # until the last row
            if row is None:
                break
            else:
                print(row)
                a_node = row[3]
                b_node = row[4]

                for item in ATG_PROT_LIST:
                    if item == a_node:
                        # Add edge to L1 to database
                    elif item == b_node:
                        # Add edge to L1 database
                    else:
                        # Do nothing




