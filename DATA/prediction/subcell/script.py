"""
 Parses the Integrated subcellular localization dataset files of 3 species (human, drosophila, celegans) from:
 http://comppi.linkgroup.hu/downloads term txt 2nd és 4th columns
"""

# Imports
import sqlite3
import logging

# Constants
FILE_LIST = [
    '../prediction/subcell/files/comppi--proteins_locs--tax_celegans_loc_all.txt',
    '../prediction/subcell/files/comppi--proteins_locs--tax_dmelanogaster_loc_all.txt',
    '../prediction/subcell/files/comppi--proteins_locs--tax_hsapiens_loc_all.txt'
]
GO_FILE = '../prediction/subcell/files/go_daily-termdb-tables/term.txt'
builder = 'SLK3_layers.db'

def subcell_prediction(logger):
    conn = sqlite3.connect(builder)

    with conn:
        c = conn.cursor()

        mapdict = {}

        # Parsing data file
        print(f'Get information from file: {GO_FILE}')
        with open(GO_FILE) as gofile:
            for line in gofile:
                line = line.strip().split('\t')
                godescr = line[1]
                goid_fetched = line[3]
                # Creating dictionary with id as key and description as value
                mapdict[goid_fetched] = godescr

        for file in FILE_LIST:

            print(f'Get information from file: {file}')
            with open(file) as data:
                # Skipping the header
                data.readline()

                for line in data:
                    line = line.strip().split('\t')

                    uniprot = f'Uniprot:{line[0]}'
                    majorloc = line[3].split('|')
                    minorloc = line[4].split('|')

                    # Mapping GO ids
                    minorloc_mapped = []
                    for item in minorloc:
                        if item in mapdict.keys():
                            minorloc_mapped.append(''.join([item, '(', mapdict[item], ')']))
                        else:
                            minorloc_mapped.append(''.join([item, '(', item, ')']))

                    # Concatenating all major and minor locs for a protein
                    minorloc_final = 'minorloc:' + '|minorloc:'.join(minorloc_mapped)

                    # Creating string with major and minor loc separated by |
                    majorloc_names = []
                    for elem in majorloc:
                        majorloc_name = elem.split(':')[0]
                        majorloc_names.append(majorloc_name)

                    majorloc_add = 'majorloc:' + '|majorloc:'.join(majorloc_names)
                    majorloc_add = '|' + majorloc_add

                    if "'" in minorloc_final:
                        minorloc_final = minorloc_final.replace("'", "")
                    minorloc_add = '|' + minorloc_final

                    # Adding data to merger
                    c.execute(f"UPDATE node SET topology = node.topology || '{majorloc_add + minorloc_add}'"
                              f"WHERE node.name = '{uniprot}'")


if __name__ == '__main__':
    subcell_prediction(logger=None)
