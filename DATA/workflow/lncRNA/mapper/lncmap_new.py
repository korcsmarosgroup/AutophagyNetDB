"""
 Translates lncRNA and miRNA ids to RNACentral ids
"""

# Imports
import sqlite3
import itertools as it
from collections import defaultdict

# Defining constants
SQL_SEED = '../../../../SLKlib/SQLiteDBApi/network-db-seed.sql'
noncode_file = 'noncode.tsv'
ens_file = 'ensembl.tsv'
mir_file = 'mirbase.dat'
mir_file_2 = 'mirbase.tsv'
lncipedia_file = 'lncipedia.tsv'
hgnc_file = 'hgnc.tsv'
DB_DESTINATION = '../../workflow/lncrna_map.db'


def main(logger):
    # Building the output SQL table
    conn = sqlite3.connect(DB_DESTINATION)
    with conn:
        c = conn.cursor()
        c.execute('DROP table IF EXISTS mapper')
        c.execute('CREATE table mapper(id INTEGER PRIMARY KEY,'
                  'source_type TEXT,'
                  'orig_ac TEXT,'
                  'mapped_ac TEXT)')

        # Opening map files
        # NONCODE
        with open(noncode_file) as noncode:
            for line in noncode:
                line = line.strip().split('\t')
                # If our organism
                if line[3] in ['9606', '7227', '6239', '7955']:
                    c.execute('INSERT INTO mapper(source_type, orig_ac, mapped_ac) VALUES(?,?,?)',
                              (line[1], line[5].split('.')[0].lower(), line[0]))
        # Ensembl
        with open(ens_file) as ensembl:
            for line in ensembl:
                line = line.strip().split('\t')
                # If our organism
                if line[3] in ['9606', '7227', '6239', '7955']:
                    c.execute('INSERT INTO mapper(source_type, orig_ac, mapped_ac) VALUES(?,?,?)',
                              (line[1], line[2].lower(), line[0]))

        # Mirbase
        nonvalid= []
        valid = 0

        def parse_mirbase_file(path):
            """ Extract the id and ac from mirbase.dat file """

            with open(path, 'r') as file:
                mirna = []

                for line in file:
                    if line.startswith('ID'):
                        setup_id = []
                        mirna.append(setup_id)

                    mirna[-1].append(line)
                result = []
                for group in mirna:
                    name = group[0][5:22].strip()
                    identifier = group[2][3:-2].strip()

                    extracted_data = {'id': name,
                                      'ac': identifier}

                    result.append(extracted_data)

            return result

        nametoac = parse_mirbase_file(mir_file)

        # Mapping mirids to urs ids (URS00000011DF)
        with open(mir_file_2) as mirbase2:
            for line in mirbase2:
                line = line.strip().split('\t')
                # If our organism
                if line[3] in ['9606', '7227', '6239', '7955']:
                    for item in nametoac:
                        if line[2] in item.values():
                            c.execute('INSERT INTO mapper(source_type, orig_ac, mapped_ac) VALUES(?,?,?)',
                                          (line[1], item['id'].lower(), line[0]))
                    # Logging ids which names were not mapped
                    else:
                        nonvalid.append(line[2])

        print(nonvalid)

        # Lncipedia
        with open(lncipedia_file) as lncipedia:
            for line in lncipedia:
                line = line.strip().split('\t')
                # If our organism
                if line[3] in ['9606', '7227', '6239', '7955']:
                    c.execute('INSERT INTO mapper(source_type, orig_ac, mapped_ac) VALUES(?,?,?)',
                              (line[1], line[2].lower(), line[0]))

        # HGNC
        with open(hgnc_file) as hgnc:
            for line in hgnc:
                line = line.strip().split('\t')
                # If our organism
                if line[3] in ['9606', '7227', '6239', '7955']:
                    c.execute('INSERT INTO mapper(source_type, orig_ac, mapped_ac) VALUES(?,?,?)',
                              (line[1], line[5].lower(), line[0]))


if __name__=='__main__':
    main(logger=None)
