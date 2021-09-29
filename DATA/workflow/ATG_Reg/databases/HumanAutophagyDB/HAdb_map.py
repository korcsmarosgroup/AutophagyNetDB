"""
  :param: protlist_file: data copied manually from Human autophagy database: http://autophagy.lu/clustering/
  :param: atg_file: Autophagy database, data from: http://www.tanpaku.org/autophagy/download/atg_gene.txt
  :param: hadb_out: human autophagy database mapped output file
  :param: atg_out: autophagy database mapped output file

"""

# Imports
import sqlite3, os

# Defining constants
# protlist_file = 'HAdb.txt'
# atg_file = 'atg_gene.txt'
#
# hadb_out = 'HAdb_mapped.txt'
# atg_out = 'atg_mapped.txt'
#
# mapper_db = '../mapper.db'
#
# db_destination = '../all_output/HumanAutophagyDB.db'


def hadb_parse(hadb_infile, hadb_outfile):
    conn = sqlite3.connect(mapper_db)

    # Starting connection to mapper dataset
    with conn:
        c = conn.cursor()
        # Parsing file, getting gene IDs
        with open(hadb_infile, encoding='ISO-8859-1') as infile1, open(hadb_outfile, 'w') as outfile1:
            # Header
            outfile1.write('Gene ID\tUniprot ID\n')
            # Human autophagy database
            for line in infile1:
                line = line.strip().split('\t')
                geneid = line[2]
                # Mapping gene IDs to uniprot IDs
                c.execute('''SELECT UNIPROT_AC.uniprot_ac FROM UNIPROT_AC JOIN MAPP ON MAPP.uniprot_ac=UNIPROT_AC.id
                          WHERE MAPP.foreign_id='%s' GROUP BY MAPP.foreign_id''' % geneid)
                while True:
                    allrows = c.fetchone()
                    if allrows is None:
                        break
                    else:
                        uniprot = 'uniprot:' + allrows[0]

                        # Writing data to outfile
                        outfile1.write('\t'.join(['HumanAutophagyDB:' + geneid, uniprot]) + '\n')


def atgdb_parse(atg_infile, atg_outfile):
    conn = sqlite3.connect(mapper_db)

    with conn:
        c2 = conn.cursor()
        with open(atg_infile, encoding='ISO-8859-1') as infile2, open(atg_outfile, 'w') as outfile2:
            outfile2.write('Gene ID\tUniprot ID\n')

            # Autophagy database
            for line in infile2:
                line = line.strip().split('\t')
                geneid = line[1]
                # Mapping gene IDs to uniprot IDs
                c2.execute('''SELECT UNIPROT_AC.uniprot_ac FROM UNIPROT_AC JOIN MAPP ON MAPP.uniprot_ac=UNIPROT_AC.id
                          WHERE MAPP.foreign_id='%s' GROUP BY MAPP.foreign_id''' % geneid)
                while True:
                    allrows = c2.fetchone()
                    if allrows is None:
                        break
                    else:
                        uniprot = 'uniprot:' + allrows[0]

                        # Writing data to outfile
                        outfile2.write('\t'.join(['AutophagyDB' + geneid, uniprot]) + '\n')


def main(logger):

    # Remove output file if it already exists
    try:
        os.remove(hadb_out)
        os.remove(atg_out)
    except OSError:
        pass

    # Calling parsing functions
    hadb_parse(protlist_file, hadb_out)
    atgdb_parse(atg_file, atg_out)

    # Setting up node table
    dest_conn = sqlite3.connect(db_destination)

    with dest_conn:
        c = dest_conn.cursor()
        c.execute("DROP TABLE IF EXISTS node")
        c.execute('''CREATE TABLE node (id INT,
                                           geneID CHAR(100),
                                           uniprot CHAR(100))''')
        # Getting data
        with open(hadb_out) as hadb, open(atg_out) as atg:
            hadb.readline()
            for line in hadb:
                line = line.strip().split('\t')
                gene = line[0]
                uniprot = line[1]
                # Adding node data to sql table
                c.execute("INSERT INTO node VALUES (?,?)", (gene, uniprot))
            atg.readline()
            for line in atg:
                line = line.strip().split('\t')
                gene = line[0]
                uniprot = line[1]
                # Adding node data to sql table
                c.execute("INSERT INTO node VALUES (?,?)", (gene, uniprot))


if __name__ == '__main__':
    main(None)