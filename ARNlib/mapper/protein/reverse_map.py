"""
 Maps from Uniprot ids to external names
"""

# Imports
import sqlite3, io

merger_conn = sqlite3.connect('../../../DATA/workflow/merger.db')

# Connecting to mapper database
# Read database to tempfile
map_conn = sqlite3.connect('../../../DATA/workflow/test.db')
tempfile = io.StringIO()
for line in map_conn.iterdump():
    tempfile.write('%s\n' % line)
map_conn.close()
tempfile.seek(0)

# Create a database in memory and import from tempfile
map_db = sqlite3.connect(":memory:")
with map_db:
    map_db.cursor().executescript(tempfile.read())

# Fetching merger node data
with merger_conn:
    c = merger_conn.cursor()
    c.execute("SELECT * FROM node")
    while True:
        row = c.fetchone()
        if row is None:
            break
        else:
            # Assigning ids that we want to be mapped
            foreign_id = row[1].split(':')[1]
            # Connecting to mapper database
            with map_db:
                c2 = map_db.cursor()
                # Getting gene names and external protein names where uniprot id matches
                c2.execute(
                    "SELECT MAPP.gene_name, MAPP.prot_full_name FROM MAPP JOIN UNIPROT_AC ON UNIPROT_AC.id=MAPP.uniprot_ac "
                    "WHERE UNIPROT_AC.uniprot_ac='%s' group by UNIPROT_AC.uniprot_ac"
                    % foreign_id
                )
                while True:
                    allrows = c2.fetchone()
                    if allrows is None:
                        break
                    else:
                        gene_name = allrows[0]
                        prot_ext_name = allrows[1]
                        #print(gene_name, foreign_id)
