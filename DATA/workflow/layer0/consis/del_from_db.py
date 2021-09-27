import sqlite3, csv

merger = '../../merger.db'
conn = sqlite3.connect(merger)

# Assigning ranks to each database
rank_dict = {
    'signor' : 8,
    'slk3' : 7,
    'slk21' : 6,
    'slk2' : 5,
    'innatedb' : 4,
    'reactome' : 3,
    'tcr' : 2,
    'acsn' : 1
}

nonpair = csv.reader(open('non_supp_pair.tsv'), delimiter ='\t')

next(nonpair) # Skipping the header

for line in nonpair:
    with conn:
        c = conn.cursor()

        # Deleting interaction data from the weaker source
        if rank_dict[line[2].split(':')[1]] > rank_dict[line[4].split(':')[1]]:
            c.execute(
                "DELETE FROM EDGE WHERE interactor_a_node_name = 'uniprot:%s' AND interactor_b_node_name = 'uniprot:%s' AND source_db = 'source database:%s'" % (
                    line[0].split(':')[1],
                    line[1].split(':')[1],
                    line[4].split(':')[1]))
        else:
            c.execute(
                "DELETE FROM EDGE WHERE interactor_a_node_name = 'uniprot:%s' AND interactor_b_node_name = 'uniprot:%s' AND source_db = 'source database:%s'" % (
                    line[0].split(':')[1],
                    line[1].split(':')[1],
                    line[2].split(':')[1]))

nonsingle = csv.reader(open('non_supporting.tsv'), delimiter ='\t')

next(nonsingle) # Skipping the header

for row in nonsingle:
    with conn:
        c = conn.cursor()

        # Deleting interaction data originating from both sources
        c.execute(
            "DELETE FROM EDGE WHERE interactor_a_node_name = 'uniprot:%s' AND interactor_b_node_name = 'uniprot:%s' AND source_db = 'source database:%s'" % (
                row[0].split(':')[1],
                row[1].split(':')[1],
                row[3]))