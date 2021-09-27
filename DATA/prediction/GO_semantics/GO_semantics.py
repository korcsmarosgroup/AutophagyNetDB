import sqlite3
import subprocess

merge_conn = sqlite3.connect('../../workflow/merger.db', isolation_level=None)
outfile = 'GO_input.tsv'
go_output_score_files = [
#    ["C", "CC", "GO_score_CC.tsv"],
    ["P", "BP", "GO_score_BP.tsv"],
#    ["F", "MF", "GO_score_MF.tsv"],
]

#with merge_conn:
#    c = merge_conn.cursor()
#    with open(outfile, 'w') as out:
#        c.execute("SELECT interactor_a_node_name, interactor_b_node_name from EDGE")
#        while True:
#            allrows = c.fetchone()
#            if allrows is None:
#                break
#            else:
#                nodea = allrows[0].split(":")[1]
#                nodeb = allrows[1].split(":")[1]
#                out.write(f'{nodea}\t{nodeb}\n')

#for go in go_output_score_files:
#    subprocess.call(["java", "-jar", "HVSM.jar", "-org", "human", "-db", go[0], "-gene", "-i", outfile, "-o", go[2]])

# Adding scores to SQL
for go in go_output_score_files:

    GOdict = []

    with open(go[2], "r") as score_file, open(outfile, 'r') as infile:
        for nodes, score in zip (infile, score_file):
            score = float(score.strip())
            node = nodes.strip().split("\t")
            if score > 0:
                GOdict.append([node[0], node[1], score])

    with merge_conn:
        c = merge_conn.cursor()
        #c.execute('PRAGMA journal_mode=OFF')
        #c.execute('PRAGMA synchronous=OFF')
        #c.execute('CREATE INDEX edge_node_a_b_index ON edge (interactor_a_node_name, interactor_b_node_name);')
        sum = len(GOdict)
        num = 0
        for GOscore in GOdict:
            c.execute("UPDATE edge SET confidence_scores = edge.confidence_scores || %s"
                        "WHERE edge.interactor_a_node_name = '%s' AND edge.interactor_b_node_name = '%s'"
                        % (f'\"|GO_semantics_{go[1]}:{GOscore[2]}|\"', f'Uniprot:{GOscore[0]}', f'Uniprot:{GOscore[1]}'))
            num += 1
            if num % 100 == 0:
                print(f'{num}/{sum}')
        #merge_conn.commit()
