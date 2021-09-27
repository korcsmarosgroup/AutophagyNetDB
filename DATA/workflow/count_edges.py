"""
  Counts edges in all of the input files from each of the imported databases
  Counts edges in all of the output sql files of the imported databases
"""
import json, os, sqlite3
from collections import OrderedDict

DB_DICT = json.load(open('sources.json'), object_pairs_hook=OrderedDict)

# Counting the number of edges from the raw source files that should be inserted into the edge table
with open('edge_stats.txt', 'w+') as statfile:
    statfile.write('database\tedge number\n')
    for layer in DB_DICT.keys():
        for db in DB_DICT[layer]:
            print(db)
            counter = 0
            path = layer + '/databases/' + db + '/files/'

            if db == 'acsn':
                for file in os.listdir(path):
                    if file == 'acsn_v2_cur_prot_list.txt':
                        with open(path + file, encoding="ISO-8859-1") as infile:
                            infile.readline()
                            for line in infile:
                                # if not line.strip():
                                counter += 1

            elif db == 'innatedb':
                for file in os.listdir(path):
                    if '.mitab' in file:
                        with open(path + file, encoding="ISO-8859-1") as infile:
                            infile.readline()
                            for line in infile:
                                # if not line.strip():
                                line = line.split('\t')
                                if line[9] == 'taxid:9606(Human)' and line[10] == 'taxid:9606(Human)':
                                    counter += 1

            elif db == 'reactome':
                for file in os.listdir(path):
                    if file == 'reactome_v3.3_datafile.txt':
                        with open(path + file, encoding="ISO-8859-1") as infile:
                            infile.readline()
                            for line in infile:
                                # if not line.strip():
                                    if '9606' in line or '7227' in line \
                                            or '6239' in line or '7955' in line:
                                        counter += 1

            elif db == 'signor':
                for file in os.listdir(path):
                    with open(path + file, encoding="ISO-8859-1") as infile:
                        infile.readline()
                        for line in infile:
                            # if not line.strip():
                            cells = line.split(';')
                            type_a = cells[1].lower()
                            type_b = cells[5].lower()

                            taxids = cells[12].split(';')[0]

                            if type_a == 'protein' and type_b == 'protein' and taxids == '9606':

                                    counter += 1

            elif db == 'slk2' or db == 'slk3' or db == 'slk21' or db == 'tcr' \
                    or db == 'PSP' or db == 'OmniPath' or db == 'ComPPI' \
                    or db == 'TarBase' or db == 'miRSponge':
                for file in os.listdir(path):
                    with open(path + file, encoding="ISO-8859-1") as infile:
                        infile.readline()
                        for line in infile:
                            # if not line.strip():
                            counter += 1

            elif db == 'PhosphoSite':
                for file in os.listdir(path):
                    with open(path + file, encoding="ISO-8859-1") as infile:
                        infile.readline()
                        infile.readline()
                        infile.readline()
                        infile.readline()
                        for line in infile:
                            # if not line.strip():
                                line = line.split('\t')
                                if line[3] == 'human' and line[8] == 'human':
                                    counter += 1
            elif db == 'PTMCode2':
                for file in os.listdir(path):
                    with open(path + file, encoding="ISO-8859-1") as infile:
                        infile.readline()
                        infile.readline()
                        infile.readline()
                        infile.readline()
                        for line in infile:
                            # if not line.strip():
                            columns = line.split('\t')
                            if columns[2] == 'Homo sapiens' or columns[2] == 'Drosophila melanogaster' \
                                    or columns[2] == 'Danio rerio' \
                                    or columns[2] == 'Caenorhabditis elegans':
                                    counter += 1

            elif db == 'biogrid':
                for file in os.listdir(path):
                    with open(path + file, encoding="ISO-8859-1") as infile:
                        infile.readline()
                        for line in infile:
                            # if not line.strip():
                            counter += 1

            elif db == 'HPRD':
                path = layer + '/databases/' + db + '/files//HPRD_Release9_062910/'
                for file in os.listdir(path):
                    with open(path + file, encoding="ISO-8859-1") as infile:
                        for line in infile:
                            # if not line.strip():
                                counter += 1

            elif db == 'IntAct':
                for file in os.listdir(path):
                    with open(path + file, encoding="ISO-8859-1") as infile:
                        infile.readline()
                        for line in infile:
                            # if not line.strip():
                            columns = line.split('\t')
                            if columns[9] == 'taxid:9606(human)' and columns[10] == 'taxid:9606(human)':
                                    counter += 1

            elif db == 'miRDeathDB':
                for file in os.listdir(path):
                    with open(path + file, encoding="ISO-8859-1") as infile:
                        infile.readline()
                        for line in infile:
                            # if not line.strip():
                            counter += 1

            elif db == 'miR2Disease':
                for file in os.listdir(path):
                    with open(path + file, encoding="ISO-8859-1") as infile:
                        infile.readline()
                        infile.readline()
                        infile.readline()
                        for line in infile:
                            # if not line.strip():
                                counter += 1

            elif db == 'miRecords':
                for file in os.listdir(path):
                    with open(path + file, encoding="ISO-8859-1") as infile:
                        infile.readline()
                        for line in infile:
                            # if not line.strip():
                                line = line.split('\t')
                                if line[1] == 'Homo sapiens' and line[5] == 'Homo sapiens' \
                                    or line[1] == 'Drosophila melanogaster' and line[5] == 'Drosophila melanogaster' \
                                    or line[1] == 'Caenorhabditis elegans' and line[5] == 'Caenorhabditis elegans' \
                                    or line[1] == 'Danio rerio' and line[5] == 'Danio rerio':
                                    counter += 1

            elif db == 'starBase' or db == 'starbase':
                for file in os.listdir(path):
                    with open(path + file, encoding="ISO-8859-1") as infile:
                        infile.readline()
                        infile.readline()
                        infile.readline()
                        infile.readline()
                        for line in infile:
                            # if not line.strip():
                            counter += 1

            elif db == 'TFBS':
                for file in os.listdir(path):
                    with open(path + file, encoding="ISO-8859-1") as infile:
                        for line in infile:
                            # if not line.strip():
                                counter += 1

            elif db == 'TFBS_Orsi':
                for file in os.listdir(path):
                    with open(path + file, encoding="ISO-8859-1") as infile:
                        for line in infile:
                            # if not line.strip():
                                if 'Homo sapiens' in line \
                                    or 'Drosophila melanogaster' in line \
                                    or 'Caenorhabditis elegans' in line \
                                    or 'Danio rerio' in line:
                                    counter += 1

            elif db == 'lncRInter':
                for file in os.listdir(path):
                    with open(path + file, encoding="ISO-8859-1") as infile:
                        infile.readline()
                        for line in infile:
                            # if not line.strip():
                                if 'Homo sapiens' in line:
                                    counter += 1

            elif db == 'NPInter':
                for file in os.listdir(path):
                    with open(path + file, encoding="ISO-8859-1") as infile:
                        infile.readline()
                        for line in infile:
                            # if not line.strip():
                                if 'Homo sapiens' in line \
                                        or 'Drosophila melanogaster' in line \
                                        or 'Caenorhabditis elegans' in line \
                                        or 'Danio rerio' in line:
                                    counter += 1

            statfile.write(db + '\t' + str(counter) + '\n')

sql_output_path = 'all_output/'

# Counting the number of edges inserted into the database SQL table
with open('edge_stats.txt', 'r') as statfile, open('edge_stats_parsed.txt', 'w') as statfile_parsed:
    statfile_parsed.write('database\tedge number\tparsed edge number\n')
    statfile.readline()

    for line in statfile:
        line = line.strip().split('\t')
        dbname = line[0] + '.db'
        conn = sqlite3.connect(sql_output_path + dbname)
        with conn:
            c = conn.cursor()
            num = c.execute("SELECT COUNT(*) FROM edge")
            edges = num.fetchone()[0]

            statfile_parsed.write('\t'.join(line) + '\t' + str(edges) + '\n')
