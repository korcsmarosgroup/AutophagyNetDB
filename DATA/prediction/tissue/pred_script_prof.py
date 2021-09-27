"""
 Tissue expression data from: https://bgee.org/?page=download&action=proc_values#id5
"""

# Imports
import logging, os, sqlite3, json, io

# Initiating logger
logger = logging.getLogger()
handler = logging.FileHandler('../../workflow/SLK3.log')
logger.setLevel(logging.DEBUG)
handler.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)

# Getting file names
RNA_SEQ_LIST = []
AFF_LIST = []
for directory in os.listdir('files'):
    if directory.endswith('.zip'):
        pass
    else:
        if 'RNA' in directory:
            for file in os.listdir('files/' + directory):
                RNA_SEQ_LIST.append('files/' + directory + '/' + file)
        elif 'Affymetrix' in directory:
            for file in os.listdir('files/' + directory):
                AFF_LIST.append('files/' + directory + '/' + file)

# Initiating mapper
map_db = sqlite3.connect('../../workflow/mapper.db')

#logger.debug('Initiating mapper')
#map_conn = sqlite3.connect('../../workflow/mapper.db')

#tempfile = io.StringIO()
#for line in map_conn.iterdump():
#     tempfile.write('%s\n' % line)
#map_conn.close()
#tempfile.seek(0)

# Create a database in memory and import from tempfile
#map_db = sqlite3.connect(":memory:")
#with map_db:
#  map_db.cursor().executescript(tempfile.read())

# Getting data
logger.debug('Started extracting experiment data')
exp_dict = {}
counter = 0
i=0
with map_db:
    c2 = map_db.cursor()
    for file in RNA_SEQ_LIST:
        # Limit TODO remove limit
        counter += 1
        logger.debug(counter)
        if counter == 2:
            break
        # Opening the library file
        i = 0
        with open(file) as lib_file:
            lib_file.readline()
            for line in lib_file:
                i += 1
                print(i)
                if i == 10:
                    break
                line = line.strip().split('\t')
                c2.execute(
                    "SELECT UNIPROT_AC.uniprot_ac FROM UNIPROT_AC JOIN MAPP ON MAPP.uniprot_ac=UNIPROT_AC.id "
                    "WHERE MAPP.foreign_id='%s' GROUP BY MAPP.foreign_id"
                    % line[3]
                )
                while True:
                    allrows = c2.fetchone()
                    if allrows is None:
                        break
                    else:
                        gene_name = allrows[0]
                        if gene_name not in exp_dict:
                            exp_dict[gene_name] = {'value': [],
                                                   'url': [],
                                                   'db': 'Uber Anatomy Ontology',
                                                   'searchable': True,
                                                   }
                        else:
                            # Adding data if it's not already in the dictionary
                            if line[5].replace('"', '') not in exp_dict[gene_name]['value'] \
                                and line[4] not in exp_dict[gene_name]['url'] \
                                    and line[7].replace('"','') not in exp_dict[gene_name]['url'] \
                                    and line[12] not in exp_dict[gene_name]['url'] \
                                    and line[13] not in exp_dict[gene_name]['url']:
                                        exp_dict[gene_name]['value'].append(line[5].replace('"',''))
                                        exp_dict[gene_name]['url'].append('|uberon:' + line[4])
                                        exp_dict[gene_name]['url'].append('|stage:' + line[7].replace('"', ''))
                                        exp_dict[gene_name]['url'].append('|FPKM:' + line[12])
                                        exp_dict[gene_name]['url'].append('|presence:' + line[13])
                            print(exp_dict[gene_name])
