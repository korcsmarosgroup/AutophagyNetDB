"""
  Parses disease data from:
        https://cancer.sanger.ac.uk/census
        http://www.disgenet.org/web/DisGeNET/menu
  Adding it as a node attribute
"""

# Imports
import sqlite3
import logging
import json

# Constants
SANGER_DATA = 'files/Census_allFri Aug 10 18_20_49 2018.tsv'
DISGEN_DATA = 'files/all_gene_disease_associations.tsv'

merger = '../../workflow/merger.db'  # builder !!!
mapper = '../../workflow/mapper.db'

# Initiating logger
logger = logging.getLogger()
handler = logging.FileHandler('SLK3.log')
logger.setLevel(logging.DEBUG)
handler.setLevel(logging.DEBUG)
formatter = logging.Formatter(
    '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
# Establishing database connections
merge_conn = sqlite3.connect(merger)
map_conn = sqlite3.connect(mapper)


def main(logger):
    dis_dict = {}
    sang_dict = {}
    with map_conn:
        c = map_conn.cursor()
        with open(DISGEN_DATA, encoding="ISO-8859-1") as disgen:
            disgen.readline()
            for line in disgen:
                line = line.strip().split('\t')
                gen_name = line[1]
                dis_id = line[2]
                dis_name = line[3]
                dis = u':'.join((dis_id, dis_name)).encode(
                    'utf8').strip().decode('utf8')
                # Mapping gene names to uniprot ids
                c.execute(
                    "SELECT UNIPROT_AC.uniprot_ac FROM UNIPROT_AC JOIN MAPP ON MAPP.uniprot_ac=UNIPROT_AC.id "
                    "WHERE MAPP.foreign_id='%s' GROUP BY MAPP.foreign_id"
                    % gen_name)
                while True:
                    allrows = c.fetchone()
                    if allrows is None:
                        break
                    else:
                        uni_name = allrows[0]
                        if uni_name not in dis_dict:
                            dis_dict[uni_name] = [{'value': dis,
                                                   'url': 'http://www.disgenet.org/web/DisGeNET/menu',
                                                   'db': 'DisGeNET',
                                                   'searchable': True,
                                                   }]
                        else:
                            if dis not in dis_dict[uni_name][0]['value']:
                                dis_dict[uni_name].append({'value': dis,
                                                           'url': 'http://www.disgenet.org/web/DisGeNET/menu',
                                                           'db': 'DisGeNET',
                                                           'searchable': True,
                                                           })

        # Sanger data
        with open(SANGER_DATA) as sandata:
            sandata.readline()
            for line in sandata:
                line = line.strip().split('\t')
                gen_name = line[0]
                cancer = line[9]
                role = line[11]
                # Assembling data
                can = u'|'.join((cancer, role)).encode(
                    'utf8').strip().decode('utf8')
                # TODO: add 'cancer:' in front of string

                # Mapping gene names to uniprot ids
                c.execute(
                    "SELECT UNIPROT_AC.uniprot_ac FROM UNIPROT_AC JOIN MAPP ON MAPP.uniprot_ac=UNIPROT_AC.id "
                    "WHERE MAPP.foreign_id='%s' GROUP BY MAPP.foreign_id"
                    % gen_name)
                while True:
                    allrows = c.fetchone()
                    if allrows is None:
                        break
                    else:
                        uni_name = allrows[0]
                        if uni_name not in sang_dict:
                            sang_dict[uni_name] = [{'value': can,
                                                    'url': 'https://cancer.sanger.ac.uk/census',
                                                    'db': 'DisGeNET',
                                                    'searchable': True,
                                                    }]
                        else:
                            if can not in sang_dict[uni_name][0]['value']:
                                sang_dict[uni_name].append({'value': can,
                                                            'url': 'https://cancer.sanger.ac.uk/census',
                                                            'db': 'DisGeNET',
                                                            'searchable': True,
                                                            })

    # Adding data to merger
    logger.debug('Adding data to the nodes of merger')
    with merge_conn:
        c2 = merge_conn.cursor()
        for key, value in dis_dict.items():
            c2.execute("UPDATE node SET topology = node.topology || '%s' WHERE node.name = '%s'"
                       % ('|disease:' + json.dumps(value), 'uniprot:' + key))
        c3 = merge_conn.cursor()
        for key, value in sang_dict.items():
            c3.execute("UPDATE node SET topology = node.topology || '%s' WHERE node.name = '%s'"
                       % ('|disease:' + json.dumps(value), 'uniprot:' + key))
    logger.debug('Disease data done')


if __name__ == '__main__':
    main(logger=logger)
