
import sqlite3, json, io, logging, functools
mapper_location = '../DATA/workflow/mapper.db'

json_mapper_file = '../DATA/workflow/uniprot_id_mapping.json'


def map_uniprot_to_external(json_mapper_file):
    """
    Maps primary uniprot ids to different external source ids based on the mapping table
    Also retrievs the database which the id is from (id type)
    :param uniprot: uniprot id that we want to fing external references for
    :return: external reference id and database
    """
    mapdict = {}
    with open(json_mapper_file) as mapfile:
        for line in mapfile:
            line = json.loads(line)
            extrefdb = line['to_id_type']
            extrefid = line['to_id'].upper()
            if line['from_id_type'] == 'uniprotac':
                uniprot = line['from_id'].upper()
                if extrefdb == 'uniprotac' or \
                        extrefdb == 'hgnc' or \
                        extrefdb == 'ensembl':
                    if uniprot not in mapdict:
                        mapdict[uniprot] = {extrefdb: [extrefid]}
                    else:
                        if extrefdb not in mapdict[uniprot]:
                            mapdict[uniprot][extrefdb] = [extrefid]
                        else:
                            mapdict[uniprot][extrefdb].append(extrefid)
    return mapdict


mapdict = map_uniprot_to_external(json_mapper_file=json_mapper_file)

for key, value in mapdict['Q7Z494'].items():
    print(key)

