# -*- coding: utf-8 -*-

"""
    Creating database for innatedb data
    :argument: DATA_FILE: PPI files in mitab format: http://www.innatedb.com/download/interactions/innatedb_ppi.mitab.gz
    :argument: DB_TYPE: name of the source database
    :argument: EXPORT_DB_LOCATION saving location of the created database files
    :argument: XGMML_LIST: list of .xgmml files for each pathway
"""

from ARNlib.SQLiteDBApi.sqlite_db_api import PsimiSQL
import xml.etree.ElementTree as ET
import csv

#Defining constants
SQL_SEED = '../../../../../ARNlib/SQLiteDBApi/network-db-seed.sql'
DATA_FILE = 'files/innatedb_ppi.mitab'
DB_TYPE = 'InnateDB'
EXPORT_DB_LOCATION = '../../../layer0/output/innatedb'
XGMML_LIST = ['files/cytosolic_dna_sensing.xgmml',
              'files/nod_like.xgmml',
              'files/rig1.xgmml',
              'files/jakstat.xgmml',
              'files/mapk.xgmml',
              'files/chemokine.xgmml',
              'files/mtor.xgmml',
              'files/toll_like.xgmml' ]

FILENAME_TO_PATHWAY_MAP = {
     'files/cytosolic_dna_sensing.xgmml' : "JAK/STAT",
     'files/nod_like.xgmml' : "Innate immune pathways",
     'files/rig1.xgmml' : "Innate immune pathways",
     'files/jakstat.xgmml'	: "JAK/STAT",
     'files/mapk.xgmml' : "Receptor tyrosine kinase",
     'files/chemokine.xgmml' : "Innate immune pathways",
     'files/mtor.xgmml' : "Receptor tyrosine kinase",
     'files/toll_like.xgmml' : "Toll-like receptor"
    }


def main(logger):
    # Getting pathway data
    ens_id = {}
    ids = {}

    for file in XGMML_LIST:
        for node in ET.parse(file).findall('.//{http://www.cs.rpi.edu/XGMML}node'):
            ids[node.get('id')] = FILENAME_TO_PATHWAY_MAP[file.replace('layer0/databases/innatedb/', '')]
            for att in node.findall('.//{http://www.cs.rpi.edu/XGMML}att'):
                if att.get('name') == 'Cross-references' and att.get('value') != '':
                    ens_id[(att.get('value')).split('|')[2]] = FILENAME_TO_PATHWAY_MAP[file.replace('layer0/databases/innatedb/', '')]

    def insert_or_get_node_dict(id, alt_id, alias, taxid, node_names_to_id, db_api, node):

        if id in ens_id:
            pathway = ens_id[id]
        elif alt_id.replace('innatedb:', '') in ids:
            pathway = ids[alt_id.replace('innatedb:', '')]
        else:
            pathway = None

        node_dict = {
            "name": node,
            "tax_id": taxid,
            "alt_accession": alt_id,
            'pathways': pathway,
            "aliases": alias,
            "topology": None
        }

        if node_dict['name'] in node_names_to_id:
            node_dict['id'] = node_names_to_id[node_dict['name']]
        else:
            db_api.insert_unique_node(node_dict)
            node_names_to_id[node_dict['name']] = node_dict['id']

        return node_dict

    # Initiating the parser
    db_api = PsimiSQL(SQL_SEED)

    # Uncomment if using run_all_auto script
    # data = DATA_FILE
    # next(data)

    with open(DATA_FILE) as data:
        # Skipping header
        data.readline()

        node_names_to_id = {}

        for line in data:
            line = str(line).split('\t')
            if line[9] == 'taxid:9606(Human)' and line[10] == 'taxid:9606(Human)':

                tax_id = "taxid:9606"
                source_node_id = line[2].split(":")[1]
                target_node_id = line[3].split(":")[1]
                source_node = f'Ensembl:{line[2].split(":")[1]}'
                target_node = f'Ensembl:{line[3].split(":")[1]}'

                # If sourcedb is biogrid, we get rid of it, it's interaction quality is not good
                if 'biogrid' in line[12].lower():
                    continue
                if 'intact' in line[12].lower():
                    continue
                if 'dip' in line[12].lower():
                    continue
                if 'transcompel' in line[12].lower():
                    continue
                if 'transfac' in line[12].lower():
                    continue
                if 'mint' in line[12].lower():
                    continue
                if 'bind' in line[12].lower():
                    continue

                # Creating the node dicts, if the node is already in the db assigning that to the node dict
                source_dict = insert_or_get_node_dict(source_node_id, line[0], line[4], tax_id, node_names_to_id, db_api, source_node)
                target_dict = insert_or_get_node_dict(target_node_id, line[1], line[5], tax_id, node_names_to_id, db_api, target_node)

                # Removing 'psimi:' from the front of the expression, removing quoutes
                intdetmethod = line[6].replace('psi-mi:', '').replace('"', '')

                # Interaction types
                # Removing 'psimi:' from the front of the expression, removing quoutes

                inttype_final = line[11].replace('psi-mi:', '').replace('"', '')
                interaction_types = "%s|is_directed:%s|is_direct:%s" \
                                    % (inttype_final, "true", "false")

                pubmed_ids = line[8].strip().split("|")
                if "pubmed:23180781" not in pubmed_ids:
                    pubmed_ids.append("pubmed:23180781")

                scores = line[14].strip().split("|")
                new_scores = []

                if line[14] != '':

                    for s in scores:
                        if s != '':
                            new_score = f'InnateDB {s}'
                            if new_score not in new_scores:
                                new_scores.append(new_score)

                edge_dict = {
                    'publication_ids': "|".join(pubmed_ids),
                    'layer': '0',
                    'source_db': 'InnateDB',  # ontology database citation
                    'interaction_identifiers': None,
                    'confidence_scores': "|".join(new_scores),  # if available
                    'interaction_detection_method': intdetmethod,
                    'interaction_types': interaction_types,
                    'first_author': None
                }

                db_api.insert_edge(source_dict, target_dict, edge_dict)

    # Saving the to a DB_TYPE.db files
    db_api.save_db_to_file(EXPORT_DB_LOCATION)


if __name__ == '__main__':
    print("Parsing database...")
    main(logger = None,)
    print("Parsing database is completed. SQLite database is saved to: " + EXPORT_DB_LOCATION)
