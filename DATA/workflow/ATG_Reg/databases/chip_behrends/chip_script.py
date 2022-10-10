"""
  Parses data collected from article: https://www.nature.com/articles/nature09204
"""

# Imports
import csv, logging
from ARNlib.SQLiteDBApi.sqlite_db_api import PsimiSQL

# Constants
SQL_SEED = '../../../../ARNlib/SQLiteDBApi/network-db-seed.sql'
DATA_FILE = 'ARN_chip-ms_Behrends.csv'
EXPORT_DB_DESTINATION = '../output/chip_behrends'
DB_TYPE = 'manual curation'

# Initiating logger
logger = logging.getLogger()
handler = logging.FileHandler('../../SLK3.log')
logger.setLevel(logging.DEBUG)
handler.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)


def main(logger):
    def get_node_a(id, taxid, pathway, alias, topology, psi_mi_to_sql_object):
        """
        This function sets up a node dict and returns it.
        If the node is already in the SQLite database it fetches that node from the db, so it won't be inserted multiple times.
        """

        # Testing if the node is already in the database

        node_dict = psi_mi_to_sql_object.get_node(id, node_tax_id=taxid)

        if not node_dict:
            node_dict = {
                "name": id,
                "tax_id": taxid,
                "alt_accession": None,
                'pathways': pathway,
                "aliases": alias,
                "topology": topology
            }

        return node_dict

    def get_node_b(id, taxid, pathway, alias, topology, psi_mi_to_sql_object):
        """
        This function sets up a node dict and returns it. If the node is already in the SQLite database it fetches that node from the db, so it won't be inserted multiple times.

        """

        # Testing if the node is already in the database
        node_dict = psi_mi_to_sql_object.get_node(id, node_tax_id=taxid)

        if not node_dict:
            node_dict = {
                "name": id,
                "tax_id": taxid,
                "alt_accession": None,
                'pathways': pathway,
                "aliases": alias,
                "topology": topology
            }

        return node_dict

    # Initiating the parser
    db_api = PsimiSQL(SQL_SEED)

    # Parsing data file
    with open(DATA_FILE) as data:
        # Skipping the header
        data.readline()

        for line in data:
            line = line.strip().split(';')
            # Taxid
            if line[2] == '9606':
                taxid_source = 'taxid:9606'
            else:
                taxid_source = line[2]
            if line[10] == '9606':
                taxid_target = 'taxid:9606'
            else:
                taxid_target = line[10]

            # Pathways
            source_ptw_list = []
            source_ptw_line = line[7].split(',')
            for ptw in source_ptw_line:
                ptw_new = ptw.strip().split('(')[0]
                source_ptw_list.append(ptw_new)

            source_ptw = '|'.join(source_ptw_list)

            target_ptw_list = []
            target_ptw_line = line[15].split(',')
            for ptw in target_ptw_line:
                ptw_new = ptw.strip().split('(')[0]
                target_ptw_list.append(ptw_new)

            target_ptw = '|'.join(target_ptw_list)

            # Topology
            source_topol = '|'.join(line[4].strip().split(','))

            target_topol = '|'.join(line[12].strip().split(','))



            # Creating the node dicts, if the node is already in the db assigning that to the node dict
            source_dict = get_node_a('Uniprot:' + line[1], taxid_source, source_ptw, line[0], source_topol, db_api)
            target_dict = get_node_b('Uniprot:' + line[9], taxid_target, target_ptw, line[8], target_topol, db_api)

            # Nodes are inserted to the db if they are not in it yet
            if not 'id' in source_dict:
                db_api.insert_node(source_dict)

            if not 'id' in target_dict:
                db_api.insert_node(target_dict)

            # Mapping layer descriptions to abbreviations
            layer_dict = {
                'Post-translational regulation': '2',
                'Interaction between autophagy proteins': '0',
                'Autophagy regulators': '1'
            }

            # Is directed
            directed_map = {
                'PPI directed': 'true',
                'PPI undirected': 'false'
            }

            # Is direct
            direct_map = {
                'direct': 'true'
            }
            is_direct = direct_map[line[18]]

            # Effect
            effect_map = {
                'stimulation': 'MI:0624(stimulation)'
            }

            if line[19] != 'unknown':
                effect = effect_map[line[19]]
                # Constructing interaction data line
                int_types = '|'.join([effect, 'is_directed:' + directed_map[line[17]],
                                      'is_direct:' + is_direct])
            else:
                # Constructing interaction data line
                int_types = '|'.join(['is_directed:' + directed_map[line[17]],
                                      'is_direct:' + is_direct])

            # Publications
            pubs = '|pubmed:'.join(line[20].split('|'))

            # Sourcedb mapping
            sourcedb_map = {
                'BioGRID': 'TheBiogrid',
                'Behrends et Al. 2010': 'Behrends',
                'direction is predicted': 'Behrends predicted'
            }
            dblist = []
            for db in line[21].split(','):
                sourcedb = db.strip().split('(')[0]
                if 'pmid' not in sourcedb:
                    if sourcedb in sourcedb_map.keys():
                        mysourcedb = sourcedb_map[sourcedb]
                    else:
                        mysourcedb = sourcedb

                    dblist.append(mysourcedb)

            final_source = '|'.join(dblist)
            if 'is_directed:false' in int_types:
                edge_dict = {
                    'publication_ids': 'pubmed:' + pubs,
                    'layer': '2',
                    'source_db': final_source,
                    'interaction_identifiers': None,
                    'confidence_scores': None,
                    'interaction_detection_method': None,
                    'interaction_types': int_types,
                    'first_author': None
                }
            else:
                edge_dict = {
                    'publication_ids': 'pubmed:' + pubs,
                    'layer': layer_dict[line[16]],
                    'source_db': final_source,
                    'interaction_identifiers': None,
                    'confidence_scores': None,
                    'interaction_detection_method': None,
                    'interaction_types': int_types,
                    'first_author': None
                }

            db_api.insert_edge(source_dict, target_dict, edge_dict)

            # Saving the to a DB_TYPE.db file
        db_api.save_db_to_file(EXPORT_DB_DESTINATION)


if __name__ == '__main__':
    print("Parsing database...")
    main(logger=None)
    print("Parsing database is completed. SQLite database is saved to: " + EXPORT_DB_DESTINATION)
