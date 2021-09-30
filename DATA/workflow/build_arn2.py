from collections import OrderedDict
import json
import logging
import pandas
import sqlite3


from ARNlib.mapper.protein import create_mapping_db
from ARNlib.mapper.protein import create_mapping_db_casesense
from ARNlib.mapper.protein import molecular_id_mapper as mol
from DATA.workflow.lncRNA.mapper import lncmap_new
from ARNlib.merger import merge_layer
from ARNlib import build_new as builder
from ARNlib import noconn_check as noconn
from ARNlib import sort_data_arn as sorter

# predictions
from DATA.prediction.tissue.pred_script_new import tissue_prediction
from DATA.prediction.subcell.script import subcell_prediction
from DATA.prediction.direction.pred_script import DirScore
from DATA.prediction.RNAipred.pred_script import main
from DATA.prediction.drugtarget.DrugBank.script import drugbank

DB_DICT = json.load(open('sources.json'), object_pairs_hook=OrderedDict)

# Initiating logger
logger = logging.getLogger()
handler = logging.FileHandler('ARN2_build.log')
logger.setLevel(logging.DEBUG)
handler.setLevel(logging.DEBUG)
formatter = logging.Formatter(
    '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)

# # Creating mapping db
# # For proteins
# logger.debug('Creating mapping dbs')
# MDB = create_mapping_db.CreateMappingDB(mappingDBfile='mapper.db', debug=False)
# MDB.add_species('9606')

#
# # For miRNAs
# lncmap_new.SQL_SEED = '../../ARNlib/SQLiteDBApi/network-db-seed.sql'
# lncmap_new.noncode_file = 'lncRNA/mapper/noncode.tsv'
# lncmap_new.ens_file = 'lncRNA/mapper/ensembl.tsv'
# lncmap_new.mir_file = 'lncRNA/mapper/mirbase.dat'
# lncmap_new.mir_file_2 = 'lncRNA/mapper/mirbase.tsv'
# lncmap_new.lncipedia_file = 'lncRNA/mapper/lncipedia.tsv'
# lncmap_new.hgnc_file = 'lncRNA/mapper/hgnc.tsv'
# lncmap_new.DB_DESTINATION = 'lncrna_map.db'
#
# lncmap_new.main(logger=logger)

# # Mapping
# mapping_results = []
# for layer in DB_DICT.keys():
#     for mol.db in DB_DICT[layer]:
#         logger.debug("Started mapping %s" % mol.db)
#         mol.m = mol.MolecularIDMapper(mol.db, PROT_DBname='mapper.db',
#                                       LNCRNAMAP_DBname='lncrna_map.db',
#                                       layer=layer)
#         mol.m.main()
#         # (db_name, total_edges, failed_edges) = mol.m.main()
#         # mapping_results.append((db_name, total_edges, failed_edges))
# logger.debug('Mapping done')
# print("\n\n================= mapping statistics ================")
# for db_name, total_edges, failed_edges in mapping_results:
#     rate = "{0:.2f}%".format(round(100.0 * (failed_edges / total_edges),2))
#     print("db: %s     total edges: %d    edges failed to map: %d    failure rate: %s" % (db_name, total_edges, failed_edges, rate))
# print("\n\n")

# Merge layers
merge_layer.SQL_SEED_LOCATION = '../../ARNlib/SQLiteDBApi/network-db-seed.sql'
merge_layer.DESTINATION = 'merger'
path = '../../ARNlib/mapper/protein/output/'
dblist = []
for layer in DB_DICT.keys():
    for db in DB_DICT[layer]:
        dblist.append(path + db + '_mapped.db')
merge_layer.SOURCE_DB_FILE_LIST = dblist
logger.debug('Started merging')
merge_layer.main(log=logger)
logger.debug('Merging done')

# Deleting nodes from merger that have no connections
logger.debug('Deleting noconn nodes')
noconn.main(logger=logger, merger_path='merger.db')

# Building all layers
logger.debug('Started building')
builder.main(log=logger, path='merger.db')
logger.debug('Building done')

# Predictions
# Tissue
tissue_prediction(logger)
# Direction
test = DirScore()
logger.debug('Creating test set')
test.test_scores()
logger.debug('Adding scores to dataset')
test.apply_to_db()
logger.debug('Direction prediction done')
# RNAi
main(logger)
# Drug-target
drugbank(logger)
# Subcellular localization
subcell_prediction(logger)

# Creating case sensitive mapping db for gene name and protein name mapping
# For proteins
logger.debug('Creating mapping dbs')
MDB = create_mapping_db_casesense.CreateMappingDB(mappingDBfile='casesense_mapper.db', debug=False)
MDB.add_species('9606')

# Sorting data into json format
logger.debug('Started sorting')
sorter.mapper_location = 'mapper.db'
sorter.json_mapper_file = 'uniprot_id_mapping.json'
sorter.logger_output_location = 'SLK3_build.log'
sorter.get_node_data(builder_location='ARN2_layers.db')
sorter.get_edge_data(builder_location='ARN2_layers.db')
sorter.get_attribute_data(builder_location='ARN2_layers.db')
logger.debug('Sorting done')
