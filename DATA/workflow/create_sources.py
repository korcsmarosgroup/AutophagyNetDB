"""
Runs all SLK3 scripts
"""
# Imports
import logging
from collections import OrderedDict
import json


# Defining executing functions for each layer, where possible using auto download from the database website
def run_layer0(db, log, path):

    if db == 'acsn':
        pathway_file_loc = path + 'acsn_v1_pathway_file.gmt'
        cur_prot_loc = path + 'acsn_v2_cur_prot_list.txt'
        all_edge_loc = path + 'acsn_v1_all_edge.gmt'

        from DATA.workflow.SLK_Core.databases.acsn import script
        script.SQL_SEED = '../../ARNlib/SQLiteDBApi/network-db-seed.sql'
        script.PATHWAY_FILE_LOCATION = pathway_file_loc
        script.EXPORT_DB_LOCATION = 'all_output/acsn'
        script.CURATED_PROTEIN_LIST_FILE_LOCATION = cur_prot_loc
        script.ALL_EDGE_FILE_LOCATION = all_edge_loc

        return script.main(logger=log)

    elif db == 'innatedb':
        data_file_loc = path + 'innatedb_v5.4_datafile.mitab'

        from DATA.workflow.SLK_Core.databases.innatedb import script
        script.SQL_SEED = '../../ARNlib/SQLiteDBApi/network-db-seed.sql'
        script.DATA_FILE = data_file_loc
        script.EXPORT_DB_LOCATION = 'all_output/innatedb'
        script.XGMML_LIST = [path + 'cytosolic_dna_sensing.xgmml',
                             path + 'nod_like.xgmml',
                             path + 'rig1.xgmml',
                             path + 'jakstat.xgmml',
                             path + 'mapk.xgmml',
                             path + 'chemokine.xgmml',
                             path + 'mtor.xgmml',
                             path + 'toll_like.xgmml']

        return script.main(logger=log)

    elif db == 'reactome':
        data_file_loc = path + 'reactome_v3.3_datafile.txt'
        uni2path_loc = path + 'rectome_v3.3_uni2pathway.txt'

        from DATA.workflow.SLK_Core.databases.reactome import script
        script.SQL_SEED = '../../ARNlib/SQLiteDBApi/network-db-seed.sql'
        script.DATA_FILE = data_file_loc
        script.PATHWAY_FILE_LOCATION = path + 'pathway_map.txt'
        script.UNI_TO_PATHWAY = uni2path_loc
        script.EXPORT_DB_LOCATION = 'all_output/reactome'

        return script.main(logger=log)

    elif db == 'signor':

        from DATA.workflow.SLK_Core.databases.signor import script
        script.SQL_SEED = '../../ARNlib/SQLiteDBApi/network-db-seed.sql'

        script.DB_DESTINATION = 'all_output/signor'
        script.CSV_LIST = [path + 'SIGNOR_WNT.csv',
                           path + 'SIGNOR_TOLLR.csv',
                           path + 'SIGNOR_TGFb.csv',
                           path + 'SIGNOR_SAPK_JNK.csv',
                           path + 'SIGNOR_P38.csv',
                           path + 'SIGNOR_NOTCH.csv',
                           path + 'SIGNOR_NFKBNC.csv',
                           path + 'SIGNOR_NFKBC.csv',
                           path + 'SIGNOR_MTOR.csv',
                           path + 'SIGNOR_MCAPO.csv',
                           path + 'SIGNOR_INSR.csv',
                           path + 'SIGNOR_IL1R.csv',
                           path + 'SIGNOR_IAPO.csv',
                           path + 'SIGNOR_AMPK.csv',
                           path + 'SIGNOR_BMP.csv',
                           path + 'SIGNOR_DR.csv',
                           path + 'SIGNOR_EGF.csv',
                           path + 'SIGNOR_HPP.csv',
                           path + 'SIGNOR_Hedgehog.csv']
        return script.main(logger=log)

    # elif db == 'signafish':
    #     data_file_loc = path + 'signafish-L0.csv'
    #
    #     from DATA.workflow.SLK_Core.databases.signafish import script
    #     script.SQL_SEED = '../../ARNlib/SQLiteDBApi/network-db-seed.sql'
    #     script.DATA_FILE = data_file_loc
    #     script.EXPORT_DB_LOCATION = 'all_output/signafish'
    #
    #     return script.main(logger=log)

    elif db == 'slk2':
        from DATA.workflow.SLK_Core.databases.slk2 import script
        script.SQL_SEED = '../../ARNlib/SQLiteDBApi/network-db-seed.sql'

        script.FILE_LOCATION = path + '04122017-signalink-1uBbAI.csv'
        script.DB_DESTINATION = 'all_output/slk2'
        return script.main(logger=log)

    elif db == 'slk21':
        from DATA.workflow.SLK_Core.databases.slk21 import script
        script.SQL_SEED = '../../ARNlib/SQLiteDBApi/network-db-seed.sql'

        script.DB_DESTINATION = 'all_output/slk21'
        script.SLK_21_FILE_LOCATION = path + 'SLK21.txt'
        script.TAX_ID_MAP_FILE_LOCATION = path + 'entrytotaxid.txt'
        return script.main(logger=log)

    elif db == 'slk3':
        from DATA.workflow.SLK_Core.databases.slk3 import script
        script.SQL_SEED = '../../ARNlib/SQLiteDBApi/network-db-seed.sql'

        script.DESTINATION = 'all_output/slk3'
        script.TSV_LIST = [path + 'SLK3_update_cele_kitti_betti.txt',
                           path + 'SLK3_update_dmel_kitti_viktor.txt',
                           path + 'SLK3_update_hsap_kitti_d.txt']
        return script.main(logger=log)

    elif db == 'tcr':
        from DATA.workflow.SLK_Core.databases.tcr import script
        script.SQL_SEED = '../../ARNlib/SQLiteDBApi/network-db-seed.sql'

        script.TCR_DATA_LOC = path + 'T_cell_activation_signaling_pathway.txt'
        script.DESTINATION = 'all_output/tcr'
        return script.main(logger=log)


def run_layer1(db, log, path):
    if db == 'PSP':
        data_file_loc = path + 'psp_2009.11.003_datafile.txt'

        from DATA.workflow.layer1.databases.PSP import script
        script.SQL_SEED = '../../ARNlib/SQLiteDBApi/network-db-seed.sql'
        script.DATA_FILE = data_file_loc
        script.EXPORT_DB_LOCATION = 'all_output/PSP'

        return script.main(logger=log)

    elif db == 'SLK2_endo':
        data_file_loc = path + 'SLK2_L1_endocytosis.csv'

        from DATA.workflow.layer1.databases.SLK2_endo import script
        script.SQL_SEED = '../../ARNlib/SQLiteDBApi/network-db-seed.sql'
        script.DATA_FILE = data_file_loc
        script.EXPORT_DB_LOCATION = 'all_output/SLK2_endo'

        return script.main(logger=log)

    elif db == 'SLK2_scaffold':
        data_file_loc = path + 'SLK2_L1_scafford.csv'

        from DATA.workflow.layer1.databases.SLK2_scaffold import script
        script.SQL_SEED = '../../ARNlib/SQLiteDBApi/network-db-seed.sql'
        script.DATA_FILE = data_file_loc
        script.EXPORT_DB_LOCATION = 'all_output/SLK2_scaffold'

        return script.main(logger=log)

def run_PTM(db, log, path):
    if db == 'PhosphoSite':
        from DATA.workflow.PTM.databases.PhosphoSite import script
        script.SQL_SEED = '../../ARNlib/SQLiteDBApi/network-db-seed.sql'
        script.DATA_FILE = path + 'Kinase_Substrate_Dataset.txt'
        script.EXPORT_DB_LOCATION = 'all_output/PhosphoSite'

        return script.main(logger=log)

    # elif db == 'ELMpred':
    #     classes_loc = path + 'elm_27082018_classes.tsv'
    #     interactions_loc = path + 'elm_30012018_interactions.tsv'
    #
    #     from DATA.workflow.PTM.databases.ELMpred import pred_new as script
    #     script.SQL_SEED = '../../ARNlib/SQLiteDBApi/network-db-seed.sql'
    #     script.ELMS_FILE = classes_loc
    #     script.INT_DOMAINS_FILE = interactions_loc
    #     script.PFAM_FILE_LIST = [path + 'homo_pfam.txt',
    #                              path + 'drosi_pfam.txt',
    #                              path + 'celegans_pfam.txt',
    #                              path + 'danio_pfam.txt']
    #     script.PDB_LIST = [path + 'homo_pdb.txt',
    #                        path + 'drosi_pdb.txt',
    #                        path + 'danio_pdb.txt',
    #                        path + 'celegans_pdb.txt']
    #     script.file_list = []
    #     script.pred_db = 'PTM/databases/ELMpred/ELM_pred.db'
    #     script.EXPORT_DB_LOCATION = 'all_output/ELM'
    #
    #     return script.main(logger=log)

    if db == 'PTMCode2':
        from DATA.workflow.PTM.databases.PTMCode2 import script
        script.SQL_SEED = '../../ARNlib/SQLiteDBApi/network-db-seed.sql'
        script.DATA_FILE = path + 'PTMcode2_associations_between_proteins.txt'
        script.EXPORT_DB_LOCATION = 'all_output/PTMCode2'

        return script.main(logger=log)


def run_ATG_Reg(db, log, path):
    if db == 'chip_behrends':
        from DATA.workflow.ATG_Reg.databases.chip_behrends import chip_script
        chip_script.SQL_SEED = '../../ARNlib/SQLiteDBApi/network-db-seed.sql'
        chip_script.DATA_FILE = path + 'ARN_chip-ms_Behrends.csv'
        chip_script.EXPORT_DB_DESTINATION = 'all_output/chip_behrends'
        return chip_script.main(logger=log)

    # TODO finish import script
    # if db == 'HumanAutophagyDB':
    #     from DATA.workflow.ATG_Reg.databases.HumanAutophagyDB import script
    #     script.SQL_SEED = '../../ARNlib/SQLiteDBApi/network-db-seed.sql'
    #     script.HA_DATA_FILE = path + 'HAdb_mapped.txt'
    #     script.ATG_DATA_FILE = path + 'atg_mapped.txt'
    #     script.mergerdb = '../../merger.db'

    if db == 'manual_curation':
        from DATA.workflow.ATG_Reg.databases.manual_curation import script
        script.SQL_SEED =  '../../ARNlib/SQLiteDBApi/network-db-seed.sql'
        script.DATA_FILE_LIST = [path + 'files/Autofágia Regulatory Network - TD-nek_2013. 09. 26..txt',
                          path + 'files/Autofágia Regulatory Network - TD-nek_v2.txt']
        script.EXPORT_DB_DESTINATION = '.all_output/manualcur'

    if db == 'biogrid':
        from DATA.workflow.ATG_Reg.databases.biogrid import new_script
        new_script.SQL_SEED = '../../ARNlib/SQLiteDBApi/network-db-seed.sql'
        new_script.RAW_FILE_LIST = [path + 'BIOGRID-ALL-3.5.170.mitab.txt']
        new_script.EXPORT_DB_LOCATION = 'all_output/biogrid'
        return new_script.main(logger=log)

    elif db == 'ComPPI':
        from DATA.workflow.ATG_Reg.databases.ComPPI import script
        script.SQL_SEED = '../../ARNlib/SQLiteDBApi/network-db-seed.sql'
        script.DATA_FILE_LIST = [path + 'comppi--interactions--tax_celegans_loc_all.txt',
                                 path + 'comppi--interactions--tax_hsapiens_loc_all.txt',
                                 path + 'comppi--interactions--tax_dmelanogaster_loc_all.txt']
        script.EXPORT_DB_LOCATION = 'all_output/ComPPI'
        return script.main(logger=log)

    elif db == 'HPRD':
        from DATA.workflow.ATG_Reg.databases.HPRD import script
        script.SQL_SEED = '../../ARNlib/SQLiteDBApi/network-db-seed.sql'
        script.INTERACTION_FILE = path + 'HPRD_Release9_062910/BINARY_PROTEIN_PROTEIN_INTERACTIONS.txt'
        script.EXPORT_DB_LOCATION = 'all_output/HPRD'
        return script.main(logger=log)

    elif db == 'IntAct':
        from DATA.workflow.ATG_Reg.databases.IntAct import script
        script.SQL_SEED = '../../ARNlib/SQLiteDBApi/network-db-seed.sql'
        script.EXPORT_DB_LOCATION = 'all_output/IntAct'
        script.DATA_FILE = path + 'intact-micluster.txt'

        return script.main(logger=log)

    elif db == 'OmniPath':
        data_file_loc = path + 'omnipath_v0.7.111_datafile.txt'

        from DATA.workflow.ATG_Reg.databases.OmniPath import script
        script.SQL_SEED = '../../ARNlib/SQLiteDBApi/network-db-seed.sql'
        script.DATA_FILE = data_file_loc
        script.EXPORT_DB_LOCATION = 'all_output/OmniPath'

        return script.main(logger=log)


def run_miRNA(db, log, path):
    if db == 'miR2Disease':
        from DATA.workflow.miRNA.databases.miR2Disease import script
        script.SQL_SEED = '../../ARNlib/SQLiteDBApi/network-db-seed.sql'
        script.DATA_FILE = path + 'miRtar.txt'
        script.DB_DESTINATION = 'all_output/miR2Disease'

        return script.main(logger=log)

    elif db == 'miRDeathDB':
        from DATA.workflow.miRNA.databases.miRDeathDB import script
        script.SQL_SEED = '../../ARNlib/SQLiteDBApi/network-db-seed.sql'
        script.DATA_FILE = path + 'miRDeathDB_all_data.txt'
        script.DB_DESTINATION = 'all_output/miRDeathDB'
        return script.main(logger=log)

    elif db == 'miRecords':
        data_file_loc = path + 'mirecords_v4_datafile.txt'

        from DATA.workflow.miRNA.databases.miRecords import script
        script.SQL_SEED = '../../ARNlib/SQLiteDBApi/network-db-seed.sql'
        script.DATA_FILE = data_file_loc
        script.DB_DESTINATION = 'all_output/miRecords'

        return script.main(logger=log)

    elif db == 'TarBase':
        celegans_file_loc = path + 'tarbase_v7_celegans.txt'
        human_file_loc = path + 'tarbase_v7_human.txt'
        danio_file_loc = path + 'tarbase_v7_danio.txt'
        drosi_file_loc = path + 'tarbase_v7_drosi.txt'

        from DATA.workflow.miRNA.databases.TarBase import script
        script.SQL_SEED = '../../ARNlib/SQLiteDBApi/network-db-seed.sql'
        script.DATA_FILE_LIST = [celegans_file_loc,
                                 human_file_loc,
                                 danio_file_loc,
                                 drosi_file_loc]
        script.DB_DESTINATION = 'all_output/TarBase'

        return script.main(logger=log)

    elif db == 'starBase':
        from DATA.workflow.miRNA.databases.starBase import script
        script.SQL_SEED = '../../ARNlib/SQLiteDBApi/network-db-seed.sql'

        script.DATA_FILE_LIST = [path + 'starbase_v3_miRNAmRNA.txt',
                                 path + 'starbase_v3_degradome_human.txt',
                                 path + 'starbase_v3_degradome_worm.txt',
                                 path + 'starbase_v3_miRNA_valid.txt']
        script.FILE_TO_TAXID = {
            path + 'starbase_v3_miRNAmRNA.txt': "taxid:9606",
            path + 'starbase_v3_degradome_human.txt': 'taxid:9606',
            path + 'starbase_v3_degradome_worm.txt': 'taxid:6239',
            path + 'starbase_v3_miRNA_valid.txt': 'taxid:9606'
        }
        script.file_to_detmet = {
            path + 'starbase_v3_miRNAmRNA.txt': "MI:1110(predicted interaction)",
            path + 'starbase_v3_degradome_human.txt': 'MI:0045(experimental interaction detection)',
            path + 'starbase_v3_degradome_worm.txt': 'MI:0045(experimental interaction detection)',
            path + 'starbase_v3_miRNA_valid.txt': 'MI:2321(high throughput sequencing)'
        }
        script.DB_DESTINATION = 'all_output/starBase'

        return script.main(logger=log)

    elif db == 'miRSponge':
        data_file_loc = path + 'mirsponge_v1_datafile.txt'

        from DATA.workflow.miRNA.databases.miRSponge import script
        script.SQL_SEED = '../../ARNlib/SQLiteDBApi/network-db-seed.sql'
        script.DATA_FILE = data_file_loc
        script.DB_DESTINATION = 'all_output/miRSponge'
        return script.main(logger=log)


def run_TF(db, log, path):
    if db == 'TFBS':
        from DATA.workflow.TF.databases.TFBS import script
        script.SQL_SEED = '../../ARNlib/SQLiteDBApi/network-db-seed.sql'
        script.DATA_FILE = path + 'HITSallbodyresults.txt'
        script.EXPORT_DB_LOCATION = 'all_output/TFBS'

        return script.main(logger=log)

    elif db == 'TFBS_Orsi':
        from DATA.workflow.TF.databases.TFBS_Orsi import script2 as script
        script.SQL_SEED = '../../ARNlib/SQLiteDBApi/network-db-seed.sql'
        script.DATA_FILE = path + 'allInteractionData_v3_fixed.tsv'
        script.EXPORT_DB_LOCATION = 'all_output/TFBS_Orsi'

        return script.main(logger=log)


def run_lncRNA(db, log, path):
    if db == 'lncRInter':
        from DATA.workflow.lncRNA.databases.lncRInter import script
        script.SQL_SEED = '../../ARNlib/SQLiteDBApi/network-db-seed.sql'
        script.DATA_FILE = path + 'result_all.txt'
        script.DB_DESTINATION = 'all_output/lncRInter'
        return script.main(logger=log)


    elif db == 'starbase':
        from DATA.workflow.lncRNA.databases.starbase import script
        script.SQL_SEED = '../../ARNlib/SQLiteDBApi/network-db-seed.sql'
        script.DATA_FILE_LIST = [path + 'starbase_v3_miRNAlncRNA.txt',
                  path + 'starbase_v3_ncRNA_degradome_human.txt',
                  path + 'starbase_v3_ncRNA_degradome_worm.txt',
                  path + 'starbase_v3_lncRNA_valid.txt']
        script.FILE_TO_TAXID = {
            path + 'starbase_v3_miRNAlncRNA.txt': "taxid:9606",
            path + 'starbase_v3_ncRNA_degradome_human.txt': 'taxid:9606',
            path + 'starbase_v3_ncRNA_degradome_worm.txt': 'taxid:6239',
            path + 'starbase_v3_lncRNA_valid.txt': 'taxid:9606'
        }
        script.file_to_detmet = {
            path + 'starbase_v3_miRNAlncRNA.txt': "MI:1110(predicted interaction)",
            path + 'starbase_v3_ncRNA_degradome_human.txt': 'MI:0045(experimental interaction detection)',
            path + 'starbase_v3_ncRNA_degradome_worm.txt': 'MI:0045(experimental interaction detection)',
            path + 'starbase_v3_lncRNA_valid.txt': 'MI:2321(high throughput sequencing)'
        }
        script.DB_DESTINATION = 'all_output/starbase'

        return script.main(logger=log)

    elif db == 'NPInter':
        data_file_loc = path + 'interaction_NPInter[v3.0].txt'

        from DATA.workflow.lncRNA.databases.NPInter import script
        script.SQL_SEED = '../../ARNlib/SQLiteDBApi/network-db-seed.sql'
        script.DATA_FILE = data_file_loc
        script.DB_DESTINATION = 'all_output/NPInter'

        return script.main(logger=log)


DB_DICT = json.load(open('test_sources.json'), object_pairs_hook=OrderedDict)

# Initiating logger
logger = logging.getLogger()
handler = logging.FileHandler('SLK3.log')
logger.setLevel(logging.DEBUG)
handler.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)

# Running all scripts
for layer in DB_DICT.keys():
    for db in DB_DICT[layer]:
        logger.debug("Started parsing %s" % db)
        print("Started parsing %s" % db)
        runpath = layer + '/databases/' + db + '/files/'

        if layer == 'SLK_Core':
            run_layer0(db, log=logger, path=runpath)
        elif layer == 'layer1':
            run_layer1(db, log=logger, path=runpath)
        elif layer == 'PTM':
            run_PTM(db, log=logger, path=runpath)
        elif layer == 'ATG_Reg':
            run_ATG_Reg(db, log=logger, path=runpath)
        if layer == 'miRNA':
            run_miRNA(db, log=logger, path=runpath)
        elif layer == 'TF':
            run_TF(db, log=logger, path=runpath)
        elif layer == 'lncRNA':
            run_lncRNA(db, log=logger, path=runpath)


