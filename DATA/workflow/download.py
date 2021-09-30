"""
Runs all SLK3 scripts
"""
# Imports
import logging, requests, gzip, xlrd, json, os, zipfile, tarfile
from io import BytesIO, StringIO
from collections import OrderedDict


# Defining executing functions for each layer, where possible using auto download from the database website
def dwld_layer0(db, log, path):
    if db == 'acsn':

        # Assigning output files
        pathway_file_loc = path + 'acsn_v1_pathway_file.gmt'
        cur_prot_loc = path + 'acsn_v2_cur_prot_list.txt'
        all_edge_loc = path + 'acsn_v1_all_edge.gmt'

        # Removing files if they exist
        if os.path.exists(pathway_file_loc):
            os.remove(pathway_file_loc)
        if os.path.exists(cur_prot_loc):
            os.remove(cur_prot_loc)
        if os.path.exists(all_edge_loc):
            os.remove(all_edge_loc)

        # Download data and write to output files
        with open(pathway_file_loc, 'w') as pathway_file:
            f = requests.get('https://acsn.curie.fr/files/acsn_master_curated.gmt')
            pathway_file.write(f.text)
        with open(cur_prot_loc, 'w') as cur_prot:
            f = requests.get('https://acsn.curie.fr/ACSN2/downloads/ACSN2_binary_relations_between_proteins_with_PMID.txt')
            cur_prot.write(f.text)
        with open(all_edge_loc, 'w') as all_edge:
            f = requests.get('https://acsn.curie.fr/files/acsn_names.gmt')
            all_edge.write(f.text)

    elif db == 'innatedb':

        # Assigning output file
        data_file_loc = path + 'innatedb_v5.4_datafile.mitab'

        # Removing file if exists
        if os.path.exists(data_file_loc):
            os.remove(data_file_loc)

        # Download data and write to output files
        with open(data_file_loc, 'w') as data_file:
            f = requests.get('https://www.innatedb.com/download/interactions/all.mitab.gz', stream=True)
            writedata = gzip.open(BytesIO(f.content), 'rt')
            data_file.write(writedata.read())

    elif db == 'reactome':

        # Assigning output file
        data_file_loc = path + 'reactome_v3.3_datafile.txt'
        uni2path_loc = path + 'rectome_v3.3_uni2pathway.txt'

        # Removing file if exists
        if os.path.exists(data_file_loc):
            os.remove(data_file_loc)
        if os.path.exists(uni2path_loc):
            os.remove(uni2path_loc)

        # Download data and write to output files
        with open(data_file_loc, 'w') as data_file:
            f = requests.get(
                'https://reactome.org/download/current/interactors/reactome.all_species.interactions.psi-mitab.txt')
            writefile = f.text
            data_file.write(writefile)
        with open(uni2path_loc, 'w') as uni2pathfile:
            f = requests.get('https://reactome.org/download/current/UniProt2Reactome_All_Levels.txt')
            writefile = f.text
            uni2pathfile.write(writefile)


def dwld_layer1(db, log, path):
    if db == 'PSP':

        # Assigning output file
        data_file_loc = path + 'psp_2009.11.003_datafile.txt'

        # Removing file if exists
        if os.path.exists(data_file_loc):
            os.remove(data_file_loc)

        # Download data and write to output files
        with open(data_file_loc, 'w') as datafile:
            f = requests.get('https://www.cell.com/cms/10.1016/j.tcb.2009.11.003/attachment/a551f0e2-ab8e-4932-85e3-1a9b33aa5c81/mmc2.xls')
            myfile = xlrd.open_workbook(file_contents=f.content)
            mysheet = myfile.sheet_by_index(0)
            # Write the rows
            for rownum in range(mysheet.nrows):
                columns = mysheet.row_values(rownum)
                line = []
                for elem in columns:
                    line.append(str(elem))
                datafile.write('\t'.join(line) + '\n')


def dwld_PTM(db, log, path):

    if db == 'ELMpred':
        # Assigning output file
        classes_loc = path + 'elm_27082018_classes.tsv'
        interactions_loc = path + 'elm_30012018_interactions.tsv'

        # Removing file if exists
        if os.path.exists(classes_loc):
            os.remove(classes_loc)
        if os.path.exists(interactions_loc):
            os.remove(interactions_loc)

        # Download data and write to output files
        with open(classes_loc, 'w') as classes:
            f = requests.get('http://elm.eu.org/elms/elms_index.tsv')
            classes.write(f.text)
        with open(interactions_loc, 'w') as interactions:
            f = requests.get('http://elm.eu.org/infos/browse_elm_interactiondomains.tsv')
            interactions.write(f.text)


def dwld_ATG_Reg(db, log, path):

    if db == 'biogrid':
        # Assigning output file
        data_file_loc = path + 'biogrid_v3.5.170_datafile.mitab'

        # Removing file if exists
        if os.path.exists(data_file_loc):
            os.remove(data_file_loc)

        # Download data and write to output files
        with open(data_file_loc, 'w') as datafile:
            f = requests.get('https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-3.5.170/BIOGRID-ALL-3.5.170.mitab.zip')
            with zipfile.ZipFile(BytesIO(f.content)) as zip_ref:
                zip_ref.extractall(path)

    elif db == 'OmniPath':
        # Assigning output file
        data_file_loc = path + 'omnipath_v0.7.111_datafile.txt'

        # Removing file if exists
        if os.path.exists(data_file_loc):
            os.remove(data_file_loc)

        # Download data and write to output files
        with open(data_file_loc, 'w') as datafile:
            f = requests.get('http://omnipathdb.org/interactions/?fields=sources&fields=references')
            datafile.write(f.text)


def dwld_miRNA(db, log, path):

    if db == 'miRecords':
        # Assigning output file
        data_file_loc = path + 'mirecords_v4_datafile.txt'

        # Removing file if exists
        if os.path.exists(data_file_loc):
            os.remove(data_file_loc)

        # Download data and write to output files
        with open(data_file_loc, 'w') as datafile:
            f = requests.get('http://c1.accurascience.com/miRecords/download_data.php?v=4')
            myfile = xlrd.open_workbook(file_contents=f.content)
            mysheet = myfile.sheet_by_index(0)
            # Write the rows
            for rownum in range(mysheet.nrows):
                columns = mysheet.row_values(rownum)
                line = []
                for elem in columns:
                    line.append(str(elem))
                datafile.write('\t'.join(line) + '\n')

    elif db == 'starBase':
        # Assigning output file
        miRNA_mRNA_data = path + 'starbase_v3_miRNAmRNA.txt'
        degradome_miRNA_mRNA_human = path + 'starbase_v3_degradome_human.txt'
        degradome_miRNA_mRNA_worm = path + 'starbase_v3_degradome_worm.txt'
        valid_miRNA_RNA_human = path + 'starbase_v3_miRNA_valid.txt'

        # Removing file if exists
        if os.path.exists(miRNA_mRNA_data):
            os.remove(miRNA_mRNA_data)
        if os.path.exists(degradome_miRNA_mRNA_human):
            os.remove(degradome_miRNA_mRNA_human)
        if os.path.exists(degradome_miRNA_mRNA_worm):
            os.remove(degradome_miRNA_mRNA_worm)
        if os.path.exists(valid_miRNA_RNA_human):
            os.remove(valid_miRNA_RNA_human)

        # Download data and write to output files
        with open(miRNA_mRNA_data, 'w') as datafile:
            f = requests.get('http://starbase.sysu.edu.cn/moduleDownload.php?source=agoClipRNA&type=txt&value=hg19;mRNA;all;1;0;0;1;None;PDCD4')
            datafile.write(f.text)

        with open(degradome_miRNA_mRNA_human, 'w') as degradome_human:
            f = requests.get('http://starbase.sysu.edu.cn/moduleDownload.php?source=degradomeRNA&type=txt&value=hg19;mRNA;hsa-let-7a-2-3p;1;all')
            degradome_human.write(f.text)
        with open(degradome_miRNA_mRNA_worm, 'w') as degradome_worm:
            f = requests.get(
                'http://starbase.sysu.edu.cn/moduleDownload.php?source=degradomeRNA&type=txt&value=ce11;mRNA;cel-let-7-3p;1;all')
            degradome_worm.write(f.text)
        with open(valid_miRNA_RNA_human, 'w') as valid:
            f = requests.get('http://starbase.sysu.edu.cn/moduleDownload.php?source=rnaRNA&type=txt&value=hg19;hsa-let-7a-5p;1;1')
            valid.write(f.text)

    elif db == 'TarBase':
        # Assigning output file
        celegans_file_loc = path + 'tarbase_v7_celegans.txt'
        human_file_loc = path + 'tarbase_v7_human.txt'
        danio_file_loc = path + 'tarbase_v7_danio.txt'
        drosi_file_loc = path + 'tarbase_v7_drosi.txt'

        # Removing file if exists
        if os.path.exists(celegans_file_loc):
            os.remove(celegans_file_loc)
        if os.path.exists(human_file_loc):
            os.remove(human_file_loc)
        if os.path.exists(danio_file_loc):
            os.remove(danio_file_loc)
        if os.path.exists(drosi_file_loc):
            os.remove(drosi_file_loc)

        # Download data and write to output files
        with open(celegans_file_loc, 'w') as celegansfile:
            celegans = requests.get('http://mirtarbase.mbc.nctu.edu.tw/cache/download/7.0/cel_MTI.xls')
            myfile = xlrd.open_workbook(file_contents=celegans.content)
            mysheet = myfile.sheet_by_index(0)
            # Write the rows
            for rownum in range(mysheet.nrows):
                columns = mysheet.row_values(rownum)
                line = []
                for elem in columns:
                    line.append(str(elem))
                celegansfile.write('\t'.join(line) + '\n')
        with open(human_file_loc, 'w') as humanfile:
            human = requests.get('http://mirtarbase.mbc.nctu.edu.tw/cache/download/7.0/hsa_MTI.xlsx')
            myfile = xlrd.open_workbook(file_contents=human.content)
            mysheet = myfile.sheet_by_index(0)
            # Write the rows
            for rownum in range(mysheet.nrows):
                columns = mysheet.row_values(rownum)
                line = []
                for elem in columns:
                    line.append(str(elem))
                humanfile.write('\t'.join(line) + '\n')
        with open(danio_file_loc, 'w') as daniofile:
            danio = requests.get('http://mirtarbase.mbc.nctu.edu.tw/cache/download/7.0/dre_MTI.xls')
            myfile = xlrd.open_workbook(file_contents=danio.content)
            mysheet = myfile.sheet_by_index(0)
            # Write the rows
            for rownum in range(mysheet.nrows):
                columns = mysheet.row_values(rownum)
                line = []
                for elem in columns:
                    line.append(str(elem))
                daniofile.write('\t'.join(line) + '\n')
        with open(drosi_file_loc, 'w') as drosifile:
            drosi = requests.get('http://mirtarbase.mbc.nctu.edu.tw/cache/download/7.0/dme_MTI.xls')
            myfile = xlrd.open_workbook(file_contents=drosi.content)
            mysheet = myfile.sheet_by_index(0)
            # Write the rows
            for rownum in range(mysheet.nrows):
                columns = mysheet.row_values(rownum)
                line = []
                for elem in columns:
                    line.append(str(elem))
                drosifile.write('\t'.join(line) + '\n')


def dwld_lncRNA(db, log, path):

    if db == 'NPInter':
        # Assigning output file
        data_file_loc = path + 'npinter_v3_datafile.txt'

        # Removing file if exists
        if os.path.exists(data_file_loc):
            os.remove(data_file_loc)

        # Download data and write to output files
        with open(data_file_loc, 'w') as datafile:

            f = requests.get('https://www.bioinfo.org/NPInter/datadownload/interaction_NPInter[v3.0].txt.zip', stream=True)
            with zipfile.ZipFile(BytesIO(f.content)) as zip_ref:
                zip_ref.extractall(path)

    elif db == 'starbase':
        # Assigning output file
        miRNA_lncRNA_data = path + 'starbase_v3_miRNAlncRNA.txt'
        degradome_miRNA_ncRNA_human = path + 'starbase_v3_ncRNA_degradome_human.txt'
        degradome_miRNA_ncRNA_worm = path + 'starbase_v3_ncRNA_degradome_worm.txt'
        valid_lncRNA_RNA_human = path + 'starbase_v3_lncRNA_valid.txt'

        # Removing file if exists
        if os.path.exists(miRNA_lncRNA_data):
            os.remove(miRNA_lncRNA_data)
        if os.path.exists(degradome_miRNA_ncRNA_human):
            os.remove(degradome_miRNA_ncRNA_human)
        if os.path.exists(degradome_miRNA_ncRNA_worm):
            os.remove(degradome_miRNA_ncRNA_worm)
        if os.path.exists(valid_lncRNA_RNA_human):
            os.remove(valid_lncRNA_RNA_human)

        # Download data and write to output files
        with open(miRNA_lncRNA_data, 'w') as datafile:
            f = requests.get('http://starbase.sysu.edu.cn/moduleDownload.php?source=agoClipRNA&type=txt&value=hg19;lncRNA;all;1;0;0;1;None;MALAT1')
            datafile.write(f.text)

        with open(degradome_miRNA_ncRNA_human, 'w') as degradome_human:
            f = requests.get('http://starbase.sysu.edu.cn/moduleDownload.php?source=degradomeRNA&type=txt&value=hg19;ncRNA;hsa-let-7a-2-3p;1;all')
            degradome_human.write(f.text)
        with open(degradome_miRNA_ncRNA_worm, 'w') as degradome_worm:
            f = requests.get(
                'http://starbase.sysu.edu.cn/moduleDownload.php?source=degradomeRNA&type=txt&value=ce11;ncRNA;cel-let-7-3p;1;all')
            degradome_worm.write(f.text)
        with open(valid_lncRNA_RNA_human, 'w') as valid:
            f = requests.get('http://starbase.sysu.edu.cn/moduleDownload.php?source=rnaRNA&type=txt&value=hg19;H19;1;1')
            valid.write(f.text)


DB_DICT = json.load(open('test_sources.json'), object_pairs_hook=OrderedDict)

#%% Initiating logger
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
        print(db)
        logger.debug("Started downloading %s" % db)
        dwldpath = layer + '/databases/' + db + '/files/'

        if layer == 'layer0':
            dwld_layer0(db, log=logger, path=dwldpath)
        elif layer == 'layer1':
            dwld_layer1(db, log=logger, path=dwldpath)
        elif layer == 'PTM':
            dwld_PTM(db, log=logger, path=dwldpath)
        elif layer == 'ATG_Reg':
            dwld_ATG_Reg(db, log=logger, path=dwldpath)
        if layer == 'miRNA':
            dwld_miRNA(db, log=logger, path=dwldpath)
        elif layer == 'lncRNA':
            dwld_lncRNA(db, log=logger, path=dwldpath)

logger.debug('Downloading done.')
