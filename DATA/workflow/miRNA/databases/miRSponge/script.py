'''
    miRNA-lncRNA interactions
    :argument: DATA_FILE: http://www.bio-bigdata.com/miRSponge/apps/download1.jsp?fileName=Experimentally_validated_miRNA_targets_in_miRSponge.txt
'''

# Imports
from SLKlib.SQLiteDBApi.sqlite_db_api import PsimiSQL

# Defining constants
SQL_SEED = '../../../../../SLKlib/SQLiteDBApi/network-db-seed.sql'
DATA_FILE = 'files/Experimentally_validated_miRNA_targets_in_miRSponge.txt'
DB_DESTINATION = '../../output/miRSponge'


def get_node_lncrna(name, id, taxid, psi_mi_to_sql_object):
    """
    This function sets up a node dict and returns it. If the node is already in the SQLite database it fetches that node from the db, so it won't be inserted multiple times.
    """

    # Testing if the node is already in the database
    node_dict = psi_mi_to_sql_object.get_node(name, node_tax_id=taxid)

    if not node_dict:
        node_dict = {
            "name" : name,
            "tax_id": taxid,
            "alt_accession": id,
            'pathways': None,
            "aliases": None,
            "topology": None
        }

    return node_dict


def get_node_mirna(name, taxid, psi_mi_to_sql_object):
    """
    This function sets up a node dict and returns it. If the node is already in the SQLite database it fetches that node from the db, so it won't be inserted multiple times.
    """

    # Testing if the node is already in the database
    node_dict = psi_mi_to_sql_object.get_node(name, node_tax_id=taxid)

    if not node_dict:
        node_dict = {
            "name": 'miR:hsa-' + name,
            "tax_id": taxid,
            "alt_accession": None,
            'pathways': None,
            "aliases": None,
            "topology": None
        }

    return node_dict

def main(logger):
    # Declaring variables and constants
    inserted_nodes = {}

    # Initiating the parser
    db_api = PsimiSQL(SQL_SEED)
    with open(DATA_FILE, encoding='ISO-8859-1') as data:
        # Skipping the header
        data.readline()
    # data = DATA_FILE.split('\n')
        for line in data:
            columns = line.strip().split('\t')
            taxid = 'taxid:9606'
            if len(columns) > 4:

                mirna_name = '-'.join([columns[0].split('-')[0], columns[0].split('-')[1]])
                # Creating the node dicts, if the node is already in the db assigning that to the node dict
                source_dict = get_node_mirna(mirna_name, taxid, db_api)
                target_dict = get_node_lncrna('entrez:' + columns[1], columns[2], taxid, db_api)

                # Nodes are inserted to the db if they are not in it yet
                if not 'id' in source_dict:
                    db_api.insert_node(source_dict)

                if not 'id' in target_dict:
                    db_api.insert_node(target_dict)

                interaction_types = "effect:%s|is_directed:%s|is_direct:%s" \
                                    % ('MI:0256(rna interference)', 'directed', 'unknown')
                detmap = {
                    'qRT-PCR': 'MI:1196(quantitative reverse transcription pcr)',
                    'Luciferase reporter assay': 'MI:2285(miRNA interference luciferase reporter assay)',
                    'Western blot': 'MI:0113(western blot)',
                    'GFP reporter assay': 'MI:0045(experimental interaction detection)',
                    'In situ hybridization': 'MI:0045(experimental interaction detection)',
                    'Northern blot': 'MI:0929(northern blot)',
                    'Reporter assay': 'MI:0045(experimental interaction detection)',
                    'Other': 'MI:0045(experimental interaction detection)',
                    'Microarray': 'MI:0008(array technology)',
                    'Immunohistochemistry': 'MI:1198(immunohistochemistry)',
                    'Immunocytochemistry': 'MI:1200(immunocytochemistry)',
                    'Immunoblot': 'MI:0045(experimental interaction detection)',
                    '5RACE': 'MI:0045(experimental interaction detection)',
                    'phenotypic sensor assay': 'MI:0045(experimental interaction detection)',
                    'real-time RT-PCR': 'MI:1196(quantitative reverse transcription pcr)',
                    'in situ hybridization': 'MI:0045(experimental interaction detection)',
                    'FACS': 'MI:0045(experimental interaction detection)',
                    'ELISA': 'MI:0045(experimental interaction detection)',
                    'Flow': 'MI:0045(experimental interaction detection)',
                    'ChIP-seq': 'MI:0402(chromatin immunoprecipitation assay)',
                    'Immunofluorescence': 'MI:0045(experimental interaction detection)',
                    'GFP Reporter Assay': 'MI:0045(experimental interaction detection)',
                    'HITS-CLIP': 'MI:2191(clip)',
                    'PAR-CLIP': 'MI:2191(clip)',
                    'intrarenal expression': 'MI:0045(experimental interaction detection)',
                    'Proteomics': 'MI:0045(experimental interaction detection)',
                    'ChIP immunoprecipitation': 'MI:0402(chromatin immunoprecipitation assay)',
                    'Luciferase assay': 'MI:2285(miRNA interference luciferase reporter assay)',
                    'QRTPCR': 'MI:1196(quantitative reverse transcription pcr)',
                    'Next Generation Sequencing (NGS)': 'MI:0078(nucleotide sequence identification)',
                    'RNA-binding protein immunoprecipitation': 'MI:1017(rna immunoprecipitation)',
                    'immunohistochemistry': 'MI:1198(immunohistochemistry)',
                    'Sequencing': 'MI:0078(nucleotide sequence identification)',
                    'CLASH': 'MI:2195(clash)',
                    'immunoprecipitaion': 'MI:1017(rna immunoprecipitation)',
                    'Quantitative proteomic approach': 'MI:0045(experimental interaction detection)',
                    'ChIP': 'MI:0402(chromatin immunoprecipitation assay)',
                    'TRAP': 'MI:0045(experimental interaction detection)',
                    'Immunoprecipitaion': 'MI:1017(rna immunoprecipitation)',
                    'LacZ reporter assay': 'MI:0045(experimental interaction detection)',
                    'flow': 'MI:0045(experimental interaction detection)',
                    'EMSA': 'MI:0045(experimental interaction detection)',
                    'Communoprecipitaion': 'MI:1017(rna immunoprecipitation)',
                    'pSILAC': 'MI:0045(experimental interaction detection)',
                    'RTPCR': 'MI:1196(quantitative reverse transcription pcr)',
                    'proteomics analysis': 'MI:0045(experimental interaction detection)',
                    'immunoblot': 'MI:0045(experimental interaction detection)',
                    'ASO assay': 'MI:0045(experimental interaction detection)',
                    'semi-qRT-PCR': 'MI:1196(quantitative reverse transcription pcr)',
                    'mice xenograft': 'MI:0045(experimental interaction detection)',
                    'Chip': 'MI:0402(chromatin immunoprecipitation assay)',
                    'Flow cytometry': 'MI:0045(experimental interaction detection)',
                    'Immuohistochemistry': 'MI:0045(experimental interaction detection)',
                    'Chromatin immunoprecipitation': 'MI:0402(chromatin immunoprecipitation assay)',
                    'microarray': 'MI:0008(array technology)',
                    'Western blotting': 'MI:0113(western blot)',
                    'TaqMan miRNA assay/RT-PCR': 'MI:0045(experimental interaction detection)|MI:1196(quantitative reverse transcription pcr)',
                    'TaqMan miRNA assay': 'MI:0045(experimental interaction detection)',
                    'QRTPCRWestern blot': 'MI:1196(quantitative reverse transcription pcr)|MI:0113(western blot)',
                    'Gluc assay': 'MI:0045(experimental interaction detection)',
                    'Real time PCR': 'MI:0045(experimental interaction detection)',
                    "3'LIFE": 'MI:0045(experimental interaction detection)',
                    'Annexin V-FITC': 'MI:0045(experimental interaction detection)',
                    "5\\'RACE": 'MI:0045(experimental interaction detection)',
                    'Real time RT-PCR': 'MI:1196(quantitative reverse transcription pcr)',
                    'Luciferase assay/RT-PCR': 'MI:2285(miRNA interference luciferase reporter assay)|MI:1196(quantitative reverse transcription pcr)',
                    'Westren blot': 'MI:0113(western blot)',
                    '2DGE': 'MI:0045(experimental interaction detection)',
                    'Mass spectrometry': 'MI:0943(detection by mass spectrometry)',
                    'EGFP reporter assay': 'MI:0045(experimental interaction detection)',
                    ' Western blot': 'MI:0113(western blot)',
                    'AGO2 binding RNA immunoprecipitation qRT-PCR': 'MI:1196(quantitative reverse transcription pcr)',
                    'B-globin reporter assay': 'MI:0045(experimental interaction detection)',
                    'RISC-IP': 'MI:1017(rna immunoprecipitation)',
                    'Western Blotting': 'MI:0113(western blot)',
                    'Immunoprecipitation': 'MI:1017(rna immunoprecipitation)',
                    'GFP reporter': 'MI:0045(experimental interaction detection)',
                    'pMIR-REPORT': 'MI:0045(experimental interaction detection)',
                    'LacZ assay': 'MI:0045(experimental interaction detection)',
                    "5'RACE": 'MI:0045(experimental interaction detection)',
                    'Western blog': 'MI:0113(western blot)',
                    'Western blo': 'MI:0113(western blot)',
                    'western blot': 'MI:0113(western blot)',
                    'Reverse-phase protein array': 'MI:0008(array technology)',
                    'Western Blot': 'MI:0113(western blot)',
                    'MTT assay': 'MI:0045(experimental interaction detection)',
                    'Immunofluorescence staining': 'MI:0045(experimental interaction detection)',
                    'Immunoblotting': 'MI:0045(experimental interaction detection)',
                    'SILAC (Stable Isotope Labeling of Amino acids in Culture)': 'MI:0045(experimental interaction detection)',
                    'Western blot, luciferase assay': 'MI:0113(western blot)|MI:2285(miRNA interference luciferase reporter assay)',
                    'DNA methylation analysis': 'MI:1189(methylation interference assay)',
                    'Wetsern blot': 'MI:0113(western blot)',
                    'Immunohistochemistry analysis': 'MI:1198(immunohistochemistry)',
                    'ChIP-PCR': 'MI:0402(chromatin immunoprecipitation assay)',
                    'luciferase reporter assays': 'MI:2285(miRNA interference luciferase reporter assay)',
                    'PCR array': 'MI:0008(array technology)',
                    'Western': 'MI:0113(western blot)',
                    'immunostaining': 'MI:0422(immunostaining)',
                    'Caspase-GloÂ® 3/7 assay': 'MI:0045(experimental interaction detection)',
                    'Cell proliferation assay': 'MI:0045(experimental interaction detection)',
                    'safranin o staining/GAGs contents assay': 'MI:0045(experimental interaction detection)',
                    'wound healing assays': 'MI:0045(experimental interaction detection)',
                    'transwell insert': 'MI:0045(experimental interaction detection)',
                    'anoikis assay': 'MI:0045(experimental interaction detection)',
                    'Gluc reporter assay': 'MI:0045(experimental interaction detection)',
                    'GUS reporter assay': 'MI:0045(experimental interaction detection)',
                    'Zymography': 'MI:0512(zymography)',
                    'Motility assay': 'MI:0045(experimental interaction detection)',
                    'CAM assay': 'MI:0045(experimental interaction detection)',
                    'Colony formation assay': 'MI:0045(experimental interaction detection)',
                    'Alizarin red S staining': 'MI:0045(experimental interaction detection)',
                    'mRNA decay': 'MI:0045(experimental interaction detection)',
                    'Cell proliferation': 'MI:0045(experimental interaction detection)',
                    'apoptosis': 'MI:0045(experimental interaction detection)',
                    'cell cycle assays': 'MI:0045(experimental interaction detection)',
                    'colony formation': 'MI:0045(experimental interaction detection)',
                    'Immunoflourescence': 'MI:0045(experimental interaction detection)',
                    'Micorarray': 'MI:0008(array technology)',
                    'Westren Blot': 'MI:0113(western blot)',
                    'Luciferase reporter assay/Western blot': 'MI:2285(miRNA interference luciferase reporter assay)|MI:0113(western blot)',
                    'Immunohistochemical (IHC) staining': 'MI:1198(immunohistochemistry)',
                    'Luciferase reporter assay/qRT-PCR': 'MI:2285(miRNA interference luciferase reporter assay)|MI:1196(quantitative reverse transcription pcr)',
                    '5"RACE': 'MI:0045(experimental interaction detection)',
                    'Immunofluorescence analysis': 'MI:0045(experimental interaction detection)',
                    'luciferase reporter assay': 'MI:2285(miRNA interference luciferase reporter assay)',
                    'Wstern blot': 'MI:0113(western blot)',
                    'Coimmunoprecipitation': 'MI:1017(rna immunoprecipitation)',
                    'Immunofluorescence microscopy': 'MI:0045(experimental interaction detection)',
                    '/Western blot': 'MI:0113(western blot)',
                    'Luciferase reporter assay/QRTPCR': 'MI:2285(miRNA interference luciferase reporter assay)|MI:1196(quantitative reverse transcription pcr)',
                    'MTT': 'MI:0045(experimental interaction detection)',
                    'immunofluorescence assays': 'MI:0045(experimental interaction detection)',
                    'qRT_PCR': 'MI:1196(quantitative reverse transcription pcr)',
                    '2-D Gel Electrophoresis (2DGE)': 'MI:0982(electrophoretic mobility-based method)',
                    'RISC analysis': 'MI:0045(experimental interaction detection)',
                    'silico analysis': 'MI:0045(experimental interaction detection)',
                    'Microarray/In situ hybridization': 'MI:0008(array technology)',
                    'Western blot ': 'MI:0113(western blot)',
                    'Genotyping': 'MI:0045(experimental interaction detection)',
                    'Weastern blot': 'MI:0113(western blot)',
                    'YFP expression': 'MI:0045(experimental interaction detection)',
                    'To test if miR-141 directly targets the PR transcript, we analyzed four predicted miR-141-binding sites (Figure 4c)': 'MI:0045(experimental interaction detection)',
                    ' three within the 3â€˛ untranslated region (UTR) as identified through Targetscan (http:': 'MI:0045(experimental interaction detection)',
                    'www.targetscan.org/) and one in the la': 'MI:0045(experimental interaction detection)',
                    'qRT-PCR/Luciferase reporter assay': 'MI:1196(quantitative reverse transcription pcr)',
                    'Luciferase reporter assay and western blot': 'MI:2285(miRNA interference luciferase reporter assay)|MI:0113(western blot)',
                    'TOPď¬‚ash/FOPď¬‚ash reporter assay': 'MI:0045(experimental interaction detection)',
                    'dual-luciferase reporter assay': 'MI:0045(experimental interaction detection)',
                    'RNA immunoprecipitation assay (RIP)': 'MI:1017(rna immunoprecipitation)',
                    'Chromogenic in situ hybridization': 'MI:0045(experimental interaction detection)',
                    'Luciferase reporter assa': 'MI:2285(miRNA interference luciferase reporter assay)',
                    'Immunoprecipitaionă„ĄLuciferase reporter assay': '|MI:2285(miRNA interference luciferase reporter assay)',
                    'Luciferase reporter assay/Microarray': 'MI:2285(miRNA interference luciferase reporter assay)|MI:0008(array technology)',
                    'q-PCR': 'MI:1196(quantitative reverse transcription pcr)',
                    'AGO2 Immunoprecipitation': 'MI:1017(rna immunoprecipitation)',
                    'Cell proliferation assays': 'MI:0045(experimental interaction detection)',
                    'LC-MS/MS': 'MI:0943(detection by mass spectrometry)',
                    'Chromatin Immunoprecipitation': 'MI:0402(chromatin immunoprecipitation assay)',
                    'Co-immunoprecipitation': 'MI:1017(rna immunoprecipitation)',
                    'IlluminaExpressionArrays': 'MI:0008(array technology)',
                    'Protein Immunoblot Analyses': 'MI:0045(experimental interaction detection)',
                    'miR PCR array system': 'MI:0008(array technology)',
                    'mtt': 'MI:0045(experimental interaction detection)',
                    'RNA immunopercipitation': 'MI:1017(rna immunoprecipitation)',
                    'TOP/FOP luciferase assay': 'MI:2285(miRNA interference luciferase reporter assay)',
                    'miRNA-masking antisense ODN (miR-Mask) assay': 'MI:0045(experimental interaction detection)',
                    'enzyme-linked immunosorbent assay': 'MI:0045(experimental interaction detection)',
                    'Ago2-IP/IgG-IP': 'MI:1017(rna immunoprecipitation)',
                    'EGFR reporter assay': 'MI:0045(experimental interaction detection)',
                    'immunoblot analysis': 'MI:0045(experimental interaction detection)',
                    'Immunohistochemical analysis': 'MI:1198(immunohistochemistry)',
                    'CC tissues and cells (C33A, HeLa, CaSki, SiHa, and ME-180)': 'MI:0045(experimental interaction detection)',
                    'Immuno-precipitation': 'MI:1017(rna immunoprecipitation)',
                    'Luciferase reporter assayMTT': 'MI:2285(miRNA interference luciferase reporter assay)',
                    'Immunostaining': 'MI:0422(immunostaining)',
                    'immunosorbent': 'MI:0411(enzyme linked immunosorbent assay)',
                    'Immunofluorescent Assay': 'MI:0045(experimental interaction detection)',
                    'YFP reporter assay': 'MI:0045(experimental interaction detection)',
                    'CLIP-seq': 'MI:2191(clip)',
                    'RNAi': 'MI:0045(experimental interaction detection)',
                    '': 'MI:0045(experimental interaction detection)',
                    'microarray assay': 'MI:0008(array technology)',
                    'luciferase reporter': 'MI:2285(miRNA interference luciferase reporter assay)',
                    'miRNA Antagonists - qRT - PCR - Western Luciferase assay - mutant MRE': 'MI:0045(experimental interaction detection)',
                    'activity assay': 'MI:0045(experimental interaction detection)',
                    'real - time RT - PCR.': 'MI:1196(quantitative reverse transcription pcr)',
                    'Overexpression of miR - 155(pSILAC) down 50 - 25 %': 'MI:0045(experimental interaction detection)',
                    'ectopic expression - In silico prediction - qRT-PCR -Western - DNA methylation assays Luciferase assay - scrambled miRNA': 'MI:2285(miRNA interference luciferase reporter assay)',
                    'qRT-PCR Luciferase assay - mutant MRE': 'MI:2285(miRNA interference luciferase reporter assay)',
                    'In silico prediction Luciferase assay': 'MI:2285(miRNA interference luciferase reporter assay)',
                    'In silico prediction - miRNA Antagonist Luciferase assay - mutant MRE': 'MI:2285(miRNA interference luciferase reporter assay)',
                    'ectopic expression - miRNA Antagonist - In silico prediction Luciferase assay - mutant MRE': 'MI:2285(miRNA interference luciferase reporter assay)',
                    'In silico prediction - Dicer depletion Luciferase assay': 'MI:2285(miRNA interference luciferase reporter assay)',
                    ' Luciferase assay - mutant MRE': 'MI:2285(miRNA interference luciferase reporter assay)',
                    ' Luciferase assay': 'MI:2285(miRNA interference luciferase reporter assay)',
                    'RT-PCR': 'MI:1196(quantitative reverse transcription pcr)',
                    'Overexpression of let-7b (pSILAC) down 50-25%': 'MI:0045(experimental interaction detection)',
                    'Knock-down of let-7b (pSILAC) up 25-50%': 'MI:0045(experimental interaction detection)',
                    'RT-qPCR': 'MI:1196(quantitative reverse transcription pcr)',
                    'Overexpression of miR-1 (pSILAC) down 75-50%': 'MI:0045(experimental interaction detection)',
                    'Overexpression of let-7b (pSILAC) down <25%': 'MI:0045(experimental interaction detection)',
                    'miRNA Antagonists - qRT-PCR - Western Luciferase assay - mutant MRE': 'MI:2285(miRNA interference luciferase reporter assay)',
                    'real-time RT-PCR.': 'MI:1196(quantitative reverse transcription pcr)',
                    'Overexpression of miR-155 (pSILAC) down 50-25%': 'MI:0045(experimental interaction detection)',
                    'Real-time RT-PCR assay': 'MI:1196(quantitative reverse transcription pcr)',
                    'In silico prediction Luciferase assay - mutant MRE': 'MI:2285(miRNA interference luciferase reporter assay)',
                    'real time RT-PCR': 'MI:1196(quantitative reverse transcription pcr)',

                }

                detlist = []
                for method in columns[4].split('//'):
                    for real_method in method.split(';'):
                        if real_method not in detmap and 'MI:0045(experimental interaction detection)' not in detlist:
                            detlist.append('MI:0045(experimental interaction detection)')
                        else:
                            if detmap[real_method] not in detlist:
                                detlist.append(detmap[real_method])

                # Inserting edges
                edge_dict = {
                    'publication_ids': 'pubmed:' + columns[3],
                    'layer': '5',
                    'source_db': 'miRSponge',
                    'interaction_identifiers': None,
                    'confidence_scores': None,
                    'interaction_detection_method': '|'.join(detlist),
                    'interaction_types': interaction_types,
                    'first_author': None
                }

                db_api.insert_edge(source_dict, target_dict, edge_dict)


# Saving the to a DB_TYPE.db files
    db_api.save_db_to_file(DB_DESTINATION)


if __name__ == '__main__':
    print("Parsing database...")
    main(logger = None)
    print("Parsing database is completed. SQLite database is saved to: " + DB_DESTINATION)


