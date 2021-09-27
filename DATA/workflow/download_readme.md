
Layer0
ACSN (v2)
	Latest date of download: 2019. 03. 07.
	
    PATHWAY_FILE: PPI interactions in .gmt format https://acsn.curie.fr/files/acsn_master_curated.gmt
    ALL_EDGE_FILE_LOCATION: Correspondance between ACSN entities and HUGO names in .gmt format: https://acsn.curie.fr/files/acsn_names.gmt
    CURATED_PROTEIN_LIST_FILE_LOCATION: PPI interaction in .sif format https://acsn.curie.fr/files/acsn_ppi.sif	

InnateDB (v5.4)
	Latest date of download: 2019. 03. 07.

    DATA_FILE: PPI files in mitab format: http://www.innatedb.com/download/interactions/innatedb_ppi.mitab.gz
    XGMML_LIST: list of .xgmml files for each pathway
	
	Statistics
		19800 human interactions

Reactome (v3.3)
	Latest date of download: 2019. 03. 07.

    DATA_FILE: human PPI pairs in psi-mitab format: http://www.reactome.org/download/current/homo_sapiens.mitab.interactions.txt.gz
    PATHWAY_FILE_LOCATION: files with 2 columns, source pathway id and signalink id
    PUBMED_INTERACTION_FILE:human PPI pairs in tab deliminated format: http://www.reactome.org/download/current/homo_sapiens.interactions.txt.gz
    UNI_TO_PATHWAY: files with uniprot ids and their pathways from: https://reactome.org/download/current/UniProt2Reactome_All_Levels.txt
	
	Statistics	
		Human: 
			10792 proteins (not including isoforms)
			12278 reactions
		Drosophila:
			9888 proteins
			3990 reactions
		Danio
			14508 proteins
			6425 reactions
		C. elegans
			5211 proteins
			2727 reactions

Signor (v2.0)
	Latest date of download: 2017. 03. 24.

    CSV_LIST: list of .csv files of each signalink pathway
    FILENAME_TO_PATHWAY_MAP: dictionary from files name to SLK pathway
	
	Statistics	
		Human: 

Layer1
PSP
	Latest date of download: 2019. 03. 07.

    DATA_FILE: from http://www.cell.com/trends/cell-biology/fulltext/S0962-8924(09)00271-2 article: potential scaffold proteins
	(usues BioMyn data)
Layer2
PhosphoSitePlus
	Latest date of download: 2019. 03. 07.

    DATA_FILE: http://www.phosphosite.org/staticDownloads.action: Kinase_Substrate_Dataset.gz files

PTMCode2 
	Latest date of download: 2019. 03. 07.

    DATA_FILE: https://ptmcode.embl.de/data/PTMcode2_associations_between_proteins.txt.gz 
ELM pred
	Latest date of download: 2019. 03. 07.

    ELMS_FILE: all ELM classes of the four used species in a .tsv files: http://elm.eu.org/elms/elms_index.tsv
    INT_DOMAINS_FILE: files containing ELM names and their interacting domain PFAM ids in a .tsv files: http://elm.eu.org/interactiondomains
    PROT_LIST: list of files for each species used, containing their whole proteomes from UniProt in .fa files

Layer3
biogrid (v3.5)
	Latest date of download: 2019. 03. 07.

    BIOGRID-SYSTEM-Positive_Genetic-3.4.145.mitab.txt
	
	Statistics
		Human: 
			491431 physical interactions (365,035 non redundant interactions)
		C. elegans:
			6360 physical interactions (5,808 non redundant interactions)
		Drosophiila:
			62165 physical interactions (52,366 non redundant interactions)
		Danio:
			223 physical interactions (196 non redundant interactions)
			
ComPPI (v2.1.1)
	Latest date of download: 2019. 03. 11.

    http://comppi.linkgroup.hu/downloads  integrated PPIs comppi--compartments--tax_celegans_loc_all.txt

HPRD (release 9)
	Latest date of download: 2019. 03. 07.

    HPRD_Release9_062910/BINARY_PROTEIN_PROTEIN_INTERACTIONS.txt
IntAct (v4.2.12)
	Latest date of download: 2019. 03. 11.

    DATA_FILE: ftp://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact-micluster.txt
OmniPath (v0.7.111)
	Latest date of download: 2019. 03. 07.

    DATA_FILE: http://omnipathdb.org/interactions/?fields=sources&fields=references

Layer5
mir2disease
	Latest date of download: 2017. 05. 07.

    DATA_FILE: miRNA-target data files http://watson.compbio.iupui.edu:8080/miR2Disease/download/miRtar.txt
mirdeathdb
	Latest date of download: 2017. 05. 04.

    DATA_FILE: http://www.rna-world.org/mirdeathdb/data/miRDeathDB_all_data.txt
mirecords (v4)
	Latest date of download: 2019. 03. 07.

    DATA_FILE: validated target data files: http://c1.accurascience.com/miRecords/download_data.php?v=4
starBase (v3.0)
	Latest date of download: 2019. 03. 11.

    DATA_FILE: human predicted miRNA-mRNA: http://starbase.sysu.edu.cn/moduleDownload.php?source=agoClipRNA&type=txt&value=hg19;mRNA;all;1;0;0;1;None;PDCD4
			   human degradome miRNA-mRNA: http://starbase.sysu.edu.cn/moduleDownload.php?source=degradomeRNA&type=txt&value=hg19;mRNA;hsa-let-7a-2-3p;1;all
			   worm degradome miRNA-mRNA: http://starbase.sysu.edu.cn/moduleDownload.php?source=degradomeRNA&type=txt&value=ce11;mRNA;cel-let-7-3p;1;all
			   human validated miRNA-RNA: http://starbase.sysu.edu.cn/moduleDownload.php?source=rnaRNA&type=txt&value=hg19;hsa-let-7a-5p;1;1
mirTarbase (v7.0)
	Latest date of download: 2019. 03. 07.

    DATA_FILE_LIST: data files of each species
                    human: http://mirtarbase.mbc.nctu.edu.tw/cache/download/6.1/hsa_MTI.xlsx
                    celegans: http://mirtarbase.mbc.nctu.edu.tw/cache/download/6.1/cel_MTI.xls
                    daino: http://mirtarbase.mbc.nctu.edu.tw/cache/download/6.1/dre_MTI.xls
                    drosi: http://mirtarbase.mbc.nctu.edu.tw/cache/download/6.1/dme_MTI.xls

Layer7
lncRinter
	Latest date of download: 2019. 03. 11.

    DATA_FILE: //bioinfo.life.hust.edu.cn/lncRInter/browse?species=&class=&level=RNA-RNA
mirsponge (v1.0)
	Latest date of download: 2019. 03. 11.

    DATA_FILE: http://www.bio-bigdata.com/miRSponge/apps/download1.jsp?fileName=Experimentally_validated_miRNA_targets_in_miRSponge.txt
NPInter(v3.0)
	Latest date of download: 2019. 03. 07.

    DATA_FILE: data files (linux version) http://www.bioinfo.org/NPInter/datadownload/interaction_NPInter[v3.0].txt.tar.gz
starbase (v3.0)
	Latest date of download: 2019. 03. 11.

    DATA_FILE_LIST: data files for each species
					human predicted miRNA-lncRNA: http://starbase.sysu.edu.cn/moduleDownload.php?source=agoClipRNA&type=txt&value=hg19;lncRNA;all;1;0;0;1;None;MALAT1
					human degradome miRNA-ncRNA: http://starbase.sysu.edu.cn/moduleDownload.php?source=degradomeRNA&type=txt&value=hg19;ncRNA;hsa-let-7a-2-3p;1;all
					worm degradome miRNA-ncRNA: http://starbase.sysu.edu.cn/moduleDownload.php?source=degradomeRNA&type=txt&value=ce11;ncRNA;cel-let-7-3p;1;all
					human validated lncRNA-RNA: http://starbase.sysu.edu.cn/moduleDownload.php?source=rnaRNA&type=txt&value=hg19;H19;1;1



