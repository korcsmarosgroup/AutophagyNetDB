"""
 Tissue expression data from: https://bgee.org/?page=download&action=proc_values#id5
"""

# Imports
import logging, os, sqlite3, json, io

# Getting files names
RNA_SEQ_LIST = []
AFF_LIST = []
for directory in os.listdir('files'):
    if directory.endswith('.zip'):
        pass
    else:
        if 'RNA' in directory:
            for file in os.listdir('files/Homo_sapiens_RNA-Seq_read_counts_TPM_FPKM/'):
                RNA_SEQ_LIST.append('files/Homo_sapiens_RNA-Seq_read_counts_TPM_FPKM/' + file)
        elif 'Affymetrix' in directory:
            for file in os.listdir('files/Homo_sapiens_Affymetrix_probesets/'):
                AFF_LIST.append('files/Homo_sapiens_Affymetrix_probesets/' + file)


gene_file = 'Amanda_ATG_tissue_infile.txt'
outfile = 'Amanda_ATG_tissue.txt'

linestowrite = []

with open(gene_file) as genes, open(outfile, 'w') as outfile:
    for line in genes:
        line = line.strip().split('\t')
        ensid = line[1]
        uniprot = line[0]
        for file in RNA_SEQ_LIST:
            with open(file) as lib_file:
                lib_file.readline()
                for line in lib_file:
                    line = line.strip().split('\t')
                    if line[3] == ensid:
                        writeline = uniprot + '\t' + line[5] + '\n'
                        if writeline not in linestowrite:
                            linestowrite.append(writeline)
        for file in AFF_LIST:
            with open(file) as chipfile:
                chipfile.readline()
                for line in chipfile:
                    line = line.strip().split('\t')
                    if line[3] == ensid:
                        writeline = uniprot + '\t' + line[5] + '\n'
                        if writeline not in linestowrite:
                            linestowrite.append(writeline)

    for elem in linestowrite:
        outfile.write(elem)
