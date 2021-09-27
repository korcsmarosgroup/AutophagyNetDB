
# Defining constants
# TODO: add all species and experiments
# TODO mapping
RNA_SEQ_LIST = ['files/Caenorhabditis_elegans_RNA-Seq_read_counts_TPM_FPKM_GSE16552.tsv']
AFF_LIST = ['files/Caenorhabditis_elegans_probesets_GSE19972_A-AFFY-60_gcRMA.tsv']
exp_dict = {}

for file in RNA_SEQ_LIST:
    # Opening the library file
    with open(file) as lib_file:
        lib_file.readline()
        for line in lib_file:
            line = line.strip().split('\t')
            if line[3] not in exp_dict:
                exp_dict[line[3]] = {'tissue_name': [],
                                     'stage': [],
                                     'FPKM': [],
                                     'presence': []}
            else:
                if line[5].replace('"', '') not in exp_dict[line[3]]['tissue_name'] \
                    and line[7].replace('"','') not in exp_dict[line[3]]['stage'] \
                        and line[12] not in exp_dict[line[3]]['FPKM'] \
                        and line[13] not in exp_dict[line[3]]['presence']:
                    exp_dict[line[3]]['tissue_name'].append(line[5].replace('"',''))
                    exp_dict[line[3]]['stage'].append(line[7].replace('"', ''))
                    exp_dict[line[3]]['FPKM'].append(line[12])
                    exp_dict[line[3]]['presence'].append(line[13])

chip_dict = {}
for file in AFF_LIST:
    with open(file) as chipfile:
        chipfile.readline()
        for line in chipfile:
            line = line.strip().split('\t')
            if line[3] not in chip_dict:
                chip_dict[line[3]] = {'tissue_name': [],
                                      'stage': [],
                                      'MAS5': [],
                                      'presence': []}
            else:
                if line[5].replace('"', '') not in chip_dict[line[3]]['tissue_name'] \
                    and line[7].replace('"','') not in chip_dict[line[3]]['stage'] \
                        and line[10] not in chip_dict[line[3]]['MAS5'] \
                        and line[11] not in chip_dict[line[3]]['presence']:
                    chip_dict[line[3]]['tissue_name'].append(line[5].replace('"',''))
                    chip_dict[line[3]]['stage'].append(line[7].replace('"', ''))
                    chip_dict[line[3]]['MAS5'].append(line[10])
                    chip_dict[line[3]]['presence'].append(line[11])

print(exp_dict['WBGene00173312'])
print(chip_dict['WBGene00000027'])

