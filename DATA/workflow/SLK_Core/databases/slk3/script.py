"""
This script parses the signalink tsv files located in it's directory and generates a Psi-Mi SQLite database files from them.
The .tsv files must be in the following format:
    1) The first line is the header line
    2) Columns separated with tab: species, apathway, agenesymbol, auniprot, apathway, atopology, bgenesymbol, buniprot, bpathway, btopology, is_direct, interaction_type, effect, reference
    :argument: DB_TYPE: name of the database
    :argument: DESTINATION: saving destination
    :argument: TSV_LIST: list of the data files
    :argument: ORGANISM_NAME_TO_MITAB_ID_MAP: species name to taxid
    :argument: IS_DIRECT_MAP: dictionary of directness to MI ids
    :argument: EFFECT_MAP: dictionary of effect to MI ids
    :argument: MOLECULAR_MAP: dictionary of molecular background to MI ids 
"""

# Imports
import csv
from ARNlib.SQLiteDBApi.sqlite_db_api import PsimiSQL

# Defining constants
SQL_SEED = '../../../../../ARNlib/SQLiteDBApi/network-db-seed.sql'
DB_TYPE = 'slk3'
DESTINATION = '../../../layer0/output/slk3'
TSV_LIST = ['files/SLK3_update_cele_kitti_betti.txt',
            'files/SLK3_update_dmel_kitti_viktor.txt',
            'files/SLK3_update_hsap_kitti_d.txt']
ORGANISM_NAME_TO_MITAB_ID_MAP = {
    "H. sapiens": "taxid:9606",
    "D. melanogaster": "taxid:7227",
    "C. elegans": "taxid:6239"
}

IS_DIRECT_MAP = {
    "direct": "MI:0407(direct interaction)",
    "dierct" : "MI:0407(direct interaction)",
    "indirect": "indirect",
    "direct/indirect": "unknown"
}

EFFECT_MAP = {
    'Unknown': 'MI:0190(interaction type)',
    'down-regulates': 'MI:2240(down-regulates)',
    "down-regulates activity": 'MI:2241(down-regulates activity)',
    "down-regulates quantity": 'MI:2242(down-regulates quantity)',
    "down-regulates quantity by destabilization": 'MI:2244(down-regulates quantity by destabilization)',
    "down-regulates quantity by repression": 'MI:2243(down-regulates quantity by repression)',
    'unknown': 'MI:0190(interaction type)',
    'up-regulates': 'MI:2235(up-regulates)',
    "up-regulates activity": 'MI:2236(up-regulates activity)',
    "up-regulates quantity by expression": 'MI:2238(up-regulates quantity by expression)',
    "up-regulates quantity by stabilization": 'MI:2239(up-regulates quantity by stabilization)',
    "inhibition" : 'MI:0623(inhibition)',
    "inhibiton" : 'MI:0623(inhibition)',
    "stimulation" : 'MI:0624(stimulation)'
}

MOLECULAR_MAP = {
    'binding': 'MI:0462(bind)',
    'transcriptional regulation' : 'MI:2247(transcriptional regulation)',
    'phosphorylation' : 'MI:0217(phosphorylation reaction)',
    'phosphorilation' : 'MI:0217(phosphorylation reaction)',
    'phosphoryltion' : 'MI:0217(phosphorylation reaction)',
    'phosphotilation' : 'MI:0217(phosphorylation reaction)',
    '' : 'MI:0190(interaction type)',
    'ubiquitination' : 'MI:0220(ubiquitination reaction)',
    'ubiquitinisation' : 'MI:0220(ubiquitination reaction)',
    'relocalization' : 'MI:2256(relocalization reaction)',
    'dephosphorylation' : 'MI:0203(dephosphorylation reaction)',
    'dephosphotilation' : 'MI:0203(dephosphorylation reaction)',
    'dephosphorilation' : 'MI:0203(dephosphorylation reaction)',
    'cleavage' : 'MI:0194(cleavage reaction)',
    'deubiquitination' : 'MI:0204(deubiquitination reaction)',
    'deubiquitinisation' : 'MI:0204(deubiquitination reaction)',
    'guanine nucleotide exchange factor' : 'MI:2252(guanine nucleotide exchange factor)',
    'sumoylation' : 'MI:0566(sumoylation reaction)',
    'dimethylation' : 'MI:0213(methylation reaction)',                #MI id needed
    'acetylation' : 'MI:0192(acetylation reaction)',
    'deacetylation' : 'MI:0197(deacetylation reaction)',
    'deacetyltion' : 'MI:0197(deacetylation reaction)',
    'autophosphorylation' : 'MI:0217(phosphorylation reaction)',        #MI id needed
    'autophosphorilation': 'MI:0217(phosphorylation reaction)',
    'methylation' : 'MI:0213(methylation reaction)',
    'monoubuquitinisation' : 'MI:0220(ubiquitination reaction)',      #MI id needed
    'polyubiquitinisation' : 'MI:0220(ubiquitination reaction)',      #MI id needed
    'autophosphorilation of mapk14' : 'MI:0217(phosphorylation reaction)'     #MI id needed

}

accepted_pathways = {
    "TCR": 'T-cell receptor',
    "BCR": 'B-cell receptor',
    "TLR": 'Toll-like receptor',
    "IIP": 'Innate immune pathways',
    "JAK/STAT": 'JAK/STAT',
    "NHR": 'Nuclear hormone receptor',
    "RTK": 'Receptor tyrosine kinase',
    "Rho/Cytoskeleton": 'Rho pathway',
    "TGF": 'TGF',
    "Notch": 'Notch',
    "GPCR": 'G-protein coupled receptor',
    "WNT/Wingless": 'WNT/Wingless',
    "WNT": 'WNT/Wingless',
    "HIPPO": 'Hippo',
    "HH": 'Hedgehog',
    "TNF/Apoptosis": 'TNF pathway',
}


def get_mitab_pathways_list(pathways):
    pathways_list = pathways.split(',')
    new_pathway_list = []
    for p in pathways_list:
        new_pathway = p.split("(")[0]
        if new_pathway == "Hedgehog":
            new_pathway_short = "HH"
            new_pathway_list.append(new_pathway_short)
        elif new_pathway == "NOTCH":
            new_pathway_short = "Notch"
            new_pathway_list.append(new_pathway_short)
        else:
            new_pathway_list.append(new_pathway)

    return new_pathway_list


def main(logger):
    db_api = PsimiSQL(SQL_SEED)

    # looping through SLK3 files
    for SLK_3_FILE_LOCATION in TSV_LIST:

        #opening each slk files and looping through it
        SLK_3_FILE = csv.reader(open(SLK_3_FILE_LOCATION, encoding="ISO-8859-1"), delimiter = '\t', quotechar = '"')
        next(SLK_3_FILE) # Skipping the header

        for line in SLK_3_FILE:

            pathways_a = get_mitab_pathways_list(line[4])
            new_pathways_a = []
            for p in pathways_a:
                pathway_a = p
                if " " in p:
                    pathway_a = p.replace(" ", "")
                elif '"' in p:
                    pathway_a = p.replace('"', "")
                new_p = accepted_pathways[pathway_a]
                new_pathways_a.append(new_p)

            new_node_a = line[3]
            if " " in line[3]:
                new_node_a = line[3].replace(" ", "")

            if line[5] == "Meditor":
                line[5] = "Mediator"

            topologies_a = set(map(lambda x: x.strip(), line[5].split(",")))

            source_dict = {
                "name" : "Uniprot:" + new_node_a,
                "alt_accession" : "gene symbol:" + line[2],
                "tax_id" : ORGANISM_NAME_TO_MITAB_ID_MAP[line[0]],
                "aliases" : '-',
                "pathways" : "|".join(new_pathways_a),
                "topology" : "|".join(topologies_a)
            }

            db_api.insert_node(source_dict)

            pathways_b = get_mitab_pathways_list(line[9])
            new_pathways_b = []
            for p in pathways_b:
                pathway_b = p
                if " " in p:
                    pathway_b = p.replace(" ", "")
                elif '"' in p:
                    pathway_b = p.replace('"', "")
                new_p = accepted_pathways[pathway_b]
                new_pathways_b.append(new_p)

            new_node_b = line[8]
            if " " in line[8]:
                new_node_b = line[8].replace(" ", "")

            topologies_b = set(map(lambda x: x.strip(), line[10].split(",")))

            target_dict = {
                "name" : "Uniprot:" + new_node_b,
                "alt_accession" : "gene symbol:" + line[7],
                "tax_id" : ORGANISM_NAME_TO_MITAB_ID_MAP[line[0]],
                "aliases" : '-',
                "pathways" : "|".join(new_pathways_b),
                "topology" : "|".join(topologies_b)
            }

            db_api.insert_node(target_dict)

            effect = EFFECT_MAP[line[14].lower()]

            molecular_background = MOLECULAR_MAP[line[13].lower()]

            inttype_final = effect + '|' + molecular_background

            is_direct = IS_DIRECT_MAP[line[12].lower()]
            if "MI:0407(direct interaction)" in is_direct:
                is_direct = "true"
            else:
                is_direct = "false"

            interaction_types = "%s|is_directed:%s|is_direct:%s" % (inttype_final, "true", is_direct)

            edge_dict = {
                'interaction_detection_method' : None,
                'first_author' : None,
                'publication_ids' : 'pubmed:' + line[15],
                'interaction_types' : interaction_types,
                'source_db' : 'SLKv3.0',
                'interaction_identifiers' : None,
                'confidence_scores' : None,
                'layer' : "3"
            }

            db_api.insert_edge(source_dict, target_dict, edge_dict)

    db_api.save_db_to_file(DESTINATION)


if __name__ == '__main__':
    main()
