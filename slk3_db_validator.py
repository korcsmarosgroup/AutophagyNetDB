#!/usr/bin/env python3

import sqlite3
import re
import argparse
import sys

accepted_id_types = "|".join([
    # proteins
    'Uniprot',
    'RefSeq',
    'BioGrid',
    'IntAct',
    'MINT',
    'Ensembl',
    'GeneID',
    'GeneCards',
    'HGNC',
    'Reactome',
    'SignaLink',
    'SIGNOR',
    'ELM',
    'Ensemb',
    'FlyBase',
    'WormBase',

    # RNAs (we can also have Ensembl and HGNC IDs for RNAs, but they are already in the list)
    'RNACentral',
    'NONCODE',
    'miRBase',
    'LNCipedia',
])

accepted_pathways = "|".join([
    'T-cell receptor',
    'B-cell receptor',
    'Toll-like receptor',
    'Innate immune pathways',
    'JAK/STAT',
    'Nuclear hormone receptor',
    'Receptor tyrosine kinase',
    'Rho pathway',
    'TGF',
    'Notch',
    'G-protein coupled receptor',
    'WNT/Wingless',
    'Hippo',
    'Hedgehog',
    'TNF pathway',
	'RTK',
	'TCR',
	'NHR'
])

accepted_topology = "|".join([
    'Co-factor',
    'Scaffold',
    'Mediator',
    'Transcription factor',
    'Receptor',
    'Endocytosis related',
    'Ligand',
])

accepted_source_databases = "|".join([
    "ACSN", "InnateDB", "Reactome", "Signor", "SLKv2\\.0", "SLKv3\\.0", "SLKv2\\.1", "TCRcuration", "SignaFish", "SignaLink3", "SLK3_core",   # layer 0
    "PSP",                                                                                            # layer 1
    "PhosphoSite", "PTMCode2",                                                                        # layer 2
    "ComPPI", "HPRD", "IntAct", "TheBiogrid", "OmniPath", 'HumanAutophagyDB',                                           # layer 3
    "miRDeathDB", "miRecords", "miR2Disease", "TarBase", "StarBase",                                  # layer 5
    "PSSMprediction", "TFlink",                                                                       # layer 6
    "NPInter", "lncRInter", "miRSponge", "StarBase",                                                   # layer 7
    "ADB", "Behrends", "manual curation", "Behrends predicted", "ARN1Core", "coremancur", "manual_curation",
    "chip_behrends", "updated_curation"
])

accepted_layers = "|".join(['0', '1', '2', '3', '5', '6', '7', '8'])

mi_pattern = "MI:\\d{4}(\\(\\w[\\(\\)/ \\w-]+\\w\\))?"
mi_pattern_or_interaction_type = "|".join([
    "(%s)" % mi_pattern,
    "(is_direct:(true|false))",
    "(is_directed:(true|false|directed))",
])
pubmed_pattern = "pubmed:\\d+"
name_pattern = "(%s):[/\\.\\w-]+" % accepted_id_types
confidence_score_pattern = "[\\./ \\w-]+:(-)?[0-9]+(\\.[0-9]+)?"
tissue_pattern = "tissue:UBERON:\\d+(\\(\\w[ \\w]+\\w\\))?"
minor_localization_pattern = "minorloc:GO:\\d+(\\(\\w[ \\w]+\\w\\))?"
major_localization_pattern = "majorloc:GO:\\d+(\\(\\w[ \\w]+\\w\\))?"
topology_or_tissue_or_minor_loc_or_major_loc_pattern = "|".join([
    "(%s)" % accepted_topology,
    "(%s)" % tissue_pattern,
    "(%s)" % minor_localization_pattern,
    "(%s)" % major_localization_pattern
])


def zero_or_more_elements_of(pattern):
    return '((%s)(\\|(%s))*)?' % (pattern, pattern)


node_patterns = {
    'name': name_pattern,                                     # like 'uniprot:P12732' or 'some_type:P1.2732-2'
    'tax_id': 'taxid:\\d{4}',                                 # like 'taxid:9606'
    'pathways': zero_or_more_elements_of(accepted_pathways),  # like 'Notch|HH'
    'topology': zero_or_more_elements_of(topology_or_tissue_or_minor_loc_or_major_loc_pattern),  # like 'Ligand|Receptor|UBERON:0002101|minorloc:GO:1234|minorloc:GO:5555(some name)|majorloc:GO:123478'
    # 'aliases': '',        # we don't use this column
    # 'alt_accession': '',  # we don't use this column
}

edge_patterns = {
    'interactor_a_node_name': name_pattern,                                        # like 'uniprot:P12732' or 'some_type:P1.2732-2'
    'interactor_b_node_name': name_pattern,                                        # like 'uniprot:P12732' or 'some_type:P1.2732-2'
    'interaction_detection_method': zero_or_more_elements_of(mi_pattern),          # like 'MI:1234|MI:2345(jdk skdf)'
    'publication_ids': zero_or_more_elements_of(pubmed_pattern),                   # like 'pubmed:123764|pubmed:8362'
    'interaction_types': zero_or_more_elements_of(mi_pattern_or_interaction_type), # like 'is_direct:true|is_directed:false|MI:1234|MI:2345(jdk skdf)'
    'source_db': zero_or_more_elements_of(accepted_source_databases),              # like 'InnateDB|Signor'
    'layer': '(%s)?' % accepted_layers,                                            # like '3' or empty string
    'confidence_scores': zero_or_more_elements_of(confidence_score_pattern),       # like 'intact-miscore:0.59257585|other score (other db):low'
    # 'interaction_identifiers': '\\d+', # we don't use this column
    # 'first_author': '\\d+',  # we don't use this column
}


def validate_table(cursor, table_name, patterns):
    cursor.execute("SELECT * FROM " + table_name)

    all_rows = 0
    all_row_errors = 0
    malformed_columns = set()

    print("\n")

    for row in cursor:
        all_rows += 1
        error = False
        for (column, pattern) in patterns.items():
            value = ""
            if row[column]:
                value = str(row[column])
            if not re.search("^%s$" % pattern, value):
                error = True
                malformed_columns.add(column)
                print(
                    "ERROR: %s (id=%s) column %s malformed. value: `%s`" % (table_name, str(row['id']), column, value))
        if error:
            all_row_errors += 1

    print("\n===== %s table: =====" % table_name)
    print("total number of rows processed: %d" % all_rows)
    print("number of rows with wrong cells: %d" % all_row_errors)
    print("malformed columns: " + ", ".join(malformed_columns))

    return all_row_errors == 0


def parse_args():
    help_text = \
        """
        === Validate SLK 3 DB files ===
        """

    parser = argparse.ArgumentParser(description=help_text)

    parser.add_argument("-i", "--input_file",
                        help="<path to the input SQLite DB file> [mandatory]",
                        type=str,
                        dest="input_file",
                        action="store",
                        required=True)

    parser.add_argument("-t", "--type",
                        help="<basic (source / mapper / merger output) or builder> [optional]",
                        type=str,
                        dest="type",
                        action="store",
                        default="basic",
                        choices = ['basic', 'builder'],
                        required=False)

    results = parser.parse_args()

    return results.input_file, results.type


def validate_db_file(db_path, db_type='basic'):
    print("opening DB file for validation: " + db_path)
    db = sqlite3.connect(db_path)
    db.row_factory = sqlite3.Row
    cursor = db.cursor()

    valid = validate_table(cursor, "node", node_patterns)

    if db_type == 'basic':
        valid = validate_table(cursor, 'edge', edge_patterns) and valid
    else:
        for layer_table_name in ['layer0', 'layer1', 'layer2', 'layer3', 'layer5', 'layer6', 'layer7']:
            valid = validate_table(cursor, layer_table_name, edge_patterns) and valid

    print("\nvalidation of %s finished\n\n" % db_path)

    return valid


def main():
    input_file, db_type = parse_args()
    valid = validate_db_file(input_file, db_type)
    if not valid:
        sys.exit(1)


if __name__ == '__main__':
    main()
