
import sqlite3, json, io, logging, functools
logger_output_location = 'SLK3.log'
mapper_location = '../../DATA/workflow/mapper.db'

json_mapper_file = '../../DATA/workflow/uniprot_id_mapping.json'

layers_list = ['layer0',
               'layer1',
               'layer2',
               'layer3',
               'layer5',
               'layer6',
               'layer7',
               'layer8']

SPEC_MAP = {
     "taxid:9606": "H. sapiens"
}
SPEC_SHORT ={
     "taxid:9606": "human"

}
SPEC_LONG = {
     "taxid:9606": "Homo sapiens"
}

# Initiating logger
logger = logging.getLogger()
handler = logging.FileHandler(logger_output_location)
logger.setLevel(logging.DEBUG)
handler.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)

# Connecting to mapper database
# Read database to tempfile
map_conn = sqlite3.connect(mapper_location)
tempfile = io.StringIO()
for line in map_conn.iterdump():
    tempfile.write('%s\n' % line)
map_conn.close()
tempfile.seek(0)

# Create a database in memory and import from tempfile
map_db = sqlite3.connect(":memory:")
with map_db:
    map_db.cursor().executescript(tempfile.read())
    # we create a few indexes on the in-memory db, to optimize the queries later
    map_db.cursor().execute("CREATE INDEX uniprot_id ON UNIPROT_AC(id);")
    map_db.cursor().execute("CREATE INDEX uniprot_uniprot ON UNIPROT_AC(uniprot_ac);")
    map_db.cursor().execute("CREATE INDEX mapp_uniprot ON MAPP(uniprot_ac);")




# Mapping function
@functools.lru_cache(maxsize=None)
def map_uniprot_to_gene(uniprot):
    """
    Maps uniprot ids to gene names
    :param uniprot: uniprot id you'd like to map
    :return: gene name of the uniprot id if it's in the searched dataset,
             the input uniprot id if it's not in the searched dataset
    """
    # Connecting to mapper database
    with map_db:
        c2 = map_db.cursor()
        # Getting gene names where uniprot id matches
        c2.execute(
            "SELECT MAPP.gene_name FROM UNIPROT_AC LEFT JOIN MAPP ON UNIPROT_AC.id=MAPP.uniprot_ac "
            "WHERE UNIPROT_AC.uniprot_ac='%s' GROUP BY UNIPROT_AC.uniprot_ac"
            % uniprot.upper()
        )
        foreign_rows = c2.fetchone()
        if foreign_rows is not None:
            if foreign_rows[0] is not None:
                return foreign_rows[0]

        return uniprot


def map_uniprot_to_protein(uniprot):
    """
        Maps uniprot ids to protein names
        :param uniprot: uniprot id you'd like to map
        :return: protein name of the uniprot id if it's in the searched dataset,
                 the input uniprot id if it's not in the searched dataset
        """
    # Connecting to mapper database
    with map_db:
        c2 = map_db.cursor()
        # Getting protein names where uniprot id matches
        c2.execute(
            "SELECT MAPP.prot_full_name FROM UNIPROT_AC LEFT JOIN MAPP ON UNIPROT_AC.id=MAPP.uniprot_ac "
            "WHERE UNIPROT_AC.uniprot_ac='%s' GROUP BY UNIPROT_AC.uniprot_ac"
            % uniprot.upper()
        )
        foreign_rows = c2.fetchone()
        if foreign_rows is not None:
            return foreign_rows[0]
        else:
            return uniprot


def map_uniprot_to_external(json_mapper_file):
    """
    Maps primary uniprot ids to different external source ids based on the mapping table
    Also retrievs the database which the id is from (id type)
    :param json_mapper_file: entire proteome from uniprot in xml converted to json by mk_mapper in workflow/uniprot_id_mapping folder
            ./mk_mapper input.gz output.json
    :param uniprot: uniprot id that we want to find external references for
    :return: external reference id and database
    """

    mapdict = {}
    with open(json_mapper_file) as mapfile:
        for line in mapfile:
            line = json.loads(line)
            extrefdb = line['to_id_type']
            extrefid = line['to_id'].upper()
            if line['from_id_type'] == 'uniprotac':
                uniprot = line['from_id'].upper()
                if extrefdb == 'uniprotac' or \
                        extrefdb == 'hgnc' or \
                        extrefdb == 'ensembl' or \
                        extrefdb == 'flybase' or \
                        extrefdb == 'wormbase' or \
                        extrefdb == 'zfin':
                    if uniprot not in mapdict:
                        mapdict[uniprot] = {extrefdb: [extrefid]}
                    else:
                        if extrefdb not in mapdict[uniprot]:
                            mapdict[uniprot][extrefdb] = [extrefid]
                        else:
                            mapdict[uniprot][extrefdb].append(extrefid)
    return mapdict


def get_node_data(builder_location):
    # Getting NODE data
    logging.debug('Started getting node data')
    nodes = []

    # building the uniprot->externalRef map
    uniprot_to_external_id_dict = map_uniprot_to_external(json_mapper_file=json_mapper_file)

    builder_conn = sqlite3.connect(builder_location)

    with builder_conn:
        builder_conn.row_factory = sqlite3.Row
        merger_cursor = builder_conn.cursor()
        counter = 0
        merger_cursor.execute("SELECT * FROM node order by name")
        while True:
            row = merger_cursor.fetchone()
            counter += 1
            node_dict = {}
            if row is None:
                break
            # Limit
            # if counter == 1000:
            #     break
            else:
                # Mapping for gene names and non-uniprot protein names
                # Assigning ids that we want to be mapped
                uniprot_id = row['name'].split(':')[1]

                gene_name = map_uniprot_to_gene(uniprot_id)
                protein_name = map_uniprot_to_protein(uniprot_id)
                #prot_ext_name = allrows[1]

                # NODES

                # Creating dictionary format
                # Gene name
                if gene_name:
                    node_dict['displayedName'] = gene_name
                # If uniprot id doesn't map to gene name, just insert uniprot id as full name
                else:
                    node_dict['displayedName'] = uniprot_id

                # Protein name
                if protein_name:
                    node_dict['fullName'] = protein_name
                # If uniprot id doesn't map to protein name, just insert uniprot id as full name
                else:
                    node_dict['fullName'] = uniprot_id

                # Uniprot id
                node_dict['name'] = uniprot_id

                # Taxon
                node_dict['taxon'] = {"name": SPEC_LONG[row['tax_id']],
                                      "shortName": SPEC_MAP[row['tax_id']],
                                      "commonName": SPEC_SHORT[row['tax_id']],
                                      "id": float(row['tax_id'].split(':')[1])}

                # Topology
                node_dict['topologicalFeatures'] = []
                if row['topology'] and row['topology'] != '':
                    for piped in row['topology'].split('|'):
                        if ':' not in piped:
                            topo = piped.strip()
                            node_dict['topologicalFeatures'].append(
                                {
                                    "value": topo,
                                    "db": "",
                                    "url": "",
                                    "searchable": False
                                })

                # Pathway
                node_dict['pathways'] = []
                slk_pathways = {
                    'TCR': 'TCR',
                    'BCR': 'BCR',
                    'TLR': 'TLR',
                    'IIP': 'IIP',
                    'JAK/STAT': 'JAK/STAT',
                    'NHR': 'NHR',
                    'RTK': 'RTK',
                    'Rho/Cytoskeleton': 'Rho/Cytoskeleton',
                    'TGF': 'TGF',
                    'Notch': 'Notch',
                    'GPCR': 'GPCR',
                    'WNT/Wingless': 'WNT/Wingless',
                    'HIPPO': 'HIPPO',
                    'Hippo': 'HIPPO',
                    'HH': 'HH',
                    'TNF/Apoptosis': 'TNF/Apoptosis',
                    'Toll-like receptor': 'TLR',
                    'Receptor tyrosine kinase': 'RTK',
                    'Innate immune pathways': 'IIP',
                    'Hedgehog': 'HH',
                    'T-cell receptor': 'TCR',
                    'B-cell receptor': 'BCR',
                    'TNF pathway': 'TNF',
                    'Nuclear hormone receptor': 'NHR',

                }
                if row['pathways'] is not None:
                    for piped in row['pathways'].split('|'):
                        if ':' in piped:
                            colond = piped.strip().split(':')[1].replace('*', '')
                        else:
                            colond = piped.strip().replace('*', '')

                        if '(' in colond:
                            colond = colond.split('(')[0]

                        if colond in slk_pathways.keys():
                            node_dict['pathways'].append({
                                "value": slk_pathways[colond],
                                "db": "",
                                "url": "",
                                "searchable": False
                            }, )


                # Tissue
                node_dict['tissues'] = []
                if row['topology'] and row['topology'] != '':
                    for piped in row['topology'].split('|'):
                        if 'tissue:' in piped:
                            tissuename = piped.replace('tissue:', '')
                            node_dict['tissues'].append(
                                {
                                "value": tissuename,
                                "db": "Uber Anatomy Ontology",
                                "url": "",
                                "searchable": False
                                })

                # Localization
                node_dict['minorCellularLocalization'] = []
                if row['topology'] and row['topology'] != '':
                    for piped in row['topology'].split('|'):
                        # Getting minorloc data
                        if 'minorloc:' in piped:
                            # Getting each minorloc
                            data = piped.replace('minorloc:', '')
                            node_dict['minorCellularLocalization'].append(
                                {
                                    "value": data,
                                    "db": 'ComPPI',
                                    "url": 'https://comppi.linkgroup.hu/protein_search/interactors/' + uniprot_id,
                                    "searchable": False
                                })

                node_dict['majorCellularLocalization'] = []
                if row['topology'] and row['topology'] != '':
                    for piped in row['topology'].split('|'):
                        # Getting majorloc data
                        if 'majorloc:' in piped:
                            # Getting each majorloc
                            data = piped.replace('majorloc:', '')
                            node_dict['majorCellularLocalization'].append(
                                {
                                    "value": data,
                                    "db": 'ComPPI',
                                    "url": 'https://comppi.linkgroup.hu/protein_search/interactors/' + uniprot_id,
                                    "searchable": False
                                })

                # Drug-target
                node_dict['drug'] = []
                if row['topology'] and row['topology'] != '':
                    for piped in row['topology'].split('|'):
                        # Getting drugbank data
                        if 'drugtarget:' in piped:
                            # Getting each drug
                            data = piped.replace('drugtarget:', '')
                            node_dict['drug'].append(
                                {
                                    "value": data,
                                    "db": 'DrugBank',
                                    "url": 'https://www.drugbank.ca/',
                                    "searchable": False
                                })
                # External ref
                node_dict['externalReferences'] = []

                # Molecule type
                # If RNA
                if 'UR' in row['name']:
                    node_dict['moleculeType'] = [{"value": "MI:0320(ribonucleic acid)",
                                                  "db": "Molecular Interactions Controlled Vocabulary",
                                                  "url": "www.ebi.ac.uk/ols/ontologies/mi/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FMI_0320",
                                                  "searchable": False
                                                  }, ]

                    # adding the RNACentral primary ID
                    node_dict['externalReferences'].append({
                        "value": uniprot_id,
                        "db": "RNACentral",
                        "url": 'rnacentral.org/rna/' + uniprot_id,
                        "searchable": True
                    })
                    # Adding other refs
                    if row['alt_accession']:
                        for piped in row['alt_accession'].split('|'):
                            if ':' in piped:
                                rnaextref = piped.strip().split(':')[1]
                                if '-' in rnaextref:
                                    node_dict['externalReferences'].append({
                                        "value": rnaextref,
                                        "db": "miRBase",
                                        "url": 'www.mirbase.org',
                                        "searchable": True
                                    })
                            else:
                                rnaextref = piped.strip()
                                if '-' in rnaextref:
                                    node_dict['externalReferences'].append({
                                        "value": rnaextref,
                                        "db": "miRBase",
                                        "url": 'www.mirbase.org',
                                        "searchable": True
                                    })

                # If protein
                else:
                    node_dict['moleculeType'] = [{"value": "MI:0326(protein)",
                                                  "db": "Molecular Interactions Controlled Vocabulary",
                                                  "url": "www.ebi.ac.uk/ols/ontologies/mi/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FMI_0326",
                                                  "searchable": False
                                                  }, ]
                    # adding the uniprot primary ID
                    node_dict['externalReferences'].append({
                        "value": uniprot_id,
                        "db": "Uniprot",
                        "url": 'www.uniprot.org/uniprot/' + uniprot_id,
                        "searchable": True
                    })

                    # adding all the aliases based on the retrieved uniprot id
                    # primary uniprot id is mapped back to its foreign ids
                    # db is the foreign id type from mapping table
                    # url is generated based on foreign id type
                    if uniprot_id in uniprot_to_external_id_dict:
                        for key, value in uniprot_to_external_id_dict[uniprot_id].items():

                            extreftype = key
                            extrefidlist = value

                            # Generating urls based on external reference type
                            if extreftype == 'uniprotac':
                                url = 'www.uniprot.org/uniprot/'
                            elif extreftype == 'ensembl':
                                url = 'www.ensembl.org/id/'
                            elif extreftype == 'hgnc':
                                url = 'www.genenames.org/data/gene-symbol-report/#!/hgnc_id/'
                            elif extreftype == 'zfin':
                                url = 'www.zfin.org/'
                            elif extreftype == 'flybase':
                                url = 'www.flybase.org/reports/'
                            elif extreftype == 'wormbase':
                                url = 'wormbase.org/search/gene/'
                            else:
                                url = ''

                            for extrefid in extrefidlist:
                                node_dict['externalReferences'].append({
                                        "value": extrefid,
                                        "db": extreftype,
                                        "url": url + extrefid,
                                        "searchable": True
                                    })

                # Adding all data
                nodes.append(node_dict)

    # Adding data to json files
    with open('nodes.json', 'w') as outfile:
        json.dump(nodes, outfile, indent=4, sort_keys=True)


# EDGES
def get_edge_data(builder_location):
    logging.debug('Started getting edge data')
    edges = []
    builder_conn = sqlite3.connect(builder_location)
    with builder_conn:
        builder_conn.row_factory = sqlite3.Row
        builder_cursor = builder_conn.cursor()
        # Getting data by layer
        for layer in layers_list:
            counter = 0
            builder_cursor.execute("SELECT * FROM %s" % layer)
            while True:
                row = builder_cursor.fetchone()
                counter += 1
                edge_dict = {}
                edge_dict['interactionDetectionMethods'] = []
                edge_dict['confidenceScore'] = []
                edge_dict['publications'] = []
                edge_dict['sourceDatabases'] = []

                # Limit
                # if counter == 1000:
                #     break

                if row is None:
                    break

                else:
                    # Mapping for gene names and non-uniprot protein names
                    # Assigning ids that we want to be mapped
                    uniprot_source = row['interactor_a_node_name'].split(':')[1]
                    uniprot_target = row['interactor_b_node_name'].split(':')[1]

                    # Mapping
                    source_name = map_uniprot_to_gene(uniprot_source)
                    target_name = map_uniprot_to_gene(uniprot_target)

                    # Setting directedness based on direction prediction
                    if row['confidence_scores']:
                        for piped in row['confidence_scores'].split('|'):
                            if 'dir_pred' in piped:
                                dir_data = piped.replace('dir_pred:', '')
                                # If the calculated score is bigger than cutoff, interaction is directed
                                if float(dir_data) > 2:
                                    edge_dict['isDirected'] = True
                                    # Uniprot
                                    edge_dict['source'] = uniprot_source
                                    edge_dict['target'] = uniprot_target
                                # If the calculated score is lower than cutoff,
                                # we change the source and target and set directedness to True
                                elif float(dir_data) < 2:
                                    edge_dict['isDirected'] = True
                                    # Uniprot
                                    edge_dict['source'] = uniprot_target
                                    edge_dict['target'] = uniprot_source
                                else:
                                    # L1, PTMs, TFs, miRNAs and lncRNAs are directed
                                    if layer == 'layer1' or layer == 'layer2' or layer == 'layer5' \
                                            or layer == 'layer6' or layer == 'layer7':
                                        # Uniprot
                                        edge_dict['source'] = uniprot_source
                                        edge_dict['target'] = uniprot_target
                                        edge_dict['isDirected'] = True
                                    else:
                                        # Uniprot
                                        edge_dict['source'] = uniprot_source
                                        edge_dict['target'] = uniprot_target
                                        edge_dict['isDirected'] = False

                    else:
                        # L1, PTMs, TFs, miRNAs and lncRNAs are directed
                        if layer == 'layer1' or layer == 'layer2' or layer == 'layer5' \
                                or layer == 'layer6' or layer == 'layer7':
                            # Uniprot
                            edge_dict['source'] = uniprot_source
                            edge_dict['target'] = uniprot_target
                            edge_dict['isDirected'] = True
                        # everything in L0 is directed except for TCR
                        elif layer == 'layer0':
                            for i in range(len(row['source_db'].split('|'))):
                                sourcedb_raw = row['source_db'].split('|')[i].replace('source database:', '')
                                if sourcedb_raw == 'tcr':
                                    edge_dict['source'] = uniprot_source
                                    edge_dict['target'] = uniprot_target
                                    edge_dict['isDirected'] = False
                                else:
                                    edge_dict['source'] = uniprot_source
                                    edge_dict['target'] = uniprot_target
                                    edge_dict['isDirected'] = True
                        else:
                            # Uniprot
                            edge_dict['source'] = uniprot_source
                            edge_dict['target'] = uniprot_target
                            edge_dict['isDirected'] = False

                    # Name
                    if map_uniprot_to_protein(uniprot_source):
                        edge_dict['sourceFullName'] = map_uniprot_to_protein(uniprot_source)
                    # If uniprot id doesn't map to protein name, just insert uniprot id as full name
                    else:
                        edge_dict['sourceFullName'] = uniprot_source
                    edge_dict['sourceDisplayedName'] = source_name
                    if map_uniprot_to_protein(uniprot_target):
                        edge_dict['targetFullName'] = map_uniprot_to_protein(uniprot_target)
                    # If uniprot id doesn't map to protein name, just insert uniprot id as full name
                    else:
                        edge_dict['targetFullName'] = uniprot_target
                    edge_dict['targetDisplayedName'] = target_name
                    # Directed
                    edge_dict['isDirect'] = True

                    # Interaction types
                    edge_dict['interactionType'] = []

                    if row['interaction_types']:
                        # effect:MI:XXXX(name)|is_direct:MI:XXXX(name)|is_directed:directed
                        for piped in row['interaction_types'].split('|'):
                            # [effect:MI:XXXX(name), is_direct:MI:XXXX(name), is_directed:directed]
                            colond = piped.split(':')
                            # colond[0] = effect or is_direct or is_directed
                            # colond[1] = MI
                            # colond[2] = XXXX(name)

                            if colond[0] == 'effect':
                                # mi_id = XXXX
                                if colond[1] == 'MI':
                                    mi_id = colond[2].split('(')[0]
                                else:
                                    continue
                                # full_mi = MI:XXXX(name)
                                full_mi = colond[1] + ':' + colond[2]

                                edge_dict['interactionType'].append(
                                    {
                                        "value": full_mi,
                                        "db": 'Molecular Interactions Controlled Vocabulary',
                                        "url": "www.ebi.ac.uk/ols/ontologies/mi/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FMI_" + mi_id,
                                        "searchable": False
                                    }
                                )
                            elif colond[0] == 'MI':
                                # mi_id = XXXX
                                mi_id = colond[1].split('(')[0]
                                # full_mi = MI:XXXX(name)
                                full_mi = colond[0] + ':' + colond[1]

                                edge_dict['interactionType'].append(
                                    {
                                        "value": full_mi,
                                        "db": 'Molecular Interactions Controlled Vocabulary',
                                        "url": "www.ebi.ac.uk/ols/ontologies/mi/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FMI_" + mi_id,
                                        "searchable": False
                                    }
                                )
                    # Adding stimulaton to effect of every interaction in TF layer
                    if layer == 'layer6':
                        edge_dict['interactionType'].append(
                            {
                                "value": 'MI:0624(stimulation)',
                                "db": 'Molecular Interactions Controlled Vocabulary',
                                "url": "www.ebi.ac.uk/ols/ontologies/mi/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FMI_" + 'MI:0624(stimulation)',
                                "searchable": False
                            }
                        )
                    # Adding inhibition to effect of every interaction in miRNA, lncRNA layer
                    elif layer == 'layer5' or layer == 'layer7':
                        edge_dict['interactionType'].append(
                            {
                                "value": 'MI:0623(inhibition)',
                                "db": 'Molecular Interactions Controlled Vocabulary',
                                "url": "www.ebi.ac.uk/ols/ontologies/mi/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FMI_" + 'MI:0623(inhibition)',
                                "searchable": False
                            }
                        )
                    # Adding data based on sign prediction
                    if row['confidence_scores']:
                        for piped in row['confidence_scores'].split('|'):
                            if 'sign_pred:' in piped:
                                data = piped.replace('sign_pred:', '')
                                # If calculated sign score is greater than cutoff (1),
                                # the interaction is predicted to be stimulation
                                if float(data) > 1:
                                    sign_data = 'MI:0624(stimulation)'
                                    edge_dict['interactionType'].append(
                                        {
                                            "value": sign_data,
                                            "db": '',
                                            "url": '',
                                            "searchable": False
                                        })
                                # If calculated sign score is lower than cutoff (/1),
                                # the interaction is predicted to be inhibition
                                elif float(data) < -1:
                                    sign_data = 'MI:0623(inhibition)'
                                    edge_dict['interactionType'].append(
                                        {
                                            "value": sign_data,
                                            "db": '',
                                            "url": '',
                                            "searchable": False
                                        })

                    # Interaction detection method
                    if row['interaction_detection_method']:
                        for piped in row['interaction_detection_method'].split('|'):
                            # MI:XXXX(methodname)
                            part = piped.strip().split('(')
                            # part[0] = MI:XXXX
                            # mi_id = XXXX
                            mi_id = part[0].split(':')[1]

                            edge_dict['interactionDetectionMethods'].append(
                                    {
                                        "value": piped,
                                        "db": 'Molecular Interactions Controlled Vocabulary',
                                        "url": "www.ebi.ac.uk/ols/ontologies/mi/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FMI_" + mi_id,
                                        "searchable": False
                                    })

                    # Publications
                    edge_dict['publications'] = []
                    if row['publication_ids']:
                        for piped in row['publication_ids'].split('|'):
                            colond = piped.strip().split(':')
                            if len(colond) != 2:
                                continue
                            else:
                                if colond[0] == 'pubmed':
                                    url = "www.ncbi.nlm.nih.gov/pubmed/%s" % colond[1]
                                else:
                                    url = ""
                                edge_dict['publications'].append(
                                    {
                                        "value": colond[1],
                                        "db": colond[0],
                                        "url": url,
                                        "searchable": False
                                    }
                                )

                    # For each database
                    # Source DB
                    edge_dict['sourceDatabases'] = []
                    # Mapping manually created source db names to displayed names
                    sourcedb_map = {
                        'acsn': 'ACSN',
                        'ACSN': 'ACSN',
                        'innatedb': 'InnateDB',
                        'InnateDB': 'InnateDB',
                        'signor': 'Signor',
                        'reactome': 'Reactome',
                        'slk2': 'SignaLink',
                        'SLKv2.0': 'SignaLink',
                        'SLKv3.0': 'SignaLink',
                        'SLKv2.1': 'SignaLink',
                        'slk3': 'SignaLink',
                        'slk21': 'SignaLink',
                        'tcr': 'SignaLink',
                        'TCRcuration': 'SignaLink',
                        'PSP': 'PSP',
                        'PhosphoSite': 'PhosphoSite',
                        'PTMCode2': 'PTMCode2',
                        'ComPPI': 'ComPPI',
                        'HPRD': 'HPRD',
                        'IntAct': 'IntAct',
                        'biogrid': 'TheBiogrid',
                        'OmniPath': 'OmniPath',
                        'miRDeathDB': 'miRDeathDB',
                        'miRecords': 'miRecords',
                        'miR2Disease': 'miR2Disease',
                        'TarBase': 'TarBase',
                        'starBase': 'StarBase',
                        'TFBS': 'PSSMprediction',
                        'TFBS_Orsi': 'TFlink',
                        'lncRInter': 'lncRInter',
                        'miRSponge': 'miRSponge',
                        'NPInter': 'NPInter',
                        'starbase': 'StarBase',
                        'SignaFish': 'SignaFish',
                        'ELM': 'ELM',
                        'ELM_pred': 'ELM'

                    }
                    for i in range(len(row['source_db'].split('|'))):
                        sourcedb_raw = row['source_db'].split('|')[i].replace('source database:', '')
                        if sourcedb_raw in sourcedb_map.values():
                            edge_dict['sourceDatabases'].append({
                                "value": sourcedb_raw.strip(),
                                "db": "",
                                "url": "",
                                "searchable": False
                            }, )
                        else:
                            edge_dict['sourceDatabases'].append({
                                "value": sourcedb_map[sourcedb_raw],
                                "db": "",
                                "url": "",
                                "searchable": False
                            },)

                    # Layer
                    edge_dict['layer'] = [
                        {
                            "value": float(row['layer']),
                            "db": "",
                            "url": "",
                            "searchable": True
                        }
                    ]

                    if row['confidence_scores'] is not None:
                        for part in row['confidence_scores'].split('|'):
                            if len(part.split(':')) > 1:
                                edge_dict['confidenceScore'].append(
                                    {
                                        "value": part.split(':')[1],
                                        "db": part.split(':')[0],
                                        "url": "",
                                        "searchable": False
                                    }
                                )

                    edges.append(edge_dict)

    #  Adding data to json files
    with open('edges.json', 'w') as outfile2:
        json.dump(edges, outfile2, indent=4, sort_keys=True)


# ATTRIBUTES
def get_attribute_data(builder_location):
    logging.debug('Started getting attribute data')
    # Getting all possibilities
    topology = []
    pathways = []
    tissues = []
    min_loc = []
    maj_loc = []
    illness = []
    cancer = []
    drug = []
    detect = []
    int_type = []
    conf_scores = []

    # Node data
    builder_conn = sqlite3.connect(builder_location)
    with builder_conn:
        builder_conn.row_factory = sqlite3.Row
        builder_cursor = builder_conn.cursor()
        builder_cursor.execute("SELECT * FROM node")
        while True:
            row = builder_cursor.fetchone()
            if row is None:
                break
            else:
                # Assigning variables
                # Pathways
                slk_pathways = {
                    'TCR',
                    'BCR',
                    'TLR',
                    'IIP',
                    'JAK/STAT',
                    'NHR',
                    'RTK',
                    'Rho/Cytoskeleton',
                    'TGF',
                    'Notch',
                    'GPCR',
                    'WNT/Wingless',
                    'HIPPO',
                    'HH',
                    'TNF/Apoptosis',
                }
                if row['pathways']:
                    for piped in row['pathways'].split('|'):
                        if ':' in piped:
                            path_var = piped.strip().split(':')[1].replace('*', '')
                        else:
                            path_var = piped.strip().replace('*', '')

                        if '(' in path_var:
                            path_var = path_var.split('(')[0]

                        if path_var in slk_pathways:
                            if path_var not in pathways:
                                pathways.append(path_var)


                # Tissue data
                if row['topology'] and row['topology'] != '':
                    piped = row['topology'].split('|')
                    for item in piped:
                        if 'tissue' in item.split(':'):
                            tiss_var = ':'.join(item.split(':')[1:])

                            if tiss_var not in tissues:
                                tissues.append(tiss_var)
                        # Subcellular localization
                        elif 'minorloc' in item.split(':'):
                                minloc_var = ':'.join(item.split(':')[1:])
                                if minloc_var not in min_loc:
                                    min_loc.append(minloc_var)

                        elif 'majorloc' in item.split(':'):
                                majloc_var = ':'.join(item.split(':')[1:])
                                if majloc_var not in maj_loc:
                                    maj_loc.append(majloc_var)
                        # Drugtarget
                        elif 'drugtarget' in item.split(':'):
                                drug_var = ':'.join(item.split(':')[1:])
                                if drug_var not in drug:
                                    drug.append(drug_var)
                        else:
                            topol_var = item

                            if topol_var not in topology:
                                topology.append(topol_var)
    # Edge data
    builder_conn = sqlite3.connect(builder_location)
    with builder_conn:
        builder_conn.row_factory = sqlite3.Row
        builder_cursor = builder_conn.cursor()
        counter = 0
        # Getting data by layer
        for layer in layers_list:
            builder_cursor.execute("SELECT * FROM %s" % layer)
            while True:
                row = builder_cursor.fetchone()
                counter += 1
                if row is None:
                    break
                # Limit
                # if counter == 1000:
                #     break
                else:
                    # Source db
                    source_var = row['source_db']

                    # Detection method
                    if row['interaction_detection_method']:
                        for piped in row['interaction_detection_method'].split('|'):
                            # MI:XXXX(methodname)
                            detect_var = piped.strip()
                            if detect_var not in detect:
                                detect.append(detect_var)

                    # Interaction type
                    int_type_var = 'N/A'
                    if row['interaction_types'] is not None:
                        for piped in row['interaction_types'].split('|'):
                            colond = piped.strip().split(':')
                            if len(colond) != 2:
                                continue
                            else:
                                if colond[0] != 'is_directed' or colond[0] != 'is_direct':
                                    int_type_var = colond[1]

                                if int_type_var not in int_type:
                                    int_type.append(int_type_var)

                    # Confidence scores
                    if row['confidence_scores']:
                        for part in row['confidence_scores'].split('|'):
                            if len(part.split(':')) > 1:
                                conf_var = part.split(':')[1]
                                if conf_var not in conf_scores:
                                    conf_scores.append(conf_var)

    attrib = []

    # Topology
    attrib.append({
        "key": "topologicalFeatures",
        "displayedName": "Topological Features",
        "dataType": "string",
        "possibleValues": [],
        "searchable": False,
        "isNode": True
    })
    # Pathways
    attrib.append({
        "key": "pathways",
        "displayedName": "Pathways",
        "dataType": "string",
        "possibleValues": pathways,
        "searchable": False,
        "isNode": True
    })
    # Tissues
    attrib.append({
        "key" : "tissues",
        "displayedName" : "Tissues",
        "dataType" : "string",
        "possibleValues" : tissues,
        "searchable" : False,
        "isNode" : True
    })
    # Minor cellular localization
    attrib.append({
        "key": "minorCellularLocalization",
        "displayedName": "Minor Cellular Localization",
        "dataType": "string",
        "possibleValues": min_loc,
        "searchable": False,
        "isNode": True
    })
    # Major cellular localization
    attrib.append({
        "key": "majorCellularLocalization",
        "displayedName": "Major Cellular Localization",
        "dataType": "string",
        "possibleValues": maj_loc,
        "searchable": False,
        "isNode": True
    })
    # Molecule type
    attrib.append({
        "key": "moleculeType",
        "displayedName": "Molecule Type",
        "dataType": "string",
        "possibleValues": [
            'MI:0326(protein)',
            'MI:0320(ribonucleic acid)',
        ],
        "searchable": False,
        "isNode": True
    })
    # Illness
    attrib.append({
        "key": "illness",
        "displayedName": "Illness",
        "dataType": "string",
        "possibleValues": illness,
        "searchable": False,
        "isNode": True
    })
    # Cancer
    attrib.append({
        "key": "cancer",
        "displayedName": "Cancer",
        "dataType": "string",
        "possibleValues": cancer,
        "searchable": False,
        "isNode": True
    })
    # Drug
    attrib.append({
        "key": "drug",
        "displayedName": "Drug",
        "dataType": "string",
        "possibleValues": drug,
        "searchable": False,
        "isNode": True
    })
    # Reference
    attrib.append({
        "key": "externalReferences",
        "displayedName": "External References",
        "dataType": "string",
        "possibleValues": [
            "uniprotac",
            "ensembl",
            "hgnc",
            "miRBase"
        ],
        "searchable": True,
        "isNode": True
    })
    # Interaction detection method
    attrib.append({
        "key": "interactionDetectionMethods",
        "displayedName": "Interaction Detection Methods",
        "dataType": "string",
        "possibleValues": detect,
        "searchable": False,
        "isNode": False
    })
    # Publications
    attrib.append({
        "key": "publications",
        "displayedName": "Publications",
        "dataType": "number",
        "possibleValues": [],
        "searchable": False,
        "isNode": False
    })
    # Source db
    attrib.append({
        "key": "sourceDatabases",
        "displayedName": "Source Databases",
        "dataType": "string",
        "possibleValues": [
            'ACSN',
            'InnateDB',
            'Reactome',
            'Signor',
            'Manual curation',
            'PSP',
            'ELM',
            'PhosphoSite',
            'BioGRID',
            'ComPPI',
            'HPRD',
            'IntAct',
            'OmniPath',
            'miR2Disease',
            'miRDeathDB',
            'miRecords',
            'starBase',
            'TarBase',
            'lncRInter',
            'miRSponge',
            'NPInter',
            'SignaFish'
        ],
        "searchable": False,
        "isNode": False
    })
    # Layer
    attrib.append({
        "key": "layer",
        "displayedName": "Layer",
        "dataType": "number",
        "possibleValues": [
            0,
            1,
            2,
            3,
            5,
            6,
            7
        ],
        "searchable": False,
        "isNode": False
    })
    # Interaction type
    attrib.append({
        "key": "interactionType",
        "displayedName": "Interaction Type",
        "dataType": "string",
        "possibleValues": int_type,
        "searchable": False,
        "isNode": False
    })
    # Confidence score
    attrib.append({
        "key": "confidenceScore",
        "displayedName": "Confidence Score",
        "dataType": "number",
        "possibleValues": [],
        "searchable": False,
        "isNode": False
    })
    # Is direct
    attrib.append({
        "key": "isDirect",
        "displayedName": "Is Direct",
        "dataType": 'undefined',
        "possibleValues": [True, False],
        "searchable": False,
        "isNode": False
    })

    # Adding data to json files
    with open('attributes.json', 'w') as outfile2:
        json.dump(attrib, outfile2, indent=4, sort_keys=True)

logging.debug('Done sorting')
