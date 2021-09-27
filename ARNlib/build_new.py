"""
    This script builds the multi-structured database. Each node is present in every layer.
    The higher layer has the edges from the previous layer and their first neighbours.
"""

import sqlite3
from collections import OrderedDict as od

# Building mock data set
DB = od([('signor', 'layer0'),
         ('PSP', 'layer1'),
         ('PhosphoSite', 'layer2'),
         ('OmniPath', 'layer3'),
         ('miRDeathDB', 'layer5'),
         ('starbase', 'layer7')])


def get_mockdata():
    with open('mockdata.tsv', 'w') as mockfile:
        # Adding header
        mockfile.write(
            'id\tinteractor_node_a_id\tinteractor_node_b_id\tinteractor_node_a_name\tinteractor_node_b_name\tinteraction_detection_method\tfirst_author\tpublication_ids\tinteraction_types\tsource_db\tinteraction_identifier\tconfidence_score\tlayer\n')
        path = 'merger.db'

        conn1 = sqlite3.connect(path)
        c1 = conn1.cursor()
        c1.execute('SELECT * FROM edge ORDER BY layer ASC')
        counter = 0
        while True:
            edge = c1.fetchone()
            counter += 1
            # Until the last row
            # For just 100 edges from each database
            if edge is None:
                break
            mockfile.write("%s\n" % '\t'.join(map(str, edge)))


# Uncomment if working with mock data
# get_mockdata()


def main(log, path):
    # Creating output table for EDGES
    # path = '../DATA/workflow/merger.db'

    conn = sqlite3.connect('SLK3_layers.db')
    # log.debug("Started connection to '%s'" % conn)
    with conn:
        c = conn.cursor()
        for layer in [0, 1, 2, 3, 5, 6, 7]:
            c.execute("DROP TABLE IF EXISTS layer%d" % layer)
            c.execute("DROP TABLE IF EXISTS node")

            c.execute('''
            CREATE TABLE layer%d (
            `interactor_a_node_name`	TEXT NOT NULL,
            `interactor_b_node_name`	TEXT NOT NULL,
            `interaction_detection_method`	TEXT,
            `first_author`	TEXT,
            `publication_ids`	TEXT NOT NULL,
            `interaction_types`	TEXT,
            `source_db`	TEXT NOT NULL,
            `interaction_identifiers`	TEXT,
            `confidence_scores`	TEXT,
            `layer` INTEGER NOT NULL)
            ''' % layer)

            c.execute('''
            CREATE TABLE node (
        	name	TEXT NOT NULL,
        	alt_accession	TEXT,
        	tax_id	INTEGER NOT NULL,
        	pathways	TEXT,
        	aliases TEXT,
        	topology TEXT )
            ''')

        conn1 = sqlite3.connect(path)
        conn1.row_factory = sqlite3.Row
        c1 = conn1.cursor()
        c1.execute('SELECT * FROM edge ORDER BY layer, id')
        counter = 0
        nodes_added_in_current_layer = set()
        nodes_added_in_previous_layers = set()
        layer_counter = {0: 0, 1: 0, 2: 0, 3: 0, 5: 0, 6: 0, 7: 0}
        layer_remaining_counter = {1: 0, 2: 0, 3: 0, 5: 0, 6: 0, 7: 0}
        current_layer = 0

        while True:
            line = c1.fetchone()
            counter += 1
            # Until the last row
            # For just 100 edges from each database
            if line is None:
                print("building finished successfully! Please find the final statistics below:")
                print("Edges processed in layer 0: %d" % layer_counter[0])
                print("Edges processed in layer 1: %d (remaining: %d)" % (layer_counter[1], layer_remaining_counter[1]))
                print("Edges processed in layer 2: %d (remaining: %d)" % (layer_counter[2], layer_remaining_counter[2]))
                print("Edges processed in layer 3: %d (remaining: %d)" % (layer_counter[3], layer_remaining_counter[3]))
                print("Edges processed in layer 5: %d (remaining: %d)" % (layer_counter[5], layer_remaining_counter[5]))
                print("Edges processed in layer 6: %d (remaining: %d)" % (layer_counter[6], layer_remaining_counter[6]))
                print("Edges processed in layer 7: %d (remaining: %d)" % (layer_counter[7], layer_remaining_counter[7]))
                break
            else:
                with conn:
                    # Extracting the data if using mock data files
                    # with open('mockdata.tsv') as mockfile:

                    # Skipping the header
                    # mockfile.readline()
                    # for line in mockfile:
                    #   line = line.strip().split('\t')
                    c = conn.cursor()
                    # Extracting L0 nodes

                    layer = line['layer']
                    layer_counter[layer] += 1

                    if layer > current_layer:
                        # we are just starting to process a new layer,
                        # let's store all the new nodes from the layer we just finished
                        nodes_added_in_previous_layers.update(nodes_added_in_current_layer)
                        nodes_added_in_current_layer = set()
                        current_layer = layer

                    if layer == 0:
                        nodes_added_in_current_layer.add(line['interactor_a_node_name'])
                        nodes_added_in_current_layer.add(line['interactor_b_node_name'])
                        insert_new_edge(c, line)

                    else:
                        if line['interactor_a_node_name'] in nodes_added_in_previous_layers or line['interactor_b_node_name'] in nodes_added_in_previous_layers:
                            layer_remaining_counter[layer] += 1
                            nodes_added_in_current_layer.add(line['interactor_a_node_name'])
                            nodes_added_in_current_layer.add(line['interactor_b_node_name'])
                            insert_new_edge(c, line)


        # also store the nodes from the last layer
        nodes_added_in_previous_layers.update(nodes_added_in_current_layer)

        c1.execute('SELECT * FROM node ORDER BY name')
        number_of_nodes_added = 0
        while True:
            line = c1.fetchone()

            if line is None:
                print("total number of unique nodes found during build: %d" % len(nodes_added_in_previous_layers))
                print("total number of unique nodes added to the output file: %d" % number_of_nodes_added)
                break

            if line['name'] in nodes_added_in_previous_layers:
                insert_new_node(c, line)
                number_of_nodes_added += 1



def insert_new_edge(c, edge_dict):

    insert_query = '''
    INSERT INTO layer%d (
        interactor_a_node_name, 
        interactor_b_node_name, 
        interaction_detection_method, 
        first_author, 
        publication_ids, 
        interaction_types, 
        source_db, 
        interaction_identifiers, 
        confidence_scores, 
        layer)
    VALUES (?,?,?,?,?,?,?,?,?,?)
        ''' % edge_dict['layer']

    c.execute(insert_query, (
        edge_dict['interactor_a_node_name'],
        edge_dict['interactor_b_node_name'],
        edge_dict['interaction_detection_method'],
        edge_dict['first_author'],
        edge_dict['publication_ids'],
        edge_dict['interaction_types'],
        edge_dict['source_db'],
        edge_dict['interaction_identifiers'],
        None,
        edge_dict['layer']
    ))

def insert_new_node(c, node_dict):
    insert_query = 'INSERT INTO node ( name, alt_accession, tax_id, pathways, aliases, topology ) VALUES (?,?,?,?,?,?)'

    c.execute(insert_query, (
        node_dict['name'],
        node_dict['alt_accession'],
        node_dict['tax_id'],
        node_dict['pathways'],
        node_dict['aliases'],
        node_dict['topology']
    ))


if __name__ == '__main__':
    main(log=None, path='../DATA/workflow/merger.db')
