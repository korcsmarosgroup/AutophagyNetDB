"""
Use case 1: Select proteins from ATG core annotated to xenophagy - select their connections in core + their regulators

Use case 2: ATG regulation in 3 tissues - intestine, skin and brain

"""

import csv
import sqlite3
import pandas as pd

# UC1
xeno_list = []
# Open ATG type file and select xenophagy proteins
with open('../../ARNlib/GO-AP-and-regulators.csv') as infile:
    infile.readline()
    for line in infile:
        line = line.strip().split(';')
        if line[1] == 'xenophagy':
            xeno_list.append(line[0])

# Get interactions of xenophagy proteins in core
xeno_int_list = []
with open('ARN2_core.csv') as corefile:
    header = corefile.readline()
    for line in corefile:
        line= line.strip().split(',')
        source = line[0].split(':')[1]
        target = line[1].split(':')[1]
        for xenoprot in xeno_list:
            if source == xenoprot or target == xenoprot:
                xeno_int_list.append(line)

with open('xeno_core.csv', 'w+', newline='') as outfile:
    outfile.write(header)
    write = csv.writer(outfile)
    write.writerows(xeno_int_list)

conn = sqlite3.connect("XENO_layers.db") # change to 'sqlite:///your_filename.db'
cur = conn.cursor()
# cur.execute("CREATE TABLE layer0 (interactor_a_node_name,interactor_b_node_name,interaction_detection_method,first_author,publication_ids,interaction_types,source_db,interaction_identifiers,confidence_scores,layer);") # use your column names here

df = pd.read_csv('xeno_core.csv')
df.to_sql('layer0', conn, if_exists='append', index=False)

def main(log, path):
    # Creating output table for EDGES
    # path = '../DATA/workflow/merger.db'

    conn = sqlite3.connect('XENO_layers.db')
    # log.debug("Started connection to '%s'" % conn)
    with conn:
        c = conn.cursor()
        for layer in [1, 2, 3, 5, 6, 7]:
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
        l3_nodes = set()
        l0_nodes = set()

        conn_xeno = sqlite3.connect('XENO_layers.db')
        conn_xeno.row_factory = sqlite3.Row
        cur_xeno = conn_xeno.cursor()
        cur_xeno.execute('SELECT * FROM layer0 ORDER BY layer')

        while True:
            line_xeno = cur_xeno.fetchone()
            # Until the last row
            # For just 100 edges from each database
            if line_xeno is None:
                print("building finished successfully! Please find the final statistics below:")
                print("Edges processed in layer 0: %d" % layer_counter[0])
                break
            else:
                layer = line_xeno['layer']
                layer_counter[layer] += 1
                if layer == 0:
                    nodes_added_in_current_layer.add(line_xeno['interactor_a_node_name'])
                    nodes_added_in_current_layer.add(line_xeno['interactor_b_node_name'])

        l0_nodes.update(nodes_added_in_current_layer)
        nodes_added_in_current_layer = set()

        while True:
            line = c1.fetchone()
            # Until the last row
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
                if line['layer'] > 0:

                    with conn:
                        c = conn.cursor()

                        # Extracting L1 nodes
                        layer = line['layer']
                        layer_counter[layer] += 1

                        if layer > current_layer:
                            # we are just starting to process a new layer,
                            # let's store all the new nodes from the layer we just finished
                            nodes_added_in_previous_layers.update(nodes_added_in_current_layer)
                            nodes_added_in_current_layer = set()
                            current_layer = layer

                        if layer == 3:
                            if line['interactor_a_node_name'] in nodes_added_in_previous_layers or line[
                                'interactor_b_node_name'] in nodes_added_in_previous_layers:
                                layer_remaining_counter[layer] += 1
                                l3_nodes.add(line['interactor_a_node_name'])
                                l3_nodes.add(line['interactor_b_node_name'])
                                insert_new_edge(c, line)

                        else:
                            if line['interactor_a_node_name'] in l0_nodes or line['interactor_b_node_name'] in l0_nodes:
                                layer_remaining_counter[layer] += 1
                                nodes_added_in_current_layer.add(line['interactor_a_node_name'])
                                nodes_added_in_current_layer.add(line['interactor_b_node_name'])
                                insert_new_edge(c, line)

        if layer == 3:
            # also store the nodes from the last layer
            c1.execute('SELECT * FROM node ORDER BY name')
            number_of_nodes_added = 0
            while True:
                line = c1.fetchone()

                if line is None:
                    print("total number of unique nodes found during build: %d" % len(l3_nodes))
                    print("total number of unique nodes added to the output file: %d" % number_of_nodes_added)
                    break

                if line['name'] in l3_nodes:
                    insert_new_node(c, line)
                    number_of_nodes_added += 1
        else:
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
    main(log=None, path='merger.db')
