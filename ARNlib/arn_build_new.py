"""
    This script builds the multi-structured database. Each node is present in every layer.
    The higher layer has the edges from the previous layer and their first neighbours.
"""

import sqlite3
from collections import OrderedDict as od

# Building mock data set
DB = od([('signor', 'SLK_Core'),
         ('PSP', 'layer1'),
         ('PhosphoSite', 'PTM'),
         ('OmniPath', 'ATG_Reg'),
         ('miRDeathDB', 'miRNA'),
         ('starbase', 'lncRNA')])


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


def build_base(log, merger_path):
    """
    Builds direct regulator layer of ARN
    Searches for connections between ATG proteins and other proteins
    :param log: logger
    :param merger_path: location of merger database
    :return: adds connections to layer1 of builder db
    """

    # Creating output table for EDGES
    # merger_path = '../DATA/workflow/merger.db'

    build_conn = sqlite3.connect('ARN_layers.db')
    # log.debug("Started connection to '%s'" % build_conn)
    with build_conn:
        build_cur = build_conn.cursor()
        for layer in [0, 1, 2, 3, 5, 6, 7, 8]:
            build_cur.execute("DROP TABLE IF EXISTS layer%d" % layer)
            build_cur.execute("DROP TABLE IF EXISTS node")

            build_cur.execute('''
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

            build_cur.execute('''
            CREATE TABLE node (
        	name	TEXT NOT NULL,
        	alt_accession	TEXT,
        	tax_id	INTEGER NOT NULL,
        	pathways	TEXT,
        	aliases TEXT,
        	topology TEXT )
            ''')

        merge_conn = sqlite3.connect(merger_path)
        merge_conn.row_factory = sqlite3.Row

        nodes_added_in_L0 = set()
        all_nodes_to_add = set()

        counter = 0
        layer_counter = {0: 0, 1: 0, 2: 0, 3: 0, 5: 0, 6: 0, 7: 0, 8: 0}
        layer_remaining_counter = {1: 0, 2: 0, 3: 0, 5: 0, 6: 0, 7: 0, 8: 0}
        current_layer = 0

        layers = [0, 1, 2, 3]

        for eachlayer in layers:

            merge_cur = merge_conn.cursor()
            merge_cur.execute('SELECT * FROM edge ORDER BY layer, id')

            while True:
                merge_line = merge_cur.fetchone()
                counter += 1
                # Until the last row
                # For just 100 edges from each database
                if merge_line is None:
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
                    with build_conn:
                        # Extracting the data if using mock data files
                        # with open('mockdata.tsv') as mockfile:

                        # Skipping the header
                        # mockfile.readline()
                        # for merge_line in mockfile:
                        #   merge_line = merge_line.strip().split('\t')

                        build_cur = build_conn.cursor()
                        # Extracting L0 nodes

                        layer = merge_line['layer']
                        layer_counter[layer] += 1

                        if layer == 0:
                            # Storing L0 nodes (ATG protieins) in a set
                            nodes_added_in_L0.add(merge_line['interactor_a_node_name'])
                            nodes_added_in_L0.add(merge_line['interactor_b_node_name'])

                            # Also adding nodes to an all nodes set
                            all_nodes_to_add.add(merge_line['interactor_a_node_name'])
                            all_nodes_to_add.add(merge_line['interactor_b_node_name'])
                            insert_new_edge(build_cur, merge_line)

                        else:
                            # If node is ATG protein we check all protein (Layers for interactions)
                            # These will be L1 interactions - direct regulators
                            # USE THIS AFTER ONLY L0 IS PARSED
                            if merge_line['interactor_a_node_name'] in nodes_added_in_L0:
                                # if 'A' node is the ATG protein, add its partner as a new node
                                all_nodes_to_add.add(merge_line['interactor_b_node_name'])
                                layer_remaining_counter[layer] += 1
                                # Add the edge to layer 1
                                layer = 1
                                insert_new_edge(build_cur, merge_line)

                            elif merge_line['interactor_b_node_name'] in nodes_added_in_L0:
                                # if 'B' node is the ATG protein, add its partner as a new node

                                all_nodes_to_add.add(merge_line['interactor_a_node_name'])
                                layer_remaining_counter[layer] += 1
                                # Add the edge to layer 1
                                layer = 1
                                insert_new_edge(build_cur, merge_line)

        # Adding nodes
        with merge_conn:
            merge_cur = merge_conn.cursor()
            merge_cur.execute('SELECT * FROM node ORDER BY name')
            number_of_nodes_added = 0
            while True:
                node_line = merge_cur.fetchone()

                if node_line is None:
                    print("total number of unique nodes found during build: %d" % len(all_nodes_to_add))
                    print("total number of unique nodes added to the output file: %d" % number_of_nodes_added)
                    break

                if node_line['name'] in all_nodes_to_add:
                    insert_new_node(build_cur, node_line)
                    number_of_nodes_added += 1


def build_whole(log, merger_path):
    """
    Searches for connections between direct regulators of ATG and other regulators
    :param log: logger
    :param merger_path: location of merger db
    :return: add connections to builder db
    """

    # Creating output table for EDGES
    # merger_path = '../DATA/workflow/merger.db'

    build_conn = sqlite3.connect('ARN_layers.db')
    merge_conn = sqlite3.connect(merger_path)
    merge_conn.row_factory = sqlite3.Row
    # log.debug("Started connection to '%s'" % build_conn)

    counter = 0
    nodes_added_in_current_layer = set()
    all_nodes_to_add = set()
    layer_counter = {0: 0, 1: 0, 2: 0, 3: 0, 5: 0, 6: 0, 7: 0, 8: 0}
    layer_remaining_counter = {1: 0, 2: 0, 3: 0, 5: 0, 6: 0, 7: 0, 8: 0}
    layers = [2, 3, 5, 6, 7]

    dirreg_dict = {}

    with build_conn:
        build_cur = build_conn.cursor()
        # Select all direct regulators from builder table
        build_cur.execute('SELECT * FROM layer1 ORDER BY layer')
        while True:
            dirreg_row = build_cur.fetchone()
            if dirreg_row is None:
                break
            else:
                # Filling up a dictionary where the keys are the node names and the
                # values are the edge lines

                dirreg_dict[dirreg_row['interactor_a_node_name']] = dirreg_row
                dirreg_dict[dirreg_row['interactor_b_node_name']] = dirreg_row

    for eachlayer in layers:
        with merge_conn:
            merge_cur = merge_conn.cursor()
            # Select everything from merger
            merge_cur.execute('SELECT * FROM edge ORDER BY layer, id')
            while True:
                merger_line = merge_cur.fetchone()
                counter += 1
                # Until the last row
                if merger_line is None:
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
                    current_layer = eachlayer
                    layer = merger_line['layer']
                    layer_counter[layer] += 1

                    if layer > current_layer:
                        # we are just starting to process a new layer,
                        # let's store all the new nodes from the layer we just finished
                        all_nodes_to_add.update(nodes_added_in_current_layer)
                        nodes_added_in_current_layer = set()
                        current_layer = layer

                    if merger_line['interactor_a_node_name'] in dirreg_dict.keys():
                        layer_remaining_counter[layer] += 1
                        nodes_added_in_current_layer.add(merger_line['interactor_a_node_name'])
                        nodes_added_in_current_layer.add(merger_line['interactor_b_node_name'])
                        insert_new_edge(build_cur, merger_line)

                    elif merger_line['interactor_b_node_name'] in dirreg_dict.keys():
                        layer_remaining_counter[layer] += 1
                        nodes_added_in_current_layer.add(merger_line['interactor_a_node_name'])
                        nodes_added_in_current_layer.add(merger_line['interactor_b_node_name'])
                        insert_new_edge(build_cur, merger_line)

    # Adding nodes
    with merge_conn:
        merge_cur = merge_conn.cursor()
        merge_cur.execute('SELECT * FROM node ORDER BY name')
        number_of_nodes_added = 0
        while True:
            node_line = merge_cur.fetchone()

            if node_line is None:
                print("total number of unique nodes found during build: %d" % len(all_nodes_to_add))
                print("total number of unique nodes added to the output file: %d" % number_of_nodes_added)
                break

            if node_line['name'] in all_nodes_to_add:
                insert_new_node(build_cur, node_line)
                number_of_nodes_added += 1


def build_pth_conns(log, merger_path):
    """
    Searches for connections between SLK core layer and ARN regulators, creating pathway connections
    :param log: logger
    :param merger_path: location of merger db file
    :return: dd connections to builder db
    """

    # Creating output table for EDGES
    # merger_path = '../DATA/workflow/merger.db'

    build_conn = sqlite3.connect('ARN_layers.db')
    merge_conn = sqlite3.connect(merger_path)
    merge_conn.row_factory = sqlite3.Row
    # log.debug("Started connection to '%s'" % build_conn)

    counter = 0
    nodes_added_in_current_layer = set()
    all_nodes_to_add = set()
    layer_counter = {0: 0, 1: 0, 2: 0, 3: 0, 5: 0, 6: 0, 7: 0, 8: 0}
    layer_remaining_counter = {1: 0, 2: 0, 3: 0, 5: 0, 6: 0, 7: 0, 8: 0}
    layers = [8]

    int_dict = {}
    with build_conn:
        build_cur = build_conn.cursor()
        # Select all edges from builder
        build_cur.execute('''SELECT * FROM layer1,
                                           layer2,
                                           layer3,
                                           layer5,
                                           layer6,
                                           layer7
                            ORDER BY layer, id''')
        while True:
            int_row = build_cur.fetchone()
            if int_row is None:
                break
            else:
                # Filling up a dictionary where the keys are the node names and the
                # values are the edge lines
                int_dict[int_row['interactor_a_node_name']] = int_row
                int_dict[int_row['interactor_b_node_name']] = int_row

    for eachlayer in layers:
        with merge_conn:
            merge_cur = merge_conn.cursor()
            # Select everything from merger
            merge_cur.execute('SELECT * FROM edge ORDER BY layer, id')
            while True:
                merger_line = merge_cur.fetchone()
                counter += 1
                # Until the last row
                if merger_line is None:
                    print("building finished successfully! Please find the final statistics below:")
                    print("Edges processed in layer 0: %d" % layer_counter[0])
                    print("Edges processed in layer 1: %d (remaining: %d)" % (
                        layer_counter[1], layer_remaining_counter[1]))
                    print("Edges processed in layer 2: %d (remaining: %d)" % (
                        layer_counter[2], layer_remaining_counter[2]))
                    print("Edges processed in layer 3: %d (remaining: %d)" % (
                        layer_counter[3], layer_remaining_counter[3]))
                    print("Edges processed in layer 5: %d (remaining: %d)" % (
                        layer_counter[5], layer_remaining_counter[5]))
                    print("Edges processed in layer 6: %d (remaining: %d)" % (
                        layer_counter[6], layer_remaining_counter[6]))
                    print("Edges processed in layer 7: %d (remaining: %d)" % (
                        layer_counter[7], layer_remaining_counter[7]))
                    break

                else:
                    current_layer = eachlayer
                    layer = merger_line['layer']
                    layer_counter[layer] += 1

                    if merger_line['interactor_a_node_name'] in int_dict.keys():
                        layer_remaining_counter[layer] += 1
                        nodes_added_in_current_layer.add(merger_line['interactor_a_node_name'])
                        nodes_added_in_current_layer.add(merger_line['interactor_b_node_name'])

                        all_nodes_to_add.add(merger_line['interactor_b_node_name'])

                        insert_new_edge(build_cur, merger_line)

                    elif merger_line['interactor_b_node_name'] in int_dict.keys():
                        layer_remaining_counter[layer] += 1
                        nodes_added_in_current_layer.add(merger_line['interactor_a_node_name'])
                        nodes_added_in_current_layer.add(merger_line['interactor_b_node_name'])

                        all_nodes_to_add.add(merger_line['interactor_a_node_name'])

                        insert_new_edge(build_cur, merger_line)

    # Adding nodes
    with merge_conn:
        merge_cur = merge_conn.cursor()
        merge_cur.execute('SELECT * FROM node ORDER BY name')
        number_of_nodes_added = 0
        while True:
            node_line = merge_cur.fetchone()

            if node_line is None:
                print("total number of unique nodes found during build: %d" % len(all_nodes_to_add))
                print("total number of unique nodes added to the output file: %d" % number_of_nodes_added)
                break

            if node_line['name'] in all_nodes_to_add:
                insert_new_node(build_cur, node_line)
                number_of_nodes_added += 1


if __name__ == '__main__':
    build_base(log=None, merger_path='../DATA/workflow/merger.db')
    build_whole(log=None, merger_path='../DATA/workflow/merger.db')
    build_pth_conns(log=None, merger_path='../DATA/workflow/merger.db')
