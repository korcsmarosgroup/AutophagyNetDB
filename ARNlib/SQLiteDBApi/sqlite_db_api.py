# -*- coding: utf-8 -*-

__author__ = 'blaise'

import sqlite3, sys


class PsimiSQL:

    def __init__(self, sql_seed_file_location):
        """
        The constuctor of the PsimiToSQL class takes one argument. This will initiate a PsimiToSQL object with SQLite db with the given table structure.
        :param sql_seed_file_location: The location of the SQL script that will be the basis of the db's schema
        :type sql_seed_file_location: str
        :return:
        """

        self.sql_seed = open(sql_seed_file_location).read()

        self.db = self.create_db(":memory:")
        self.cursor = self.db.cursor()

    def sort_attributes(self, string, separator='|', force_null_if_empty = True):
        """
        making easier the testing of the data, we order the attributes case insensitive
        """
        if string is None:
            string = ""
        output_list = string.split(separator)
        output_list = list(filter(lambda item: item != '-' and len(item.strip()) > 0, output_list))
        output_list.sort(key=lambda s: s.lower())
        result = separator.join(output_list)

        if len(result) == 0 and force_null_if_empty:
            return None
        return result

    def merge_attributes(self, string1, string2, separator='|'):
        """
        merge two attribute list strings
        """
        if string1 is None:
            string1 = ""
        if string2 is None:
            string2 = ""
        output_list = string1.split(separator)
        output_list.extend(string2.split(separator))
        output_list = list(set(output_list))
        result = separator.join(output_list)

        if len(result) == 0:
            return None
        return result

    def import_from_db_file(self, db_file_location):
        """
        This function imports the data contained in a given db files. The source files will not be modified.
        :param db_file_location: The absolute or relative location of the db files
        :type db_file_location: str
        :return:
        """

        temporary_db = sqlite3.connect(db_file_location)
        temporary_db_name = "file_db"

        params = (db_file_location, temporary_db_name)
        self.db.execute("ATTACH ? AS ? ", params)
        self.db.commit()

        self.db.execute("INSERT INTO node SELECT * FROM %s.node" % temporary_db_name)
        self.db.commit()

        self.db.execute("INSERT INTO edge SELECT * FROM %s.edge" % temporary_db_name)
        self.db.commit()

        temporary_db.close()

        self.db.execute("DETACH DATABASE %s" % temporary_db_name)
        self.db.commit()

    def create_db(self, location):
        """
        This method creates a db with the schema that was given when the class was initiated
        :param location: The location of the db that will be created. It can be a files path, or with the ":memory:" parameter it will reside in the memory
        :type location: str
        :return: A sqlite3 database object
        """
        db = sqlite3.connect(location)

        db.text_factory = str

        create_tables_query = self.sql_seed

        db.executescript(create_tables_query)
        db.commit()

        return db

    def insert_node(self, node_dict):
        """
        This function inserts a node to the node SQL table. And adds row id property to the node dict
        :param db: A reference to the db that contains the node table where the insertion will take place
        :type db: object
        :param node_dict: A dictionary that contains the node's data
        :type node_dict: dict
        :return:
        """
        # TODO: do not mutate the row dict!!!! return the row id instead. All usages must be refactored

        #Passing the db is essentially for the db.commit() function, but cursor also needs to be created because
        #the last row id is only set when a specified instance of cursor.execute() is used

        if (not 'id' in node_dict) and (not self.get_node(node_dict['name'], node_dict['tax_id'])):

            if not 'topology' in node_dict:
                node_dict['topology'] = ""

            query = "INSERT INTO node (name, alt_accession, tax_id, pathways, aliases, topology) VALUES (?,?,?,?,?,?)"

            aliases = self.sort_attributes(self.merge_attributes(node_dict['alt_accession'], node_dict['aliases']))
            self.cursor.execute(query, (
                node_dict['name'],
                aliases, # alt_accession
                node_dict['tax_id'],
                self.sort_attributes(node_dict['pathways']),
                aliases, # aliases
                self.sort_attributes(node_dict['topology'])))
            self.db.commit()

            node_dict['id'] = self.cursor.lastrowid

        elif (not 'id' in node_dict) and (self.get_node(node_dict['name'], node_dict['tax_id'])):

            node_dict['id'] = self.get_node(node_dict['name'], node_dict['tax_id'])['id']

    def insert_unique_node(self, node_dict):
        """
        This function inserts a unique node to the node SQL table. This method uses less resources, because it does not look up
        wether the entry is in the db already. If the given node is already in the db, it will be duplicated. And adds row id property to the given node dict
        :param db: A reference to the db that contains the node table where the insertion will take place
        :type db: object
        :param node_dict: A dictionary that contains the node's data
        :type node_dict: dict
        :return:
        """
        #Passing the db is essentially for the db.commit() function, but cursor also needs to be created because
        #the last row id is only set when a specified instance of cursor.execute() is used

        if not 'topology' in node_dict:
            node_dict['topology'] = ""

        query = "INSERT INTO node (name, alt_accession, tax_id, pathways, aliases, topology) VALUES (?,?,?,?,?,?)"

        aliases = self.sort_attributes(self.merge_attributes(node_dict['alt_accession'], node_dict['aliases']))
        self.cursor.execute(query, (
            node_dict['name'],
            aliases, # alt_accession
            node_dict['tax_id'],
            self.sort_attributes(node_dict['pathways']),
            aliases, # aliases
            self.sort_attributes(node_dict['topology'])))
        self.db.commit()

        node_dict['id'] = self.cursor.lastrowid


    def get_node(self,node_name, node_tax_id):
        """
        This function returns a node from the objects db as a node dict.
        The returned object contains an extra id property (that is the SQLite id of the node).
        If there is result of the query is an empty tuple, than None is returned.
        :param node_name: The name of the node in Psi-Mi format: db:name eg.: hgnc:APAF1
        :type node_name: str
        :param node_tax_id: The taxonomy id of the node in Psi-Mi format: taxid:number eg.: taxid:9606
        :return: A node dict object, with a extra id property
        """
        query = "SELECT * FROM node WHERE name = ? AND tax_id = ?"
        tup = (node_name, node_tax_id)

        self.cursor.execute(query,tup)
        self.db.commit()

        answer = self.cursor.fetchone()

        if answer:
            node_dict = {
                'id' : answer[0],
                'name' : answer[1], #primary id of the database
                'alt_accession' : answer[2], #other ids in mi: format
                'tax_id' : answer[3],
                'pathways' : answer[4], #only in 0-1st layers
                'aliases' : answer[5],  #gene name that is not an id
                'topology' : answer[6] #only in 0-1st layers
            }
            return node_dict
        else:
            return None

    def get_node_by_alt_acc(self,alt_accession):
        """
        Returns the first match for the specified alt accession. If there is no match a None type object is returned.
        :param alt_accession: The alternative accession (gene symbol, short label...) of the node.
        :type alt_accession: str
        :return:
        """
        query = "SELECT * FROM node WHERE alt_accession LIKE '%" + alt_accession + "%'"


        self.cursor.execute(query)
        self.db.commit()

        answer = self.cursor.fetchone()

        if answer:
            node_dict = {
                'id' : answer[0],
                'name' : answer[1],
                'alt_accession' : answer[2],
                'tax_id' : answer[3],
                'pathways' : answer[4],
                'aliases' : answer[5],
                'topology' : answer[6]
            }
            return node_dict
        else:
            return None

    def get_node_by_alias(self,alias):
        """
        Returns the first node that has that alias. Only usable with unique nodes.
        :param alias: The alias of the node
        :type alias: str
        :return:
        """
        query = "SELECT * FROM node WHERE aliases LIKE '%" + alias + "%'"


        self.cursor.execute(query)
        self.db.commit()

        answer = self.cursor.fetchone()

        if answer:
            node_dict = {
                'id' : answer[0],
                'name' : answer[1],
                'alt_accession' : answer[2],
                'tax_id' : answer[3],
                'pathways' : answer[4],
                'aliases' : answer[5],
                'topology' : answer[6]
            }
            return node_dict
        else:
            return None

    def get_node_by_id(self,id):
        """
        This function returns a node dict, with the given id
        :param id: The answer id of the node that needs to be fetched.
        :type id: int
        :return: A node dict.
        """
        self.cursor.execute("SELECT * FROM node WHERE id = ? ",(id,))
        answer = self.cursor.fetchone()
        if not answer:
            return None
        else:
            if answer:
                id, name, alt_accession, tax_id, pathways, aliases, topology = answer
                node_dict = {
                    "id" : id,
                    "name" : name,
                    "alt_accession" : alt_accession,
                    "tax_id" : tax_id,
                    "pathways" : pathways,
                    "aliases" : aliases,
                    "topology" : topology
                }
            else:
                return None
            return  node_dict

    def update_node(self,node_dict):
        """
        This function updates a node in the :memory: db
        :param node_dict: A node dict that contains an id property. This node dict will be the base for modifying the db node entry
        :type node_dict: dict
        :return:
        """
        for k, v in node_dict.items():
            if v == None:
                node_dict[k] = '-'

        tup = (node_dict['name'],
               self.sort_attributes(node_dict['alt_accession']),
               node_dict['tax_id'],
               self.sort_attributes(node_dict['pathways']),
               self.sort_attributes(node_dict['aliases']),
               self.sort_attributes(node_dict['topology']),
               node_dict['id'])

        query = """
            UPDATE node
            SET name = ?, alt_accession = ?, tax_id = ?, pathways = ?, aliases = ?, topology = ?
            WHERE id = ?;
        """

        self.cursor.execute(query,tup)
        self.db.commit()

    def insert_edge(self, interactor_a_dict, interactor_b_dict, edge_dict):
        """
        This function inserts an edge to the edge table to the database
        :param interactor_a_dict: A node dict conaining the the node's data
        :type interactor_a_dict: dict
        :param interactor_b_dict: A node dict conaining the the node's data
        :type interactor_b_dict: dict
        :param edge_dict:
        :return:
        """
        query = """
                INSERT INTO `edge` (
                `interactor_a_node_id`,
                `interactor_b_node_id`,
                `interactor_a_node_name`,
                `interactor_b_node_name`,
                `interaction_detection_method`,
                `first_author`,
                `publication_ids`,
                `source_db`,
                `interaction_identifiers`,
                `interaction_types`,
                `confidence_scores`,
                `layer`
                )
                VALUES ( ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """

        tup = (
            interactor_a_dict['id'],
            interactor_b_dict['id'],
            interactor_a_dict['name'],
            interactor_b_dict['name'],
            self.sort_attributes(edge_dict['interaction_detection_method']),
            self.sort_attributes(edge_dict['first_author']),
            self.sort_attributes(edge_dict['publication_ids'], force_null_if_empty=False),
            self.sort_attributes(edge_dict['source_db']),
            self.sort_attributes(edge_dict['interaction_identifiers']),
            self.sort_attributes(edge_dict['interaction_types']),
            self.sort_attributes(edge_dict['confidence_scores']),
            edge_dict['layer']
        )

        self.db.execute(query, tup)
        self.db.commit()

    def save_db_to_file(self, db_file_name):
        """
        This function saves the memory db to name.db to a specific location, in only the filename was given the files will be saved to the same directory where the script is located
        :param db_file_name: The name of the source database in lowercase, eg.: biogrid
        :type db_file_name: string
        :return:
        """
        #creating the new db files and setting up the schema


        if '.db' not in db_file_name:
            export_file = db_file_name + '.db'
        else:
            export_file = db_file_name

        file_db = self.create_db(export_file)

        db_file_name = db_file_name.split('/')[-1]

        db_name = db_file_name.replace(".db", "")

        #attaching the empty database to the memory db
        tup = (export_file,db_name)
        self.db.execute("ATTACH ? as ?", tup)
        self.db.commit()

        #copying the contents of the memory.node table to the files.node table
        attached_table = (db_name +'.node')
        self.db.execute("INSERT INTO %s SELECT * FROM node" % attached_table)
        self.db.commit()

        #copying the content of the memory.edge table to the files.edge table
        attached_table = (db_name + '.edge')
        self.db.execute("INSERT INTO %s SELECT * FROM edge" % attached_table)
        self.db.commit()

        self.db.close()

        import slk3_db_validator
        if not slk3_db_validator.validate_db_file(export_file):
            print("ERROR! invalid DB file: " + export_file)
            sys.exit(1)