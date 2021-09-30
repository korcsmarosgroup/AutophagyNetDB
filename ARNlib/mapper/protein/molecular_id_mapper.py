"""
 Maps all nodes of a parsed database to a more universal one, based on the node's molecular type type, by the previously created mapping data sets for each type.
"""

# Impports
import sqlite3
import io

from ARNlib.SQLiteDBApi.sqlite_db_api import PsimiSQL


class MolecularIDMapper:
    def __init__(self, db, layer, PROT_DBname, LNCRNAMAP_DBname):
        """
        :param db: name of the parsed source database
        :param PROT_DBname: if value is not None, database is created in memory
        :argument DICTIONARY_DB_LOCATION: location of the mapping db, output of create_mapping_db
        :argument SQL_SEED_LOCATION
        :argument SOURCE_DB_LOCATION: location of the parsed source database
        :argument DESTINATION_DB_LOCATION: location where the mapped db will be saved
        """

        # Declaring, and assigning constants
        self.DICTIONARY_DB_LOCATION = PROT_DBname
        self.SQL_SEED_LOCATION = '../../ARNlib/SQLiteDBApi/network-db-seed.sql'
        self.source_db_map = {
            "acsn": 'ACSN',
            "innatedb": 'InnateDB',
            "reactome": 'Reactome',
            "signor": 'Signor',
            "slk2": 'SLKv2.0',
            "slk3": "SLKv3.0",
            "slk21": "SLKv2.1",
            "SLK2_endo": 'SLKv2.0',
            "SLK2_scaffold": 'SLKv2.0',
            "tcr": "TCRcuration",
            "signafish": "SignaFish",
            "biogrid": 'TheBiogrid',
            "starbase_mirna": 'StarBase',
            "TFBS": "PSSMprediction",
            "TFBS_Orsi": "TFlink",
            "starbase_lncrna": 'StarBase'
        }
        if db in self.source_db_map.keys():
            self.SOURCE_DB_TYPE = self.source_db_map[db]
        else:
            self.SOURCE_DB_TYPE = db
        self.layer = layer
        # The db we want to map
        self.SOURCE_DB_LOCATION = 'all_output/' + db + '.db'
        # Saving location
        self.DESTINATION_DB_LOCATION = '../../ARNlib/mapper/protein/output/' + db + '_mapped.db'
        # Protein map db
        self.DICTIONARY_DB = sqlite3.connect(self.DICTIONARY_DB_LOCATION)
        self.DICTIONARY_DB_CURSOR = self.DICTIONARY_DB.cursor()
        # lncRNA map db
        if self.layer == 'layer7' or self.layer == 'layer5':
            self.LNCRNAMAP_DB_LOCATION = LNCRNAMAP_DBname
            self.LNCRNAMAP_DB = sqlite3.connect(self.LNCRNAMAP_DB_LOCATION)

        self.PROT_DBname = PROT_DBname
        if self.PROT_DBname is not None:
            # Read database to tempfile
            self.con = sqlite3.connect(self.PROT_DBname)
            tempfile = io.StringIO()
            for line in self.con.iterdump():
                tempfile.write('%s\n' % line)
            self.con.close()
            tempfile.seek(0)

            # Create a database in memory and import from tempfile
            self.PROT_DB = sqlite3.connect(":memory:")
            with self.PROT_DB:
                self.PROT_DB.cursor().executescript(tempfile.read())
                self.PROT_DB.cursor().execute("CREATE INDEX map_uniprot ON MAPP(uniprot_ac);")
                self.PROT_DB.cursor().execute("CREATE INDEX uniprotac_id ON UNIPROT_AC(id);")
                self.PROT_DB.cursor().execute("CREATE INDEX taxid ON SPECIES(tax_id);")
                self.PROT_DB.cursor().execute("CREATE INDEX map_foreign ON MAPP(foreign_id);")
                self.PROT_DB.cursor().execute("CREATE INDEX is_reviewed ON UNIPROT_AC(id);")
        else:
            self.PROT_DB = sqlite3.connect(self.DICTIONARY_DB_LOCATION)
            self.PROT_DB.cursor().execute("CREATE INDEX index_name ON mapp (foreign_id);")

        # For lncRNA and miRNA
        if self.layer == 'layer7' or self.layer == 'layer5':
            self.LNCRNAMAP_DBname = LNCRNAMAP_DBname
            if self.LNCRNAMAP_DBname is not None:
                # Read database to tempfile
                self.con = sqlite3.connect(self.LNCRNAMAP_DBname)
                tempfile = io.StringIO()
                for line in self.con.iterdump():
                    tempfile.write('%s\n' % line)
                self.con.close()
                tempfile.seek(0)

                # Create a database in memory and import from tempfile
                self.LNCRNAMAP_DB = sqlite3.connect(":memory:")
                with self.LNCRNAMAP_DB:
                    self.LNCRNAMAP_DB.cursor().executescript(tempfile.read())
                    self.LNCRNAMAP_DB.cursor().execute("CREATE INDEX index_name ON mapper (orig_ac);")
            else:
                self.LNCRNAMAP_DB = sqlite3.connect(self.LNCRNAMAP_DB_LOCATION)

        self.new_db = PsimiSQL(self.SQL_SEED_LOCATION)
        # iterating through the old_db's nodes
        self.source_db = sqlite3.connect(self.SOURCE_DB_LOCATION)
        self.source_db.row_factory = sqlite3.Row
        self.cur = self.source_db.cursor()

    def add_node(self, old_node_id, old_to_new_node_ids_dict, new_name,
                 new_alt, new_taxid, new_pathways, new_alias, new_topo, new_db_api):
        """
        :param old_node_id: node id from the source db's node table
        :param old_to_new_node_ids_dict: A dictionary that contains an old node id as key and new node ids as values
        :param new_name: mapped uniprot ac of the mapped node
        :param new_alt: foreign id of the mapped node
        :param new_taxid: taxid
        :param new_pathways: pathway
        :param new_alias: aliases
        :param new_topo: topology
        :param new_db_api: A PsimiSQL object
        """

        new_node_dict = {
            "name": new_name,
            "alt_accession": new_alt,
            "tax_id": new_taxid,
            "pathways": new_pathways,
            "aliases": new_alias,
            "topology": new_topo
        }

        # inserting the node to the PSI-MI SQLite
        new_db_api.insert_unique_node(new_node_dict)

        # getting the new last row id of the inserted node
        new_node_id = new_node_dict['id']

        # if the node maps to more than one swissprot uniprot id it will be inserted for every swissprot id and
        # this function will be called for every insertion
        if old_node_id not in old_to_new_node_ids_dict:
            old_to_new_node_ids_dict[old_node_id] = [new_node_id]
        else:
            old_to_new_node_ids_dict[old_node_id].append(new_node_id)

    def main(self):
        old_node_ids_dict = {}
        invalid_node_counter = 0

        # MAPPING NODES
        self.cur.execute("SELECT * FROM node")
        node_counter = 0
        while True:
            # Getting data for each node
            node_row = self.cur.fetchone()
            node_counter += 1
            # Until the last row
            if node_row is None:
                break

            # Getting the old information into a dictionary
            old_node_dict = dict(node_row)

            # For all other databases
            foreign_id = old_node_dict['name'].split(':')[1].strip()
            # Taxid
            taxid = old_node_dict['tax_id'].split(':')[1].split('(')[0]

            # miRNA and lncRNA mapping
            if self.layer == 'layer7' or self.layer == 'layer5':
                with self.LNCRNAMAP_DB:
                    c = self.LNCRNAMAP_DB.cursor()
                    for indiv_id in foreign_id.split(','):
                        indiv_id = indiv_id.replace('"', '').lower()
                        c.execute(
                            '''SELECT mapped_ac FROM MAPPER WHERE '%s' = MAPPER.orig_ac GROUP BY MAPPER.orig_ac'''
                            % indiv_id
                        )
                        while True:
                            allrows = c.fetchone()
                            if allrows is None:
                                break
                            else:
                                m.add_node(node_row['id'], old_node_ids_dict, 'RNACentral:' + allrows[0], node_row['name'],
                                           node_row['tax_id'], node_row['pathways'], node_row['aliases'], node_row['topology'],
                                           self.new_db)
                with self.PROT_DB:
                    c2 = self.PROT_DB.cursor()
                    c2.execute(
                        "SELECT UNIPROT_AC.uniprot_ac, UNIPROT_AC.uniprot_ac_alt_acc FROM UNIPROT_AC "
                        "JOIN MAPP ON MAPP.uniprot_ac=UNIPROT_AC.id "
                        "JOIN SPECIES ON SPECIES.id=UNIPROT_AC.taxon WHERE SPECIES.tax_id='%s'"
                        "AND MAPP.foreign_id='%s' AND UNIPROT_AC.is_reviewed = '1' GROUP BY MAPP.foreign_id"
                        % (taxid, foreign_id)
                    )
                    while True:
                        allrows = c2.fetchone()
                        if allrows is None:
                            break
                        else:
                            if node_row['alt_accession'] is not None \
                                 and node_row['name'] is not None and allrows[1] is not None and node_row['alt_accession'] != '-':
                                m.add_node(node_row['id'], old_node_ids_dict, 'Uniprot:' + allrows[0], '|'.join([node_row['name'], node_row['alt_accession'], "Uniprot:" + allrows[1]]),
                                           node_row['tax_id'], node_row['pathways'], node_row['aliases'], node_row['topology'],
                                           self.new_db)
                            else:
                                m.add_node(node_row['id'], old_node_ids_dict, 'Uniprot:' + allrows[0], node_row['name'],
                                           node_row['tax_id'], node_row['pathways'], node_row['aliases'], node_row['topology'],
                                           self.new_db)

            # Protein mapping
            else:
                with self.PROT_DB:
                    c = self.PROT_DB.cursor()
                    # Getting uniprot acs for each node and adding the node with new data to the new database
                    c.execute(
                        "SELECT UNIPROT_AC.uniprot_ac, UNIPROT_AC.uniprot_ac_alt_acc FROM UNIPROT_AC "
                        "JOIN MAPP ON MAPP.uniprot_ac=UNIPROT_AC.id "
                        "JOIN SPECIES ON SPECIES.id=UNIPROT_AC.taxon WHERE SPECIES.tax_id='%s'"
                        "AND MAPP.foreign_id='%s' AND UNIPROT_AC.is_reviewed = '1' GROUP BY MAPP.foreign_id"
                        % (taxid, foreign_id)
                    )
                    while True:
                        allrows = c.fetchone()
                        if allrows is None:
                            break
                        else:
                            if node_row['alt_accession'] is not None \
                                    and node_row['name'] is not None and allrows[1] is not None and node_row['alt_accession'] != '-':
                                m.add_node(node_row['id'], old_node_ids_dict, 'Uniprot:' + allrows[0], '|'.join([node_row['name'], node_row['alt_accession'], "Uniprot:" + allrows[1]]),
                                           node_row['tax_id'], node_row['pathways'], node_row['aliases'], node_row['topology'],
                                           self.new_db)
                            else:
                                m.add_node(node_row['id'], old_node_ids_dict, 'Uniprot:' + allrows[0], node_row['name'],
                                           node_row['tax_id'], node_row['pathways'], node_row['aliases'], node_row['topology'],
                                           self.new_db)

        # MAPPING EDGES

        self.cur.execute("SELECT * from EDGE")
        while True:
            edge_row = self.cur.fetchone()
            if edge_row is None:
                break
            else:
                # Since we get the old interactor id's from this query we can simply look up ther new id(s) in the old_node_ids dict
                # if both nodes mapped we add them as an edge to the new db
                # deduplicating
                for key, value in old_node_ids_dict.items():
                    old_node_ids_dict[key] = list(set(value))
                if edge_row['interactor_a_node_id'] in old_node_ids_dict and edge_row['interactor_b_node_id'] in old_node_ids_dict:
                    # looping through every new 'A' node
                    for new_node_id_a, new_node_id_b in zip(old_node_ids_dict[edge_row['interactor_a_node_id']], old_node_ids_dict[edge_row['interactor_b_node_id']]):
                        new_node_a_dict = self.new_db.get_node_by_id(new_node_id_a)
                        new_node_b_dict = self.new_db.get_node_by_id(new_node_id_b)

                        # looping through every new 'B' node for every new 'A' node and inserting them as an edge
                        #for new_node_id_b in old_node_ids_dict[edge_row[2]]:
                            #new_node_b_dict = self.new_db.get_node_by_id(new_node_id_b)
                        new_edge_dict = dict(edge_row)
                        new_edge_dict['interactor_a_node_id'] = new_node_id_a
                        new_edge_dict['interactor_b_node_id'] = new_node_id_b
                        new_edge_dict['source_db'] = self.SOURCE_DB_TYPE

                        # inserting the new node
                        self.new_db.insert_edge(new_node_a_dict, new_node_b_dict, new_edge_dict)
                else:
                    invalid_node_counter += 1

        # Saving the mapped database
        self.new_db.save_db_to_file(self.DESTINATION_DB_LOCATION)
        print(invalid_node_counter)


# if __name__ == '__main__':
#     for db in DB_LIST:
#         m = MolecularIDMapper(db, layer=None)
#         m.main()
