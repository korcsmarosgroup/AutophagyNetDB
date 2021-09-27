import sqlite3
import xml.etree.ElementTree as ET
import requests
import gzip
import os


class CreateMappingDB():
    def __init__(self, mappingDBfile = '~/.mappingDB', debug=False):
        self.mappingDBfile = mappingDBfile
        self.debug = debug
        if self.debug:
            try:
                os.remove(self.mappingDBfile)
            except FileNotFoundError:
                pass
        self.conn = sqlite3.connect(self.mappingDBfile)
        self.db_blacklist = ['PubMed', 'DOI', 'PDBsum', 'PDB', 'EMBL', 'CCDS', 'GO', 'InterPro', 'PROSITE', 'Pfam', 'MIM', 'CCDS']
        self.db_whitelist = ['ELM', 'BioGrid', 'Ensembl', 'EnsemblMetazoa', 'FlyBase',
                          'GenBank', 'GeneCards', 'GeneID', 'HGNC', 'IntAct', 'MINT',
                          'Reactome', 'SIGNOR', 'SignaLink', 'WormBase', 'RefSeq']

        #
        self.current_species = ""
        self.DB_species = {}
        self.DB_foreignIDs = {}

        with self.conn:
            cur = self.conn.cursor()
            cur.execute("SELECT name FROM sqlite_master WHERE type='table';")
            tables = [i[0] for i in cur.fetchall()]
        if 'data' not in tables:
            self.mappingDB_structure()

    def add_species(self, taxID):
        self.new_species(taxID)
        self.current_species = taxID
        uniprot_xml_file = self.uni_to_xml(taxID)
        with gzip.open(uniprot_xml_file) as input_file:
            self.process_uniprot_xml(input_file)
        os.remove(uniprot_xml_file)

    def mappingDB_structure(self):
        with self.conn:
            c = self.conn.cursor()
            c.execute('''DROP TABLE IF EXISTS species''')
            c.execute('''DROP TABLE IF EXISTS uniprot_ac''')
            c.execute('''DROP TABLE IF EXISTS id_type''')
            c.execute('''DROP TABLE IF EXISTS mapp''')
            c.execute('''CREATE TABLE species (id INTEGER PRIMARY KEY AUTOINCREMENT,
                                               tax_id CHAR(25) NOT NULL,
                                               tax_name CHAR(250) NOT NULL)''')
            c.execute('''CREATE TABLE uniprot_ac (id INTEGER PRIMARY KEY AUTOINCREMENT,
                                                  uniprot_ac CHAR(10) NOT NULL,
                                                  uniprot_ac_alt_acc CHAR(10) NOT NULL,
                                                  taxon INTEGER NOT NULL,
                                                  length INT NOT NULL,
                                                  is_reviewed INT NOT NULL)''')
            c.execute('''CREATE TABLE id_type (id INTEGER PRIMARY KEY AUTOINCREMENT,
                                               foreign_id_type CHAR(250) NOT NULL)''')
            c.execute('''CREATE TABLE mapp (id INTEGER PRIMARY KEY AUTOINCREMENT,
                                            foreign_id CHAR(100),
                                            foreign_id_type INTEGER NOT NULL,
                                            uniprot_ac INTEGER NOT NULL,
                                            gene_name CHAR(100),
                                            gene_name_synonym CHAR(100),
                                            gene_disp_name CHAR(100),
                                            prot_full_name CHAR(10000))''')

    def new_species(self, taxID, tax_name=""):
        with self.conn:
            DBcur = self.conn.cursor()
            DBcur.execute(
                'INSERT INTO species(tax_id, tax_name) VALUES(?,?)',
                (taxID, tax_name))
            self.DB_species[taxID] = DBcur.lastrowid

    def new_foreignID_type(self, foreign_id_type):
        with self.conn:
            DBcur = self.conn.cursor()
            DBcur.execute(
                'INSERT INTO id_type(foreign_id_type) VALUES(?)',
                (foreign_id_type,))
            rowid = DBcur.lastrowid
            self.DB_foreignIDs[foreign_id_type] = rowid
        return rowid

    def get_foreignID_type(self, foreign_id_type):
        try:
            return self.DB_foreignIDs[foreign_id_type]
        except KeyError:
            return self.new_foreignID_type(foreign_id_type)

    def extract_entry_data(self, entry, db_whitelist):
        entry_data = {}
        # Getting uniprot accessions
        acs = [i.text for i in entry.findall('.//{http://uniprot.org/uniprot}accession')]
        entry_data['accession'] = acs
        # Getting ids from 3rd party databases

        xrefs = []
        for i in entry.findall('.//{http://uniprot.org/uniprot}dbReference'):
            if i.attrib['type'] not in db_whitelist:
                continue
            xrefs.append((i.attrib['type'], i.attrib['id']))
            for ensgene in i.findall('{http://uniprot.org/uniprot}property'):
                xrefs.append((i.attrib['type'], ensgene.attrib['value'].upper()))
        # Also add uniprot ids to foreign ids

        if len(entry_data['accession']) > 1:
            for alt_uniprot in entry_data['accession'][1:]:
                xrefs.append(('alt_uniprot', alt_uniprot))
        else:
            xrefs.append(('Uniprot', entry_data['accession'][0]))

        entry_data['xref'] = xrefs
        # Getting dataset
        entry_data['dataset'] = entry.attrib['dataset']
        # Getting sequence length
        entry_data['length'] = entry.findall('{http://uniprot.org/uniprot}sequence')[0].attrib['length']
        # Getting primary gene name
        entry_data['gene_name'] = [i.text for name in entry.findall('{http://uniprot.org/uniprot}gene')
                                   for i in name.findall('{http://uniprot.org/uniprot}name') if i.attrib['type']=='primary']
        # Getting secondary gene name
        entry_data['gene_name_syn'] = [i.text for name in entry.findall('{http://uniprot.org/uniprot}gene')
                                       for i in name.findall('{http://uniprot.org/uniprot}name') if i.attrib['type']=='synonym']
        entry_data['gene_name_ordloc'] = [i.text for name in entry.findall('{http://uniprot.org/uniprot}gene')
                                          for i in name.findall('{http://uniprot.org/uniprot}name') if i.attrib['type']=='ordered locus']
        # Getting gene orfs
        entry_data['gene_name_orf'] = [i.text for name in entry.findall('{http://uniprot.org/uniprot}gene')
                                       for i in name.findall('{http://uniprot.org/uniprot}name') if i.attrib['type']=='ORF']
        # Getting protein names
        prot_full = [i.text for name in entry.findall('.//{http://uniprot.org/uniprot}recommendedName')
                     for i in name.findall('{http://uniprot.org/uniprot}fullName')]
        entry_data['prot_full_name'] = prot_full
        # Getting displayed name
        entry_data['gene_dispname'] = [i.text for i in entry.findall('.//{http://uniprot.org/uniprot}name')]
        return entry_data

    def insert_entry_to_sqlite(self, entry_data, DBcur):
        insert_rows = []
        # Check if uniprot id is reviewed
        is_reviewed = 1 if entry_data['dataset'] == 'Swiss-Prot' else 0
        if is_reviewed == 1:
            primary_ac = entry_data['accession'][0]
            if len(entry_data['accession']) > 1:
                alt_acc = '|'.join(entry_data['accession'][1:])
            else:
                alt_acc = '-'
            DBcur.execute('INSERT INTO uniprot_ac(uniprot_ac, uniprot_ac_alt_acc, taxon, length, is_reviewed) VALUES (?,?,?,?,?)',
                          (primary_ac, alt_acc, self.DB_species[self.current_species], entry_data['length'], is_reviewed))
            primary_ac_id = DBcur.lastrowid
            uniprot_fid = self.get_foreignID_type("Uniprot")

            if len(entry_data['gene_name_syn']) != 0:
                syn = ' | '.join(entry_data['gene_name_syn'])
            else:
                syn = None
            if len(entry_data['prot_full_name']) != 0:
                prot_name = entry_data['prot_full_name'][0]
            else:
                prot_name = None
            if len(entry_data['gene_dispname']) != 0:
                dispname = entry_data['gene_dispname'][0]
            else:
                dispname = None
            for gene in entry_data['gene_name']:
                for ac in entry_data['accession'][1:]:
                    insert_rows.append((ac, uniprot_fid, primary_ac_id, gene, syn, dispname, prot_name))
                for xref in entry_data['xref']:
                    fid = self.get_foreignID_type(xref[0])
                    insert_rows.append((xref[1], fid, primary_ac_id, gene, syn, dispname, prot_name))

            DBcur.executemany('INSERT INTO mapp(foreign_id, foreign_id_type, uniprot_ac, gene_name, gene_name_synonym,'
                              'gene_disp_name, prot_full_name) '
                              'VALUES (?,?,?,?,?,?,?)',
                              insert_rows)

    def uni_to_xml(self, taxID):
        local_filename = taxID + '.xml.gz'
        query_term = 'organism:%s' % str(taxID)
        download_params = {
            'sort':'score',
            'compress':'yes',
            'query':query_term,
            'format':'xml',
            'force':'yes'
        }
        if self.debug:
            download_params['limit'] = '1000'

        r = requests.get(
            "http://www.uniprot.org/uniprot/",
            stream=True, params=download_params)

        with open(local_filename, 'wb') as f:
            for chunk in r.iter_content(chunk_size=1024):
                if chunk:
                    f.write(chunk)

        return local_filename

    def process_uniprot_xml(self, uniprot_xml):
        with self.conn:
            DBcur = self.conn.cursor()
            n = 1
            #for entry in xml.findall('{http://uniprot.org/uniprot}entry'):
            for event, elemnt in ET.iterparse(uniprot_xml, events=('end',)):
                if elemnt.tag == '{http://uniprot.org/uniprot}entry':
                    entry_data = self.extract_entry_data(elemnt, self.db_whitelist)
                    self.insert_entry_to_sqlite(entry_data, DBcur)
                    #print(entry_data)
                    elemnt.clear()
                    if self.debug:
                        print(n)
                    n += 1


if __name__ == '__main__':
    MDB = CreateMappingDB(mappingDBfile='../../../DATA/workflow/casesense_mapper.db', debug=False)
    MDB.add_species('9606')
    MDB.add_species('7227')
    MDB.add_species('6239')
    MDB.add_species('7955')
