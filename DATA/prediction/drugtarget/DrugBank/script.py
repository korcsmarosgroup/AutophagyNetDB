"""
  Parses xml data from: https://www.drugbank.ca/releases/latest

  targetben van az interactant group
"""

# Imports
import xml.etree.ElementTree as ET
import sqlite3
import logging
import json
import io
# # Initiating logger
# logger = logging.getLogger()
# handler = logging.FileHandler('SLK3.log')
# logger.setLevel(logging.DEBUG)
# handler.setLevel(logging.DEBUG)
# formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
# handler.setFormatter(formatter)
# logger.addHandler(handler)

# Constants
XML_FILE = '../prediction/drugtarget/DrugBank/files/full_database.xml'
builder = 'SLK3_layers.db'
mapper = 'mapper.db'

#  Establishing connection
buildconn = sqlite3.connect(builder)
mapconn = sqlite3.connect(mapper)

# Setting up sqlite3 connection
# Read database to tempfile
conn2 = sqlite3.connect(mapper)
tempfile = io.StringIO()
for line in conn2.iterdump():
    tempfile.write('%s\n' % line)
conn2.close()
tempfile.seek(0)

# Create a database in memory and import from tempfile
DB = sqlite3.connect(":memory:")
with DB:
    DB.cursor().executescript(tempfile.read())
    DB.cursor().execute("CREATE INDEX mapp_uniprot ON MAPP(uniprot_ac);")
    DB.cursor().execute("CREATE INDEX uniprot_id ON UNIPROT_AC(id);")
    DB.cursor().execute("CREATE INDEX mapp_foreign ON MAPP(foreign_id);")
    DB.cursor().execute("CREATE INDEX uniprot_uniprot ON UNIPROT_AC(uniprot_ac);")
    DB.cursor().execute("CREATE INDEX gene_disp ON MAPP(gene_disp_name);")


def drugbank(logger):
    tree = ET.parse(XML_FILE)
    root = tree.getroot()

    drugdict = {}

    # Getting the name and id of the drug
    for drug in root.findall('.//{http://www.drugbank.ca}drug'):
        for name, firstid in zip(drug.findall('.//{http://www.drugbank.ca}name'),
                                 drug.findall('.//{http://www.drugbank.ca}drugbank-id')):
            drugname = name.text
            #  Getting DrugBank id in psi-mi format
            if 'primary' in firstid.attrib and firstid.attrib['primary'] == 'true':     # TODO: Do we need non primary ids?
                drugid = 'DrugBank:' + firstid.text
                # Getting name of the drug
                drugdata = ''.join([drugid, '(', drugname, ')'])
                #  Adding it to a dictionary's keys
                if drugdata not in drugdict.keys():
                    drugdict[drugdata] = ''
                # Getting names of genes targeted by the drug
                for targets in drug.findall('.//{http://www.drugbank.ca}targets'):
                    for target in targets.findall('.//{http://www.drugbank.ca}target'):
                        for idents in target.findall('.//{http://www.drugbank.ca}external-identifiers'):
                            for ident in idents.findall('.//{http://www.drugbank.ca}external-identifier'):
                                for resource in ident.findall('.//{http://www.drugbank.ca}resource'):
                                    if resource.text == 'UniProt Accession':    # TODO: do we need other databases?
                                        # Adding gene names to the dictionary
                                        for ids in ident.findall('.//{http://www.drugbank.ca}identifier'):
                                            drugdict[drugdata] = ids.text

    # Mapping gene names to uniprot ids
    with DB:
        c = DB.cursor()
        for key, value in drugdict.items():
            c.execute("SELECT UNIPROT_AC.uniprot_ac FROM MAPP "
                      "LEFT JOIN UNIPROT_AC ON MAPP.uniprot_ac=UNIPROT_AC.id "
                      "WHERE MAPP.gene_disp_name == '%s' GROUP BY MAPP.foreign_id"
                      % value)
            while True:
                row = c.fetchone()
                if row is None:
                    break
                else:
                    uniprot = row[0]
                    drugdict[key] = 'Uniprot:' + uniprot

    # Adding drugtarget data to merger
    with buildconn:
        c = buildconn.cursor()
        for key, value in drugdict.items():
            c.execute("UPDATE node SET topology = node.topology || ? WHERE node.name = ?"
                      , ('|drugtarget:' + key, value))


if __name__ == '__main__':
    drugbank(logger)
