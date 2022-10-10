"""
    Select interaction between core proteins from 'atg reg' layer databases
"""
import sqlite3

core_mapped_nodes = '../../../ARNlib/mapper/protein/output/ARN1Core_mapped.db'

# TODO iterate
# reg_mapped_db = '../../../ARNlib/mapper/protein/output/HumanAutophagyDB_mapped.db'

reg_mapped_db = ['../../../ARNlib/mapper/protein/output/HumanAutophagyDB_mapped.db',
                 '../../../ARNlib/mapper/protein/output/chip_behrends_mapped.db',
                 '../../../ARNlib/mapper/protein/output/manual_curation_mapped.db',
                 '../../../ARNlib/mapper/protein/output/OmniPath_mapped.db',
                 '../../../ARNlib/mapper/protein/output/IntAct_mapped.db',
                 '../../../ARNlib/mapper/protein/output/biogrid_mapped.db',
                 '../../../ARNlib/mapper/protein/output/HPRD_mapped.db',
                 '../../../ARNlib/mapper/protein/output/ComPPI_mapped.db']


core_uniprot = []

core_sqliteConnection = sqlite3.connect(core_mapped_nodes)
with core_sqliteConnection:
    core_cursor = core_sqliteConnection.cursor()

    for row in core_cursor.execute("SELECT name FROM node"):
        core_uniprot.append(row[0])

print(core_uniprot)
print(','.join(['?']*len(core_uniprot)))
for db in reg_mapped_db:
    reg_conn = sqlite3.connect(db)
    with reg_conn:
        reg_cur = reg_conn.cursor()
        reg_cur.execute("UPDATE edge SET layer='0' WHERE interactor_a_node_name IN ({seq})"
                        " AND interactor_b_node_name IN ({seq})".format(seq=','.join(['?']*len(core_uniprot))), [*core_uniprot, *core_uniprot])
        reg_cur.execute("SELECT changes()")
        print(reg_cur.fetchall())
