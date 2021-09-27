import sqlite3

def uniprot_to_foreign_id(my_uniprot):
    conn = sqlite3.connect('test.db')
    with conn:
        c = conn.cursor()
        c.execute("SELECT foreign_id FROM uniprot_ac JOIN mapp ON mapp.uniprot_ac=uniprot_ac.id WHERE uniprot_ac.uniprot_ac='%s'" %my_uniprot)
        allrows = c.fetchall()
        return allrows

#print(uniprot_to_foreign_id('P04637'))

def foreign_id_to_uniprot(my_foreign_id):
    conn = sqlite3.connect('test.db')
    with conn:
        c = conn.cursor()
        c.execute("SELECT uniprot_ac.uniprot_ac FROM uniprot_ac JOIN mapp ON mapp.uniprot_ac=uniprot_ac.id WHERE mapp.foreign_id='%s'" %my_foreign_id)
        allrows = c.fetchall()
        return allrows

print(foreign_id_to_uniprot('Q15086'))


def searching_in_db(expression):
    db = sqlite3.connect("test.db")
    searched_value_list = []
    with db:
        c = db.cursor()
        c.execute("SELECT uniprot_ac FROM data WHERE external_id LIKE ? ", ('%' + expression + '%', ))
        for row in c.fetchall():
            searched_value_list.append(row[0])
    return searched_value_list

print(searching_in_db('NP_000'))



