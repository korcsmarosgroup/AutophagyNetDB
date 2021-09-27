import sqlite3

def mappingDB_structure(mappingDB):
    conn = sqlite3.connect(mappingDB)
    with conn:
        c = conn.cursor()
        c.execute('''DROP TABLE IF EXISTS data''')
        c.execute('''CREATE TABLE data (id INTEGER PRIMARY KEY AUTOINCREMENT,
                                        external_id CHAR(100) NOT NULL,
                                        external_id_type CHAR(100) NOT NULL,
                                        uniprot_ac CHAR(10) NOT NULL,
                                        is_reviewed INT DEFAULT 0,
                                        is_primary INT DEFAULT 0,
                                        length INT NOT NULL)''')

def generate_mock_data(mappingDB):
    conn = sqlite3.connect(mappingDB)

    with conn:
        c = conn.cursor()
        insert_lines = []
        with open('mock_db.csv') as input_table:
            for line in input_table:
                cells = line.strip().split(',')
                insert_lines.append((cells[0],cells[1],cells[2],int(cells[3]), int(cells[4]), int(cells[5])))
        c.executemany('INSERT INTO data(external_id, external_id_type, uniprot_ac, is_reviewed, is_primary, length) '
                      'VALUES(?,?,?,?,?,?)', insert_lines)

if __name__ =='__main__':
    mappingDB_structure('test.db')
    #generate_mock_data('test.db')

