
import sqlite3
conn = sqlite3.connect('../workflow/SLK3_layers.db')
c = conn.cursor()
c.execute("SELECT count(*) FROM miRNA")
num = c.fetchone()
print(num)
