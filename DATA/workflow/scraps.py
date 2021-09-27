
import sqlite3
conn = sqlite3.connect('../workflow/SLK3_layers.db')
c = conn.cursor()
c.execute("SELECT count(*) FROM layer5")
num = c.fetchone()
print(num)
