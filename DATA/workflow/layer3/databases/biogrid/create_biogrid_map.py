"""
This scriot generates an anydbm .db files from the biogrid map (that was downloaded from biogrid)
"""
import dbm

db = dbm.open('./biogrid_mapper.db', 'c')

with open('BIOGRID-SYSTEM-3.4.145.mitab/BIOGRID-SYSTEM-Far_Western-3.4.145.mitab.txt', 'r') as f:
    while('BIOGRID_ID\tIDENTIFIER_VALUE' not in f.readline()):
        pass

    for line in f:
        biogrid_id, id_value, id_type, organism_name = line.strip('\n').split('\t')

        # logic here
        if organism_name not in ('Homo sapiens', 'Drosophila melanogaster', 'Caenorhabditis elegans'):
            continue

        if id_type != "SWISS-PROT":
            continue

        if biogrid_id not in db:
            db[biogrid_id] = id_value
        else:
            db[biogrid_id] += '|' + id_value

db.close()

#counting the keys (biogrid accessions) that have more than one value (Uniprot SWISSPROT accession)
fdb = anydbm.open('biogrid_mapper.db', 'r')
sum = 0
all_sum = 0
for k,v in fdb.iteritems():
    all_sum += 1
    if '|' in v:
        sum += 1
print(sum)
print(all_sum)