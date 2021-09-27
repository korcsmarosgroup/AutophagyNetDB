__author__ = 'blaise'

egf = open('SIGNOR_EGF.tsv')
shh = open('SIGNOR_Hedgehog.tsv')

egfs = set()
shhs = set()

egf.readline()
for line in egf:
    cells = line.split('\t')
    egfs.add(cells[2])
    egfs.add(cells[6])

for line in shh:
    cells = line.split('\t')
    shhs.add(cells[2])
    shhs.add(cells[6])

print('egf: '+ str(egfs))
print('shh: '+ str(shhs))
print('intersection: '+ str(egfs.intersection(shhs)))