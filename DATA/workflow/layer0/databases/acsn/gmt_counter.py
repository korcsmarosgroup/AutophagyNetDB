f1 = open('acsn_master_curated.gmt')
f2 = open('acsn_ppi.sif.txt')

bioinformatics_shit_set = set()

for line in f1:
  line = line.strip()
  cells = line.split('\t')
  if len(cells) > 4:
    bioinformatics_shit_set.add(cells[0])

ppi_counter = 0
shitcounter = 0

for line in f2:
  ppi_counter += 1
  cells = line.strip().split('\t')
  if cells[0] in bioinformatics_shit_set or cells[2] in bioinformatics_shit_set:
    shitcounter += 1

print(ppi_counter)
print(shitcounter)