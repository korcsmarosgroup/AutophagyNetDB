__author__ = 'blaise'

import json
import pprint

pathways_map = {}

with open("pathway_map.tsv") as pathways_file:
    for line in pathways_file:
        reactome_pathway, signalink_pathway = line.split('\t')
        pathways_map[reactome_pathway] = {
            "signalink_pathway" : signalink_pathway,
            "sum" : 0
            }

with open("homo_sapiens.mitab.interactions.txt") as mitab_file:
    #skipping header
    mitab_file.readline()

    for line in mitab_file:
        cells = line.strip().split('\t')

        if cells[9] == "taxid:9606" and cells[10] == "taxid:9606":
            #source_uniprot = cells[0]
            #target_uniprot = cells[1]
            reactome_pathway_id = cells[15].replace('reactome:','')

            if pathways_map.has_key(reactome_pathway_id):
                pathways_map[reactome_pathway_id]["sum"] += 1

result = json.dumps(pathways_map)

sum_all = 0

with open("mitab_results", "w") as mitab_result_file:
    for key, value in pathways_map.iteritems():
        if pathways_map[key]["sum"] != 0:
            mitab_result_file.write("%s %i \n" % (key, pathways_map[key]["sum"]) )
            sum_all += pathways_map[key]["sum"]
    mitab_result_file.write("Sum %i edges" % (sum_all))