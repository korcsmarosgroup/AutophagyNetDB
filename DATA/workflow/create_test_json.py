import json


data = []
with open('edges_RC3_07_16.json') as json_file:
    json_object = json.load(json_file)
    for edge in json_object:
        if edge['source'] == 'P45897' or edge['target'] == 'P45897':
            data.append(edge)

with open('slk3_test_edges.json', 'w') as outfile:
    json.dump(data, outfile)