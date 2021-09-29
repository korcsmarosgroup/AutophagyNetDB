__author__ = 'blaise'

from networkx import DiGraph

g = DiGraph()

g.add_node('A1')
g.add_node('A2')
g.add_node('A3')
g.add_edge('A1','A3',{ 'shit' : 'lol'})

g.remove_node('A1')


print(g.edges())