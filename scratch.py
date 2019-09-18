from cl1_randomized import *
from common import *
from graph import *

# get_quality("cl1_gold_standard/gold_standard/mips_3_100.txt", "original_algo_findings.txt")
name = loadData('pickles/most_recent')['title']
print(name)
x = loadData(name)
print(x.graph.hash_graph[980])
print(x.graph.hash_graph[981])
print(x.graph.hash_graph[104])
print(x.graph.hash_graph[105])
print(x.graph.hash_graph[99])
print(x.graph.hash_graph[117])
