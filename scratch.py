from cl1_randomized import *
from common import loadData
from graph import *
name = \
"pickles/pickle+gavin2006_socioaffinities_rescaled+2019-08-29_01:37:29:650865"




x = loadData(name)
print(x.quality_report)
print(x.initial_clustering)
print(x.final_clusters_stats)
print(x.found)
print(x.gsc_appearing_found_stats)
print(x.gsc_appearing_notFound_stats)
print(x.gsc_appearing_stats)
for d in [x.final_clusters_stats, x.gsc_appearing_stats, x.gsc_appearing_found_stats, x.gsc_appearing_notFound_stats]:
    print(d['average_cohesiveness'], d['average_density'], d['average_size'])
    print()

import matplotlib.pyplot as plt
import networkx as nx

G = nx.Graph()

G.add_edge('a', 'b', weight=0.6)
G.add_edge('a', 'c', weight=0.2)
G.add_edge('c', 'd', weight=0.1)
G.add_edge('c', 'e', weight=0.7)
G.add_edge('c', 'f', weight=0.9)
G.add_edge('a', 'd', weight=0.3)

elarge = [(u, v) for (u, v, d) in G.edges(data=True) if d['weight'] > 0.5]
esmall = [(u, v) for (u, v, d) in G.edges(data=True) if d['weight'] <= 0.5]

pos = nx.spring_layout(G)  # positions for all nodes

# nodes
nx.draw_networkx_nodes(G, pos, node_size=700)

# edges
nx.draw_networkx_edges(G, pos, edgelist=elarge,
                       width=6)
nx.draw_networkx_edges(G, pos, edgelist=esmall,
                       width=6, alpha=0.5, edge_color='b', style='dashed')

# labels
nx.draw_networkx_labels(G, pos, font_size=20, font_family='sans-serif')

plt.axis('off')
plt.show()