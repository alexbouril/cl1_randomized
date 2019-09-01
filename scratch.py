from cl1_randomized import *
from common import *
from graph import *
n = "gavin2006_socioaffinities_rescaled+2019-08-31_17:14:11:720422"
name = \
"pickles/pickle+"+n


x = loadData(name)
print(stringify_construction_log(x.construction_log))
# exit(0)

x = loadData(name)
my_data = loadData(name)
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
# pp.pprint(x.graph.hash_graph)

import matplotlib.pyplot as plt
import plotly.graph_objects as go
import networkx as nx

# G.add_edge('a', 'b', weight=0.6)
# G.add_edge('a', 'c', weight=0.2)
# G.add_edge('c', 'd', weight=0.1)
# G.add_edge('c', 'e', weight=0.7)
# G.add_edge('c', 'f', weight=0.9)
# G.add_edge('a', 'd', weight=0.3)
#
# elarge = [(u, v) for (u, v, d) in G.edges(data=True) if d['weight'] > 0.5]
# esmall = [(u, v) for (u, v, d) in G.edges(data=True) if d['weight'] <= 0.5]
#
# pos = nx.spring_layout(G)  # positions for all nodes
#
# # nodes
# nx.draw_networkx_nodes(G, pos, node_size=7)
#
# # edges
# nx.draw_networkx_edges(G, pos, edgelist=elarge,
#                        width=6)
# nx.draw_networkx_edges(G, pos, edgelist=esmall,
#                        width=6, alpha=0.5, edge_color='b', style='dashed')
#
# # labels
# nx.draw_networkx_labels(G, pos, font_size=20, font_family='sans-serif')
#
# plt.axis('off')
# plt.show()

edgelist = []
for source in x.graph.hash_graph:
    for target in x.graph.hash_graph[source]:
        edgelist.append(tuple([source, target, x.graph.hash_graph[source][target]]))
G = nx.Graph()
G.add_weighted_edges_from(edgelist)
'''
Create Edges
Add edges as disconnected lines in a single trace and nodes as a scatter trace
'''
pos_dict = nx.spring_layout(G)
edge_x = []
edge_y = []
for edge in G.edges():
    x0, y0 = pos_dict[edge[0]]
    x1, y1 = pos_dict[edge[1]]
    edge_x.append(x0)
    edge_x.append(x1)
    edge_x.append(None)
    edge_y.append(y0)
    edge_y.append(y1)
    edge_y.append(None)

edge_trace1 = go.Scatter(
    x=edge_x[:2000], y=edge_y[:2000],
    line=dict(width=0.5, color='#888'),
    hoverinfo='none',
    mode='lines')

edge_trace2 = go.Scatter(
    x=edge_x[2000:], y=edge_y[2000:],
    line=dict(width=0.5, color='#888', dash="dash"),
    hoverinfo='none',
    mode='lines')

# edge_styles = []
# for i in range(len(edgelist)):
#     if i%2 ==0:
#         edge_styles.append("dash")
#     else:
#         edge_styles.append("solid")
# edge_trace.line.dash = edge_styles

node_x = []
node_y = []
for node in G.nodes():
    # x, y = G.node[node]['pos']
    x, y = pos_dict[node]#['pos']
    node_x.append(x)
    node_y.append(y)

node_trace = go.Scatter(
    x=node_x, y=node_y,
    mode='markers',
    hoverinfo='text',
    marker=dict(
        showscale=True,
        # colorscale options
        #'Greys' | 'YlGnBu' | 'Greens' | 'YlOrRd' | 'Bluered' | 'RdBu' |
        #'Reds' | 'Blues' | 'Picnic' | 'Rainbow' | 'Portland' | 'Jet' |
        #'Hot' | 'Blackbody' | 'Earth' | 'Electric' | 'Viridis' |
        colorscale='YlGnBu',
        reversescale=True,
        color=[],
        size=100,
        colorbar=dict(
            thickness=15,
            title='Node Connections',
            xanchor='left',
            titleside='right'
        ),
        line_width=2))
# Color Node Points
# Color node points by the number of connections.
#
# Another option would be to size points by the number of connections i.e. node_trace.marker.size = node_adjacencies

node_adjacencies = []
node_text = []
for node, adjacencies in enumerate(G.adjacency()):
    node_adjacencies.append(len(adjacencies[1]))
    node_text.append('%s # of connections: '%my_data.graph.id_to_name[node] +str(len(adjacencies[1])))



node_colorings = []
for i in range(len(node_adjacencies)):
    if i%2 ==0:
        node_colorings.append(1)
    else:
        node_colorings.append(2)

node_trace.marker.color = node_adjacencies
node_trace.marker.size = node_adjacencies
pp.pprint(node_adjacencies)
node_trace.text = node_text
# Create Network Graph
fig = go.Figure(data=[edge_trace1,edge_trace2, node_trace],
             layout=go.Layout(
                title='<br>Network graph made with Python',
                titlefont_size=16,
                showlegend=False,
                hovermode='closest',
                margin=dict(b=20,l=5,r=5,t=40),
                annotations=[ dict(
                    text="Python code: <a href='https://plot.ly/ipython-notebooks/network-graphs/'> https://plot.ly/ipython-notebooks/network-graphs/</a>",
                    showarrow=False,
                    xref="paper", yref="paper",
                    x=0.005, y=-0.002 ) ],
                xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                yaxis=dict(showgrid=False, zeroline=False, showticklabels=False))
                )
fig.show()
time.sleep(3)
