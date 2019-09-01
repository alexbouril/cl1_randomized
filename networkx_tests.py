from cl1_randomized import *
from common import *
from graph import *
import plotly.graph_objects as go
import networkx as nx

n = "gavin2006_socioaffinities_rescaled+2019-08-31_17:14:11:720422"
name = \
"pickles/pickle+"+n


x = loadData(name)
id_to_name = x.graph.id_to_name
print(stringify_construction_log(x.construction_log))


edgelist = []
for source in x.graph.hash_graph:
    for target in x.graph.hash_graph[source]:
        edgelist.append(tuple([source, target, x.graph.hash_graph[source][target]]))


cluster_names = [key for key in x.construction_log]
current_cluster = cluster_names[0]
current_cluster_construction_log = x.construction_log[current_cluster]
'''
Create Edges
Add edges as disconnected lines in a single trace and nodes as a scatter trace
'''

##############################################
# Find all Nodes involved in this cluster's construction
##############################################
all_nodes = set()
for cs in [entry for entry in current_cluster_construction_log if type(entry)==ClusterState]:
    for id in [name for name in cs.current_cluster] + [name for name in cs.add_candidates] + [name for name in cs.remove_candidates]:
        all_nodes.add(id)
all_nodes = list(all_nodes)
print(all_nodes)
##############################################
# Find all the edges involved in this cluster's construction
##############################################
all_edges = set()
for i in all_nodes:
    for j in all_nodes:
        if j in x.graph.hash_graph[i]:
            a, b = sorted([i, j])
            all_edges.add(tuple([a, b, x.graph.hash_graph[i][j]]))
all_edges = list(all_edges)
print(all_edges)
##############################################
# Create a graph G_plus will all the edges and nodes involved in this cluster's construction
##############################################
G_plus= nx.Graph()
G_plus.add_weighted_edges_from(all_edges)

##############################################
# Create the pos_dict using G_plus
##############################################
pos_dict = nx.spring_layout(G_plus)

##############################################
##############################################
##############################################
##############################################

for entry in current_cluster_construction_log:
    if type(entry)==ClusterState:
        print('e')
cluster_states = [entry for entry in current_cluster_construction_log if type(entry)==ClusterState]

# cs = cluster_states[3]
for cs in cluster_states[:10]:
    cs_add_candidates=cs.add_candidates
    cs_remove_candidates=cs.remove_candidates
    cs_current_cluster = cs.current_cluster


    nodes_inside_cluster = [key for key in cs_current_cluster]
    nodes_outside_cluster = [key for key in cs_add_candidates]

    all_involved_edges = set()
    for i in nodes_inside_cluster+nodes_outside_cluster:
        for j in nodes_inside_cluster+nodes_outside_cluster:
            if j in x.graph.hash_graph[i]:
                a,b = sorted([i,j])
                all_involved_edges.add(tuple([a, b, x.graph.hash_graph[i][j]]))
    all_involved_edges = list(all_involved_edges)


    edges_inside_cluster_x = []
    edges_inside_cluster_y = []
    edges_outside_cluster_x = []
    edges_outside_cluster_y = []


    G = nx.Graph()
    G.add_weighted_edges_from(all_involved_edges)

    # pos_dict = nx.spring_layout(G)
    edge_x = []
    edge_y = []
    for edge in G.edges():
        x0, y0 = pos_dict[edge[0]]
        x1, y1 = pos_dict[edge[1]]
        if edge[0] in nodes_inside_cluster and edge[1] in nodes_inside_cluster:
            edges_inside_cluster_x.append(x0)
            edges_inside_cluster_x.append(x1)
            edges_inside_cluster_x.append(None)
            edges_inside_cluster_y.append(y0)
            edges_inside_cluster_y.append(y1)
            edges_inside_cluster_y.append(None)
        else:
            edges_outside_cluster_x.append(x0)
            edges_outside_cluster_x.append(x1)
            edges_outside_cluster_x.append(None)
            edges_outside_cluster_y.append(y0)
            edges_outside_cluster_y.append(y1)
            edges_outside_cluster_y.append(None)

    edge_trace_inside_cluster = go.Scatter(
        x=edges_inside_cluster_x, y=edges_inside_cluster_y,
        line=dict(width=1, color='black'),
        hoverinfo='none',
        mode='lines')

    edge_trace2_outside_cluster = go.Scatter(
        x=edges_outside_cluster_x, y=edges_outside_cluster_y,
        line=dict(width=0.5, color='red', dash="dash"),
        hoverinfo='none',
        mode='lines')


    node_x = []
    node_y = []
    for node in G.nodes():
        # x, y = G.node[node]['pos']
        pos_x, pos_y = pos_dict[node]#['pos']
        node_x.append(pos_x)
        node_y.append(pos_y)

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
        node_text.append('%s # of connections: '%id_to_name[node] +str(len(adjacencies[1])))



    node_colorings = []
    for i in range(len(node_adjacencies)):
        if i%2 ==0:
            node_colorings.append(1)
        else:
            node_colorings.append(2)

    node_trace.marker.color = node_adjacencies
    node_trace.marker.size = node_adjacencies
    # pp.pprint(node_adjacencies)
    node_trace.text = node_text
    # Create Network Graph
    fig = go.Figure(data=[edge_trace_inside_cluster,edge_trace2_outside_cluster, node_trace],
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
