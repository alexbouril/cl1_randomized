from src.COMMON.cmn import *

def initialize_complex(cl1, current_seed):
    ###########################################
    # initalize the current cluster
    ###########################################
    current_cluster = dict()
    seed_weight_to = 0
    seed_num_edges_to = 0
    seed_weight_from = sum([cl1.graph.hash_graph[current_seed][tar] for tar in cl1.graph.hash_graph[current_seed]])
    seed_num_edges_from = len(cl1.graph.hash_graph[current_seed])
    current_cluster[current_seed] = Relationship(seed_weight_to,
                                                 seed_num_edges_to,
                                                 seed_weight_from,
                                                 seed_num_edges_from)
    current_score = 0
    current_cluster_weight_in = 0
    current_cluster_weight_out = seed_weight_from

    ###########################################
    # initalize the candidates for removal
    ###########################################
    remove_candidates = dict()
    remove_candidates[current_seed] = current_cluster[current_seed].copy()

    ###########################################
    # initialize the candidates for addition
    ###########################################
    add_candidates = dict()
    for target in cl1.graph.hash_graph[current_seed]:
        target_weight_to = cl1.graph.hash_graph[current_seed][target]
        target_num_edges = len(cl1.graph.hash_graph[target])
        target_num_edges_to = 1
        target_num_edges_from = target_num_edges - 1
        target_weight_from = sum([cl1.graph.hash_graph[target][tar] for tar in cl1.graph.hash_graph[target] if
                           tar != current_seed])
        add_candidates[target] = Relationship(target_weight_to,
                                              target_num_edges_to,
                                              target_weight_from,
                                              target_num_edges_from)
    return current_cluster,   \
           remove_candidates, \
           add_candidates,    \
           current_score,     \
           current_cluster_weight_in, \
           current_cluster_weight_out






