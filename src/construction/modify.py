from src.CL1R.cl1r import CL1_Randomized
from src.COMMON.cmn import *

def add(cl1:CL1_Randomized, cs):
    #################################################################
    # update the overall weight into and out of the current_cluster #
    #################################################################
    change_vertex_in = cs.add_candidates[cs.best_change]._in
    change_vertex_out = cs.add_candidates[cs.best_change]._out
    cs.current_cluster_weight_in += change_vertex_in
    cs.current_cluster_weight_out = cs.current_cluster_weight_out - change_vertex_in + change_vertex_out

    ###################################
    # Move the change vertex from cs.add_candidates to current_cluster
    ###################################
    to_add = cs.add_candidates[cs.best_change].copy()
    del cs.add_candidates[cs.best_change]
    cs.current_cluster[cs.best_change] = to_add.copy()

    ###################################
    # Change vertex to remove_candidates if applicable
    ###################################
    if to_add.num_edges_from:
        cs.remove_candidates[cs.best_change] = to_add.copy()

    def update_v(v, edge_weight, collection):
        collection[v]._in += edge_weight
        collection[v]._out -= edge_weight
        collection[v].num_edges_to += 1
        collection[v].num_edges_from -= 1
        #######################
        # sanity check
        #######################
        thresh = -.001
        a = collection[v]._in < thresh
        b = collection[v]._out < thresh
        c = collection[v].num_edges_to < thresh
        d = collection[v].num_edges_from < thresh
        if a or b or c or d:
            print("oh no %s; %s%s%s%s"%(str(v), a, b, c, d))
            exit()

    def initialize_new_add_candidate_for_V(add_vertex, cl1, cs):
        num_edges_to = 0
        weight_to = 0
        num_edges_from = 0
        weight_from = 0
        ###################################
        # iterate over the neighbors of v
        ###################################
        for neighbor in cl1.graph.hash_graph[add_vertex]:
            weight_prime = cl1.graph.hash_graph[add_vertex][neighbor]
            if neighbor in cs.current_cluster:
                num_edges_to += 1
                weight_to += weight_prime
            else:
                num_edges_from += 1
                weight_from += weight_prime

        cs.add_candidates[add_vertex] = Relationship(weight_to,
                                         num_edges_to,
                                         weight_from,
                                         num_edges_from)
    #######################################################################
    # iterate over neighbors of change_vertex, and update each Relationship
    #######################################################################
    for v in cl1.graph.hash_graph[cs.best_change]:
        edge_weight = cl1.graph.hash_graph[v][cs.best_change]
        if v in cs.add_candidates:
            update_v(v, edge_weight, cs.add_candidates)
        if v in cs.current_cluster:
            update_v(v, edge_weight, cs.current_cluster)
        # note that v may be in both the current_cluster and in remove_candidates
        # remove_candidates is a subset of current_cluster
        if v in cs.remove_candidates:
            update_v(v, edge_weight, cs.remove_candidates)
            # Check that a candidate for removal is still on the boundary
            if cs.remove_candidates[v].num_edges_from == 0:
                del cs.remove_candidates[v]
        # handle the case that v is on the new boundary
        # add v to add_candidates
        if v not in cs.add_candidates and v not in cs.current_cluster:
            initialize_new_add_candidate_for_V(v, cl1, cs)
    return cs.current_cluster_weight_in, cs.current_cluster_weight_out



def remove(cl1: CL1_Randomized, cs:ClusterState):
    ###############################################################################################
    # Update the current_cluster 's score, and overall weight into and out of the current cluster #
    ###############################################################################################
    change_vertex_in = cs.remove_candidates[cs.best_change]._in
    change_vertex_out = cs.remove_candidates[cs.best_change]._out
    cs.current_cluster_weight_in -= change_vertex_in
    cs.current_cluster_weight_out = cs.current_cluster_weight_out - change_vertex_out + change_vertex_in
    #######################################################################
    # Remove the change vertex from remove_candidates and current_cluster #
    #######################################################################
    to_remove = cs.remove_candidates[cs.best_change].copy()
    del cs.remove_candidates[cs.best_change]
    del cs.current_cluster[cs.best_change]
    ################################################
    # Also add the change vertex to add_candidates #
    ################################################
    cs.add_candidates[cs.best_change] = to_remove
    #########################################################################
    # iterate over neighbors of change_vertex, and update each Relationship #
    #########################################################################
    def update_v(v, edge_weight, collection):
        collection[v]._in -= edge_weight
        collection[v]._out += edge_weight
        collection[v].num_edges_to -= 1
        collection[v].num_edges_from += 1

    for v in cl1.graph.hash_graph[cs.best_change]:
        edge_weight = cl1.graph.hash_graph[cs.best_change][v]
        # note that v may be in both the current_cluster and in remove_candidates
        if v in cs.remove_candidates:
            update_v(v, edge_weight, cs.remove_candidates)
        if v in cs.current_cluster:
            update_v(v, edge_weight, cs.current_cluster)
            if cs.current_cluster[v].num_edges_from == 1:
                cs.remove_candidates[v] = cs.current_cluster[v].copy()
        if v in cs.add_candidates:
            update_v(v, edge_weight, cs.add_candidates)
            if cs.add_candidates[v].num_edges_to == 0:
                del cs.add_candidates[v]
    return cs.current_cluster_weight_in, cs.current_cluster_weight_out