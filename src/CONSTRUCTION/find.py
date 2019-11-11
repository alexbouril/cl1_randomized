from src.CLUSTER_STATE.cluster_state import ClusterState
from src.COMMON.cmn import *
from src.GRAPH.graph import dfs
import multiprocessing
from functools import partial
from src.QUALITY.quality import density, num_edges_inside
def find_best_density_add(cl1, cs: ClusterState):
    cs.best_change = None
    cs.best_change_score = cs.cohesiveness
    highest_inward_count = 0
    for v in cs.add_candidates:

        numerator = cs.current_cluster_weight_in + \
                    cs.add_candidates[v]._in
        denominator = cs.current_cluster_weight_in + \
                      cs.current_cluster_weight_out + \
                      cs.add_candidates[v]._out + \
                      cl1.penalty_value_per_node * (len(cs.current_cluster) + 1)
        proposed_score = numerator / denominator
        inward = cs.add_candidates[v].num_edges_to
        if inward > highest_inward_count:
            cs.best_change = v
            cs.best_change_score = proposed_score
            highest_inward_count = inward #######################################################
    return cs.best_change, cs.best_change_score


def find_best_mixed_measure_add(cl1, cs: ClusterState):
    cs.best_change = None
    cs.best_change_score = cs.cohesiveness
    original_density = density(cl1, [v for v in cs.current_cluster])
    best_adder = 0
    for v in cs.add_candidates:

        numerator = cs.current_cluster_weight_in + \
                    cs.add_candidates[v]._in
        denominator = cs.current_cluster_weight_in + \
                      cs.current_cluster_weight_out + \
                      cs.add_candidates[v]._out + \
                      cl1.penalty_value_per_node * (len(cs.current_cluster) + 1)
        proposed_score = numerator / denominator

        density_gain = density(cl1, [x for x in cs.current_cluster]+[v]) - original_density
        cohesiveness_gain = proposed_score -cs.cohesiveness
        adder = density_gain+cohesiveness_gain
        inward = cs.add_candidates[v].num_edges_to
        if adder > best_adder:
            cs.best_change = v
            cs.best_change_score = proposed_score
            best_adder = adder
    return cs.best_change, cs.best_change_score


def find_best_mixed_measure_add2(cl1, cs: ClusterState):
    cs.best_change = None
    cs.best_change_score = cs.cohesiveness
    base_edges_inside = 0
    # TODO: modify ClusterState so that it keeps track of the number of internal edges
    for v in cs.current_cluster:
        base_edges_inside+=cs.current_cluster[v].num_edges_to############
        # base_edges_inside+=cs.current_cluster[v]._in############
    # base_edges_inside/=2
    n = len(cs.current_cluster)
    base_denominator = (n*(n-1)) /2
    new_denominator = ((n+1)*n)/2
    base_density = base_edges_inside/base_denominator
    best_adder = 0
    denom_partial = cs.current_cluster_weight_in + \
                    cs.current_cluster_weight_out + \
                    cl1.penalty_value_per_node * (len(cs.current_cluster) + 1)
    for v in cs.add_candidates:
        numerator = cs.current_cluster_weight_in + \
                    cs.add_candidates[v]._in
        denominator = denom_partial + cs.add_candidates[v]._out
        proposed_score = numerator / denominator
        inward = cs.add_candidates[v].num_edges_to#######################
        # inward = cs.add_candidates[v]._in#######################
        density_gain = (base_edges_inside+inward)/new_denominator - base_density
        cohesiveness_gain = proposed_score -cs.cohesiveness
        adder = density_gain+cohesiveness_gain
        if adder > best_adder:
            cs.best_change = v
            cs.best_change_score = proposed_score
            best_adder = adder
    return cs.best_change, cs.best_change_score


def find_best_add(cl1, cs:ClusterState):
    """
    given a ClusterState, find the add candidate whose addition to the complex will result in the greatest cohesiveness.
    If no add candidate would increase the cohesiveness, the best_change returned is None.
    """
    cs.best_change = None
    cs.best_change_score = cs.cohesiveness
    denom_partial = cs.current_cluster_weight_in+\
                    cs.current_cluster_weight_out + \
                    cl1.penalty_value_per_node * (len(cs.current_cluster) + 1)
    for v in cs.add_candidates:
        numerator = cs.current_cluster_weight_in +\
                    cs.add_candidates[v]._in
        denominator = denom_partial+ cs.add_candidates[v]._out
        proposed_score = numerator / denominator
        if proposed_score > cs.best_change_score:
            cs.best_change = v
            cs.best_change_score = proposed_score
    return cs.best_change, cs.best_change_score


def find_best_suboptimal_add(cl1, cs:ClusterState):
    """
    given a ClusterState, find the add candidate whose addition to the complex will result in the greatest cohesiveness.
    """
    cs.best_change = None
    cs.best_change_score = -10000
    denom_partial = cs.current_cluster_weight_in+\
                    cs.current_cluster_weight_out + \
                    cl1.penalty_value_per_node * (len(cs.current_cluster) + 1)
    for v in cs.add_candidates:
        numerator =  cs.current_cluster_weight_in + \
                     cs.add_candidates[v]._in
        denominator =  denom_partial + cs.add_candidates[v]._out
        proposed_score = numerator / denominator
        if proposed_score > cs.best_change_score:
            cs.best_change = v
            cs.best_change_score = proposed_score
    return cs.best_change, cs.best_change_score


def find_best_remove(cl1, cs:ClusterState):
    """
    given a ClusterState, find the remove candidate whose removal from the complex will result in the greatest cohesiveness.
    If no remove candidate would increase the cohesiveness, the best_change returned is None.
    """
    cs.best_change = None
    cs.best_change_score = cs.cohesiveness
    denom_partial = cs.current_cluster_weight_in + \
                  cs.current_cluster_weight_out + \
                  cl1.penalty_value_per_node * (len(cs.current_cluster) - 1)
    if len(cs.current_cluster) > 1:
        if cl1.care_about_cuts:
            current_cluster_membership_hashset = [vertex for vertex in cs.current_cluster]
        for v in cs.remove_candidates:
            if cl1.care_about_cuts:
                # TODO: check if there is a cut.
                #   Implement more efficiently using a Dynamic Connectivity algorithm
                is_a_cut = True
                visted = set()
                start_point = None
                for potential_start_point in current_cluster_membership_hashset:
                    if potential_start_point != v:
                        start_point = potential_start_point
                        # consider break statement here
                # check that
                #   (2) removal of vertex under consideration will not disconnect cluster
                dfs(cl1, start_point, v, current_cluster_membership_hashset, visted)
                if len(visted) == -1 + len(current_cluster_membership_hashset):
                    is_a_cut = False
                    debug("%s is NOT a CUT" % str(v))
                if is_a_cut:
                    debug("%s is a CUT!" % str(v))
            else:
                is_a_cut = False
            if not is_a_cut:
                # TODO: check that this makes sense
                numerator = cs.current_cluster_weight_in - cs.remove_candidates[v]._in
                denominator =  denom_partial - cs.remove_candidates[v]._out
                proposed_score = numerator / denominator
                if proposed_score > cs.best_change_score:
                    cs.best_change = v
                    cs.best_change_score = proposed_score
    ###################
    # sanity check: showing improvement from removing
    ##################
    # if cs.best_change:
    #     print(cs.best_change, cs.best_change_score, cs.cohesiveness)
    return cs.best_change, cs.best_change_score


def careful_find_best_2neighborhood_add(cl1, cs:ClusterState):
    cs.best_change = None
    cs.best_change_score = cs.cohesiveness
    best_proposed_score = cs.cohesiveness
    neighborhood_gain = {v: {"in": 0,
                            "out": 0}
                                        for v in cs.add_candidates}
    good_neighbors = dict()
    for v in cs.add_candidates:
        numerator = cs.current_cluster_weight_in + \
                    cs.add_candidates[v]._in
        denominator = cs.current_cluster_weight_in + \
                      cs.current_cluster_weight_out + \
                      cl1.penalty_value_per_node * (len(cs.current_cluster) + 1)+\
                      cs.add_candidates[v]._out
        actual_score = numerator/denominator
        #################################################################
        #  CONSIDER WHAT COULD BE GAINED ON THE NEWLY EXPOSED BOUNDARY  #
        #################################################################
        distance_2_neighbors = cl1.graph.hash_graph[v]
        in_and_out_for_distance_2_neighbors = dict()
        count_d2n = 0
        count_d3n = 0
        for d2n in distance_2_neighbors:
            count_d2n += 1
            if d2n in cs.current_cluster:
                '''already captured in the term current_cluster_weight_in'''
                continue
            else:
                in_and_out_for_distance_2_neighbors[d2n]={"to_current_cluster":0,
                                  "to_v":0,
                                  "to_N(v)": 0,
                                  "out":0}
                distance_3_neighbors = cl1.graph.hash_graph[d2n]
                for d3n in distance_3_neighbors:
                    count_d3n+=1
                    edge_weight = cl1.graph.hash_graph[d2n][d3n]
                    if d3n is v:
                        in_and_out_for_distance_2_neighbors[d2n]["to_v"] += 3*edge_weight
                    # Todo: this is calculated multiple times
                    if d3n in cs.current_cluster:
                        in_and_out_for_distance_2_neighbors[d2n]["to_current_cluster"] += 3 * edge_weight
                    elif d3n in distance_2_neighbors:
                        in_and_out_for_distance_2_neighbors[d2n]["to_N(v)"] = edge_weight
                    else:
                        in_and_out_for_distance_2_neighbors[d2n]["out"] = edge_weight
        d2n_to_grab = set()
        ####################################
        # get an idea of how many distance-2 and distance-3 neighbors are considered
        ####################################
        # print(count_d2n, count_d3n)
        for d2n in in_and_out_for_distance_2_neighbors:
            # TODO: figure out what to do with "to_N(v)"
            inbound = in_and_out_for_distance_2_neighbors[d2n]["to_current_cluster"] + \
                      in_and_out_for_distance_2_neighbors[d2n]["to_v"]
            outbound = in_and_out_for_distance_2_neighbors[d2n]["out"]
            if inbound>=outbound:
            # if inbound/(inbound+outbound)>actual_score:
                d2n_to_grab.add(d2n)
                neighborhood_gain[v]["in"]+=inbound
                neighborhood_gain[v]["out"]+=outbound
        good_neighbors[v] = d2n_to_grab
        # factor = len(current_cluster)*.09
        # factor = logistic.cdf(len(current_cluster)/10)
        # factor = .2* logistic.cdf(len(current_cluster))
        factor = len(cs.current_cluster)/10
        numerator+= factor*neighborhood_gain[v]["in"]
        denominator+= factor*\
                      (neighborhood_gain[v]["in"] + neighborhood_gain[v]["out"])
        proposed_score = numerator / denominator
        if proposed_score > best_proposed_score:
            cs.best_change = v
            cs.best_change_score = actual_score
            best_proposed_score = proposed_score
    return cs.best_change, cs.best_change_score, good_neighbors.get(cs.best_change, None)
