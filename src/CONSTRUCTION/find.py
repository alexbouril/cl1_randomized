from src.CLUSTER_STATE.cluster_state import ClusterState
from src.COMMON.cmn import *
from src.GRAPH.graph import dfs
import multiprocessing
from functools import partial

def find_best_add(cl1, cs:ClusterState):
    cs.best_change = None
    cs.best_change_score = cs.cohesiveness
    for v in cs.add_candidates:
        numerator = cs.current_cluster_weight_in +\
                    cs.add_candidates[v]._in
        denominator = cs.current_cluster_weight_in+\
                      cs.current_cluster_weight_out +\
                      cs.add_candidates[v]._out +\
                      cl1.penalty_value_per_node * (len(cs.current_cluster) + 1)
        proposed_score = numerator / denominator
        if proposed_score > cs.best_change_score:
            cs.best_change = v
            cs.best_change_score = proposed_score
    return cs.best_change, cs.best_change_score


def find_best_suboptimal_add(cl1, cs:ClusterState):
    cs.best_change = None
    cs.best_change_score = -10000
    for v in cs.add_candidates:
        numerator =  cs.current_cluster_weight_in + \
                     cs.add_candidates[v]._in
        denominator = cs.current_cluster_weight_in +\
                      cs.current_cluster_weight_out+ \
                      cs.add_candidates[v]._out +\
                      cl1.penalty_value_per_node * (len(cs.current_cluster) + 1)
        proposed_score = numerator / denominator
        if proposed_score > cs.best_change_score:
            cs.best_change = v
            cs.best_change_score = proposed_score
    return cs.best_change, cs.best_change_score


def find_best_remove(cl1, cs:ClusterState):
    cs.best_change = None
    cs.best_change_score = cs.cohesiveness##########TODO: uncomment
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
                denominator = cs.current_cluster_weight_in + \
                              cs.current_cluster_weight_out - \
                              cs.remove_candidates[v]._out +\
                              cl1.penalty_value_per_node * (len(cs.current_cluster) - 1)
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




def helper_careful_find_best_2neighborhood_add(cl1, cs:ClusterState, v):
    # print(list(map(type, [cl1, cs, v])))
    v = v[0]
    # print(v)
    good_neighbors = dict()
    neighborhood_gain = {v: {"in": 0,
                            "out": 0}}
    numerator = cs.current_cluster_weight_in + \
                cs.add_candidates[v]._in
    denominator = cs.current_cluster_weight_in + \
                  cs.current_cluster_weight_out + \
                  cl1.penalty_value_per_node * (len(cs.current_cluster) + 1) + \
                  cs.add_candidates[v]._out
    actual_score = numerator / denominator
    #################################################################
    #  CONSIDER WHAT COULD BE GAINED ON THE NEWLY EXPOSED BOUNDARY  #
    #################################################################
    distance_2_neighbors = cl1.graph.hash_graph[v]
    in_and_out_for_distance_2_neighbors = dict()
    for d2n in distance_2_neighbors:
        if d2n in cs.current_cluster:
            '''already captured in the term current_cluster_weight_in'''
            continue
        else:
            in_and_out_for_distance_2_neighbors[d2n] = {"to_current_cluster": 0,
                                                        "to_v": 0,
                                                        "to_N(v)": 0,
                                                        "out": 0}
            distance_3_neighbors = cl1.graph.hash_graph[d2n]
            for d3n in distance_3_neighbors:
                edge_weight = cl1.graph.hash_graph[d2n][d3n]
                if d3n is v:
                    in_and_out_for_distance_2_neighbors[d2n]["to_v"] += 3 * edge_weight
                if d3n in cs.current_cluster:
                    in_and_out_for_distance_2_neighbors[d2n]["to_current_cluster"] += 3 * edge_weight
                elif d3n in distance_2_neighbors:
                    in_and_out_for_distance_2_neighbors[d2n]["to_N(v)"] = edge_weight
                else:
                    in_and_out_for_distance_2_neighbors[d2n]["out"] = edge_weight
    d2n_to_grab = set()
    for d2n in in_and_out_for_distance_2_neighbors:
        # TODO: figure out what to do with "to_N(v)"
        inbound = in_and_out_for_distance_2_neighbors[d2n]["to_current_cluster"] + \
                  in_and_out_for_distance_2_neighbors[d2n]["to_v"]
        outbound = in_and_out_for_distance_2_neighbors[d2n]["out"]
        if inbound >= outbound:
            # if inbound/(inbound+outbound)>actual_score:
            d2n_to_grab.add(d2n)
            neighborhood_gain[v]["in"] += inbound
            neighborhood_gain[v]["out"] += outbound
    good_neighbors[v] = d2n_to_grab
    # factor = len(current_cluster)*.09
    # factor = logistic.cdf(len(current_cluster)/10)
    # factor = .2* logistic.cdf(len(current_cluster))
    factor = len(cs.current_cluster) / 10
    numerator += factor * neighborhood_gain[v]["in"]
    denominator += factor * \
                   (neighborhood_gain[v]["in"] + neighborhood_gain[v]["out"])
    proposed_score = numerator / denominator
    # if proposed_score > best_proposed_score:
    #     cs.best_change = v
    #     cs.best_change_score = actual_score
    #     best_proposed_score = proposed_score
    return proposed_score, actual_score




def multiprocessing_careful_find_best_2neighborhood_add(cl1, cs:ClusterState):
    cs.best_change = None
    cs.best_change_score = cs.cohesiveness
    best_proposed_score = cs.cohesiveness

    li = [(v,cl1, cs) for v in cs.add_candidates]
    p = multiprocessing.Pool()
    func = partial(helper_careful_find_best_2neighborhood_add, cl1, cs)
    result = p.map(func, li)
    p.close()
    p.join()

    for i, (proposed_score, actual_score) in enumerate(result):
        if proposed_score>best_proposed_score:
            best_proposed_score=proposed_score
            cs.best_change = li[i][0]
            cs.best_change_score = actual_score
    # print("_____________________-")
    # print(cs.best_change)
    # print(cs.best_change_score)
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
        for d2n in distance_2_neighbors:
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
                    edge_weight = cl1.graph.hash_graph[d2n][d3n]
                    if d3n is v:
                        in_and_out_for_distance_2_neighbors[d2n]["to_v"] += 3*edge_weight
                    if d3n in cs.current_cluster:
                        in_and_out_for_distance_2_neighbors[d2n]["to_current_cluster"] += 3 * edge_weight
                    elif d3n in distance_2_neighbors:
                        in_and_out_for_distance_2_neighbors[d2n]["to_N(v)"] = edge_weight
                    else:
                        in_and_out_for_distance_2_neighbors[d2n]["out"] = edge_weight
        d2n_to_grab = set()
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
