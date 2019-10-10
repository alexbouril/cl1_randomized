from src.cl1_randomized.cl1_randomized import *

def careful_find_best_2neighborhood_add(cl1: CL1_Randomized, cs:ClusterState):
    # if len(current_cluster)<5:
    #     return find_best_add(self, add_candidates, current_cluster, current_score, current_cluster_weight_in, current_cluster_weight_out)
    # best_change_list = find_best_add_list(self, add_candidates, current_cluster, current_score, current_cluster_weight_in, current_cluster_weight_out)

    cs.best_change = None
    best_change_score = cs.cohesiveness
    best_proposed_score = cs.cohesiveness
    neighborhood_gain = {v: {"in": 0,
                            "out": 0}
                                        for v in cs.add_candidates}
    good_neighbors = dict()
    for v in cs.add_candidates:
    # for v in best_change_list:
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
                #########################################
                #  THIS HAS ALREADY BEEN ACCOUNTED FOR  #
                #########################################
                '''
                already captured in the term current_cluster_weight_in
                '''
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
        # in_N_v_grabbed = 0
        # for a in d2n_to_grab:
        #     for b in d2n_to_grab:
        #         if b in self.graph.hash_graph[a]:
        #             in_N_v_grabbed += self.graph.hash_graph[a][b]
        # dict_gain[v]["in"]+=in_N_v_grabbed
        # print(v, d2n_to_grab)
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
        sleep_debug(.25)
    # print("-------------",best_change)
    return cs.best_change, cs.best_change_score, good_neighbors[cs.best_change]


def initialize_complex(self, current_seed):
    ###########################################
    # initalize the current cluster
    ###########################################
    current_cluster = dict()
    seed_weight_to = 0
    seed_num_edges_to = 0
    seed_weight_from = sum([self.graph.hash_graph[current_seed][tar] for tar in self.graph.hash_graph[current_seed]])
    seed_num_edges_from = len(self.graph.hash_graph[current_seed])
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
    for target in self.graph.hash_graph[current_seed]:
        target_weight_to = self.graph.hash_graph[current_seed][target]
        target_num_edges = len(self.graph.hash_graph[target])
        target_num_edges_to = 1
        target_num_edges_from = target_num_edges - 1
        target_weight_from = sum([self.graph.hash_graph[target][tar] for tar in self.graph.hash_graph[target] if
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


def find_best_add(cl1:CL1_Randomized, cs:ClusterState):
    cs.best_change = None
    best_change_score = cs.cohesiveness
    for v in cs.add_candidates:
        numerator = cs.current_cluster_weight_in +\
                    cs.add_candidates[v]._in
        denominator = cs.current_cluster_weight_in+\
                      cs.current_cluster_weight_out +\
                      cs.add_candidates[v]._out +\
                      cl1.penalty_value_per_node * (len(cs.current_cluster) + 1)
        proposed_score = numerator / denominator
        if proposed_score > best_change_score:
            best_change = v
            best_change_score = proposed_score
    return best_change, best_change_score


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


def find_best_remove(cl1:CL1_Randomized, cs:ClusterState):
    cs.best_change = None
    # check that
    #   (1) the cluster has more than one element
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
                debug("cluster: %s" % str(current_cluster_membership_hashset))
                debug("visited by DFS: %s" % str(visted))
            else:
                is_a_cut = False

            if not is_a_cut:
                # TODO: check that this makes sense
                numerator = cs.current_cluster_weight_in - cs.remove_candidates[v]._in
                denominator = cs.current_cluster_weight_in + \
                              cs.current_cluster_weight_out - \
                              cs.remove_candidates[v]._out +\
                              cl1.penalty_value_per_node * (len(cl1.current_cluster) - 1)
                proposed_score = numerator / denominator
                if proposed_score > cs.cohesiveness:
                    cs.best_change = v
                    cs.cohesiveness = proposed_score
    return cs.best_change, cs.cohesiveness()


def remove(cl1: CL1_Randomized, cs:ClusterState):
    ###############################################################################################
    # Update the current_cluster 's score, and overall weight into and out of the current cluster #
    ###############################################################################################
    change_vertex_in = cs.remove_candidates[cs.change_vertex]._in
    change_vertex_out = cs.remove_candidates[cs.change_vertex]._out
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


def find_best_suboptimal_add(cl1:CL1_Randomized, cs:ClusterState):
    cs.best_change = None
    cs.best_change_score = -10000
    for v in cs.add_candidates:
        numerator =  cs.current_cluster_weight_in + cs.add_candidates[v]._in
        denominator = cs.current_cluster_weight_in +\
                      cs.current_cluster_weight_out+ \
                      cs.add_candidates[v]._out +\
                      cl1.penalty_value_per_node * (len(cs.current_cluster) + 1)
        proposed_score = numerator / denominator
        if proposed_score > cs.best_change:
            cs.best_change = v
            cs.best_change_score = proposed_score
    return cs.best_change, cs.best_change_score


# def add_shake(self, add_candidates, remove_candidates, current_cluster, cc_weight_in, cc_weight_out, round_no, last_failed_add_round_no):
#     for i in range(self.number_of_bad_adds):
#         best_suboptimal_change, best_suboptimal_score = find_best_suboptimal_add(self, add_candidates, current_cluster,
#                                                                                  cc_weight_in, cc_weight_out)
#         if best_suboptimal_change:
#             print("adding_suboptimally")
#             sleep_debug(1)
#             round_no += 1
#             last_failed_add_round_no = -5
#             # TODO handle round_no numbers
#             cc_weight_in, cc_weight_out = add(self, add_candidates, current_cluster, remove_candidates,
#                                               best_suboptimal_change, best_suboptimal_score, cc_weight_in,
#                                               cc_weight_out)
#     return best_suboptimal_score, \
#            cc_weight_in, \
#            cc_weight_out, \
#            round_no, \
#            last_failed_add_round_no


# from graph import *
# from src.common.common import *
#
# def find_best_add_list(self, add_candidates, current_cluster, current_score, current_cluster_weight_in, current_cluster_weight_out):
#     best_change_score = current_score
#     best_change_list = []
#     for v in add_candidates:
#         numerator = current_cluster_weight_in + add_candidates[v]._in
#         denominator = current_cluster_weight_in + current_cluster_weight_out + add_candidates[
#             v]._out + self.penalty_value_per_node * (len(current_cluster) + 1)
#         proposed_score = numerator / denominator
#         if proposed_score > best_change_score:
#             best_change = v
#             best_change_score = proposed_score
#             best_change_list.append((proposed_score, best_change))
#
#     def filter_list(best_change_list, length_to_keep):
#         for index, tup in enumerate(best_change_list):
#             new_tuple = (-1 * tup[0], tup[1])
#             best_change_list[index] = new_tuple
#
#         length_to_keep = min(len(best_change_list), length_to_keep)
#         heapq.heapify(best_change_list)
#         update_list = []
#         for i in range(length_to_keep):
#             update_list.append(heapq.heappop(best_change_list))
#         best_change_list = update_list
#
#     print(best_change_list)
#     filter_list(best_change_list, 7)
#     retval = [t[1] for t in best_change_list]
#     return retval
#
#
# def careful_find_best_2neighborhood_add(self, add_candidates, current_cluster, current_score, current_cluster_weight_in, current_cluster_weight_out):
#     # if len(current_cluster)<5:
#     #     return find_best_add(self, add_candidates, current_cluster, current_score, current_cluster_weight_in, current_cluster_weight_out)
#     # best_change_list = find_best_add_list(self, add_candidates, current_cluster, current_score, current_cluster_weight_in, current_cluster_weight_out)
#
#     best_change = None
#     best_change_score = current_score
#     best_proposed_score = current_score
#     neighborhood_gain = {v: {"in": 0,
#                             "out": 0}
#                                         for v in add_candidates}
#     good_neighbors = dict()
#     for v in add_candidates:
#     # for v in best_change_list:
#         numerator = current_cluster_weight_in + \
#                     add_candidates[v]._in
#         denominator = current_cluster_weight_in + \
#                       current_cluster_weight_out + \
#                       self.penalty_value_per_node * (len(current_cluster) + 1)+\
#                       add_candidates[v]._out
#
#         actual_score = numerator/denominator
#         #################################################################
#         #  CONSIDER WHAT COULD BE GAINED ON THE NEWLY EXPOSED BOUNDARY  #
#         #################################################################
#         distance_2_neighbors = self.graph.hash_graph[v]
#         in_and_out_for_distance_2_neighbors = dict()
#         for d2n in distance_2_neighbors:
#             if d2n in current_cluster:
#                 #########################################
#                 #  THIS HAS ALREADY BEEN ACCOUNTED FOR  #
#                 #########################################
#                 '''
#                 already captured in the term current_cluster_weight_in
#                 '''
#                 continue
#             else:
#                 in_and_out_for_distance_2_neighbors[d2n]={"to_current_cluster":0,
#                                   "to_v":0,
#                                   "to_N(v)": 0,
#                                   "out":0}
#                 distance_3_neighbors = self.graph.hash_graph[d2n]
#                 for d3n in distance_3_neighbors:
#                     edge_weight = self.graph.hash_graph[d2n][d3n]
#                     if d3n is v:
#                         in_and_out_for_distance_2_neighbors[d2n]["to_v"] += 3*edge_weight
#                     if d3n in current_cluster:
#                         in_and_out_for_distance_2_neighbors[d2n]["to_current_cluster"] += 3 * edge_weight
#                     elif d3n in distance_2_neighbors:
#                         in_and_out_for_distance_2_neighbors[d2n]["to_N(v)"] = edge_weight
#                     else:
#                         in_and_out_for_distance_2_neighbors[d2n]["out"] = edge_weight
#         d2n_to_grab = set()
#         for d2n in in_and_out_for_distance_2_neighbors:
#             # TODO: figure out what to do with "to_N(v)"
#             inbound = in_and_out_for_distance_2_neighbors[d2n]["to_current_cluster"] + \
#                       in_and_out_for_distance_2_neighbors[d2n]["to_v"]
#             outbound = in_and_out_for_distance_2_neighbors[d2n]["out"]
#             if inbound>=outbound:
#             # if inbound/(inbound+outbound)>actual_score:
#                 d2n_to_grab.add(d2n)
#                 neighborhood_gain[v]["in"]+=inbound
#                 neighborhood_gain[v]["out"]+=outbound
#         good_neighbors[v] = d2n_to_grab
#         # in_N_v_grabbed = 0
#         # for a in d2n_to_grab:
#         #     for b in d2n_to_grab:
#         #         if b in self.graph.hash_graph[a]:
#         #             in_N_v_grabbed += self.graph.hash_graph[a][b]
#         # dict_gain[v]["in"]+=in_N_v_grabbed
#         # print(v, d2n_to_grab)
#         # factor = len(current_cluster)*.09
#         # factor = logistic.cdf(len(current_cluster)/10)
#         # factor = .2* logistic.cdf(len(current_cluster))
#         factor = len(current_cluster)/10
#         numerator+= factor*neighborhood_gain[v]["in"]
#         denominator+= factor*\
#                       (neighborhood_gain[v]["in"] + neighborhood_gain[v]["out"])
#
#         proposed_score = numerator / denominator
#         if proposed_score > best_proposed_score:
#             best_change = v
#             best_proposed_score = proposed_score
#             best_change_score = actual_score
#         sleep_debug(.25)
#     # print("-------------",best_change)
#     return best_change, best_change_score, #good_neighbors[best_change]
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
# def find_best_2neighborhood_add(self, add_candidates, current_cluster, current_score, current_cluster_weight_in, current_cluster_weight_out):
#     best_change = None
#     best_change_score = current_score
#     factor = .3
#     for v in add_candidates:
#         numerator = current_cluster_weight_in + add_candidates[v]._in
#         denominator = current_cluster_weight_in + current_cluster_weight_out + add_candidates[
#             v]._out + self.penalty_value_per_node * (len(current_cluster) + 1)
#         #####################
#         # CONSIDER WHAT COULD BE GAINED ON THE NEWLY EXPOSED BOUNDARY
#         #####################
#         # find the neighbors of v
#         distance_2_neighbors = self.graph.hash_graph[v]
#         for d2n in distance_2_neighbors:
#             if d2n in current_cluster:
#                 numerator += factor * self.graph.hash_graph[v][d2n]
#             else:
#                 distance_3_neighbors = self.graph.hash_graph[d2n]
#                 for d3n in distance_3_neighbors:
#                     if d3n in current_cluster:
#                         numerator += factor * self.graph.hash_graph[d2n][d3n]
#                     # elif d3n in add_candidates:
#                     #     numerator += self.graph.hash_graph[d2n][d3n]
#                     else:
#                         denominator += factor * self.graph.hash_graph[d2n][d3n]
#
#         proposed_score = numerator / denominator
#         if proposed_score > best_change_score:
#             best_change = v
#             best_change_score = proposed_score
#         debug("##################### ADD Consideration ########################")
#         debug("v: %s" % str(v))
#         debug("proposed_score: %s" % str(proposed_score))
#         debug("best_change_score: %s" % str(best_change_score))
#         debug("best_change: %s" % str(best_change))
#         debug("numerator: %s" % str(numerator))
#         debug("denominator: %s" % str(denominator))
#         debug("current_cluster_weight_in: %s" % str(current_cluster_weight_in))
#         debug("add_candidates[v]._in: %s" % str(add_candidates[v]._in))
#         debug("current_cluster_weight_out: %s" % str(current_cluster_weight_out))
#         debug("add_candidates[v]._out: %s" % str(add_candidates[v]._out))
#         debug("len(current_cluster): %s" % str(len(current_cluster)))
#         sleep_debug(.25)
#     return best_change, best_change_score
#
#
#
#
#
#
#
#
#
# def initialize_complex(self, current_seed):
#     ###########################################
#     # initalize the current cluster
#     ###########################################
#     current_cluster = dict()
#     seed_weight_to = 0
#     seed_num_edges_to = 0
#     seed_weight_from = sum([self.graph.hash_graph[current_seed][tar] for tar in self.graph.hash_graph[current_seed]])
#     seed_num_edges_from = len(self.graph.hash_graph[current_seed])
#     current_cluster[current_seed] = Relationship(seed_weight_to,
#                                                  seed_num_edges_to,
#                                                  seed_weight_from,
#                                                  seed_num_edges_from)
#     current_score = 0
#     current_cluster_weight_in = 0
#     current_cluster_weight_out = seed_weight_from
#
#     ###########################################
#     # initalize the candidates for removal
#     ###########################################
#     remove_candidates = dict()
#     remove_candidates[current_seed] = current_cluster[current_seed].copy()
#
#     ###########################################
#     # initialize the candidates for addition
#     ###########################################
#     add_candidates = dict()
#     for target in self.graph.hash_graph[current_seed]:
#         target_weight_to = self.graph.hash_graph[current_seed][target]
#         target_num_edges = len(self.graph.hash_graph[target])
#         target_num_edges_to = 1
#         target_num_edges_from = target_num_edges - 1
#         target_weight_from = sum([self.graph.hash_graph[target][tar] for tar in self.graph.hash_graph[target] if
#                            tar != current_seed])
#         add_candidates[target] = Relationship(target_weight_to,
#                                               target_num_edges_to,
#                                               target_weight_from,
#                                               target_num_edges_from)
#     return current_cluster,   \
#            remove_candidates, \
#            add_candidates,    \
#            current_score,     \
#            current_cluster_weight_in, \
#            current_cluster_weight_out
#
#
# # maybe if we add vertex that is best in terms of 2 neighborhood
# #       but may or may not be the best in terms of 1 neighborhood
# # after we add this vertex, we add the helpful neighbor vertices, without considering the whole boundary
# def find_best_2neighborhood_add(self, add_candidates, current_cluster, current_score, current_cluster_weight_in, current_cluster_weight_out):
#     best_change = None
#     best_change_score = current_score
#     factor = .3
#     for v in add_candidates:
#         numerator = current_cluster_weight_in + add_candidates[v]._in
#         denominator = current_cluster_weight_in + current_cluster_weight_out + add_candidates[
#             v]._out + self.penalty_value_per_node * (len(current_cluster) + 1)
#         #####################
#         # CONSIDER WHAT COULD BE GAINED ON THE NEWLY EXPOSED BOUNDARY
#         #####################
#         # find the neighbors of v
#         distance_2_neighbors = self.graph.hash_graph[v]
#         for d2n in distance_2_neighbors:
#             if d2n in current_cluster:
#                 numerator += factor * self.graph.hash_graph[v][d2n]
#             else:
#                 distance_3_neighbors = self.graph.hash_graph[d2n]
#                 for d3n in distance_3_neighbors:
#                     if d3n in current_cluster:
#                         numerator += factor * self.graph.hash_graph[d2n][d3n]
#                     # elif d3n in add_candidates:
#                     #     numerator += self.graph.hash_graph[d2n][d3n]
#                     else:
#                         denominator += factor * self.graph.hash_graph[d2n][d3n]
#
#         proposed_score = numerator / denominator
#         if proposed_score > best_change_score:
#             best_change = v
#             best_change_score = proposed_score
#         debug("##################### ADD Consideration ########################")
#         debug("v: %s" % str(v))
#         debug("proposed_score: %s" % str(proposed_score))
#         debug("best_change_score: %s" % str(best_change_score))
#         debug("best_change: %s" % str(best_change))
#         debug("numerator: %s" % str(numerator))
#         debug("denominator: %s" % str(denominator))
#         debug("current_cluster_weight_in: %s" % str(current_cluster_weight_in))
#         debug("add_candidates[v]._in: %s" % str(add_candidates[v]._in))
#         debug("current_cluster_weight_out: %s" % str(current_cluster_weight_out))
#         debug("add_candidates[v]._out: %s" % str(add_candidates[v]._out))
#         debug("len(current_cluster): %s" % str(len(current_cluster)))
#         sleep_debug(.25)
#     return best_change, best_change_score
#
#
# def find_best_add(self, add_candidates, current_cluster, current_score, current_cluster_weight_in, current_cluster_weight_out):
#     best_change = None
#     best_change_score = current_score
#     # best_change_list = []
#     for v in add_candidates:
#         numerator = current_cluster_weight_in +\
#                     add_candidates[v]._in
#         denominator = current_cluster_weight_in+\
#                       current_cluster_weight_out +\
#                       add_candidates[v]._out +\
#                       self.penalty_value_per_node * (len(current_cluster) + 1)
#         proposed_score = numerator / denominator
#         if proposed_score > best_change_score:
#             best_change = v
#             best_change_score = proposed_score
#         # UNFINISHED BEST_CHANGE_LIST
#         #     best_change_list = [v]
#         # elif proposed_score == best_change_score:
#         #     best_change_list.append(v)
#         debug("##################### ADD Consideration ########################")
#         debug("v: %s" % str(v))
#         debug("proposed_score: %s" % str(proposed_score))
#         debug("best_change_score: %s" % str(best_change_score))
#         debug("best_change: %s" % str(best_change))
#         debug("numerator: %s" % str(numerator))
#         debug("denominator: %s" % str(denominator))
#         debug("current_cluster_weight_in: %s" % str(current_cluster_weight_in))
#         debug("add_candidates[v]._in: %s" % str(add_candidates[v]._in))
#         debug("current_cluster_weight_out: %s" % str(current_cluster_weight_out))
#         debug("add_candidates[v]._out: %s" % str(add_candidates[v]._out))
#         debug("len(current_cluster): %s" % str(len(current_cluster)))
#         sleep_debug(.25)
#     # UNFINISHED BEST_CHANGE_LIST
#     # if len( best_change_list)>1:
#     #     print(best_change_list)
#     #     for tie in best_change_list:
#     #         print(add_candidates[tie].stringify())
#     return best_change, best_change_score
#
#
# def add(self, add_candidates, current_cluster, remove_candidates,
#         change_vertex, change_vertex_score, cc_weight_in, cc_weight_out):
#     debug("\n", "ADD: %s" % str(change_vertex), "change_vertex_score: %s" % str(change_vertex_score), "\n")
#     #################################################################
#     # update the overall weight into and out of the current_cluster #
#     #################################################################
#     change_vertex_in = add_candidates[change_vertex]._in
#     change_vertex_out = add_candidates[change_vertex]._out
#     cc_weight_in += change_vertex_in
#     cc_weight_out = cc_weight_out - change_vertex_in + change_vertex_out
#
#     ###################################
#     # Move the change vertex from add_candidates to current_cluster
#     ###################################
#     to_add = add_candidates[change_vertex].copy()
#     del add_candidates[change_vertex]
#     current_cluster[change_vertex] = to_add.copy()
#
#     ###################################
#     # Change vertex to remove_candidates if applicable
#     ###################################
#     if to_add.num_edges_from:
#         remove_candidates[change_vertex] = to_add.copy()
#
#     def update_v(v, edge_weight, collection):
#         collection[v]._in += edge_weight
#         collection[v]._out -= edge_weight
#         collection[v].num_edges_to += 1
#         collection[v].num_edges_from -= 1
#         #######################
#         # sanity check
#         #######################
#         thresh = -.001
#         a = collection[v]._in < thresh
#         b = collection[v]._out < thresh
#         c = collection[v].num_edges_to < thresh
#         d = collection[v].num_edges_from < thresh
#         if a or b or c or d:
#             print("oh no %s; %s%s%s%s"%(str(v), a, b, c, d))
#             exit()
#
#     def initialize_new_add_candidate_for_V(add_vertex, current_cluster, add_candidates):
#         num_edges_to = 0
#         weight_to = 0
#         num_edges_from = 0
#         weight_from = 0
#         ###################################
#         # iterate over the neighbors of v
#         ###################################
#         for neighbor in self.graph.hash_graph[add_vertex]:
#             weight_prime = self.graph.hash_graph[add_vertex][neighbor]
#             if neighbor in current_cluster:
#                 num_edges_to += 1
#                 weight_to += weight_prime
#             else:
#                 num_edges_from += 1
#                 weight_from += weight_prime
#
#         add_candidates[add_vertex] = Relationship(weight_to,
#                                          num_edges_to,
#                                          weight_from,
#                                          num_edges_from)
#
#     #######################################################################
#     # iterate over neighbors of change_vertex, and update each Relationship
#     #######################################################################
#     for v in self.graph.hash_graph[change_vertex]:
#         edge_weight = self.graph.hash_graph[v][change_vertex]
#         if v in add_candidates:
#             update_v(v, edge_weight, add_candidates)
#         if v in current_cluster:
#             update_v(v, edge_weight, current_cluster)
#         # note that v may be in both the current_cluster and in remove_candidates
#         # remove_candidates is a subset of current_cluster
#         if v in remove_candidates:
#             update_v(v, edge_weight, remove_candidates)
#             # Check that a candidate for removal is still on the boundary
#             if remove_candidates[v].num_edges_from == 0:
#                 del remove_candidates[v]
#         # handle the case that v is on the new boundary
#         # add v to add_candidates
#         if v not in add_candidates and v not in current_cluster:
#             initialize_new_add_candidate_for_V(v, current_cluster, add_candidates)
#
#     return cc_weight_in, cc_weight_out
#
#
# def find_best_remove(self, remove_candidates, current_cluster, current_cluster_weight_in, current_cluster_weight_out, current_score):
#     best_change = None
#     best_change_score = current_score
#     # check that
#     #   (1) the cluster has more than one element
#     if len(current_cluster) > 1:
#         current_cluster_membership_hashset = [vertex for vertex in current_cluster]
#         for v in remove_candidates:
#             if self.care_about_cuts:
#                 # TODO: check if there is a cut.
#                 #   Implement more efficiently using a Dynamic Connectivity algorithm
#                 is_a_cut = True
#                 visted = set()
#                 start_point = None
#                 for potential_start_point in current_cluster_membership_hashset:
#                     if potential_start_point != v:
#                         start_point = potential_start_point
#                         # consider break statement here
#                 # check that
#                 #   (2) removal of vertex under consideration will not disconnect cluster
#                 debug("DFS starting vertex: %s" % str(start_point), "DFS ignore vertex: %s" % str(v))
#                 dfs(self, start_point, v, current_cluster_membership_hashset, visted)
#                 debug("=============================")
#                 if len(visted) == -1 + len(current_cluster_membership_hashset):
#                     is_a_cut = False
#                     debug("%s is NOT a CUT" % str(v))
#                 if is_a_cut:
#                     debug("%s is a CUT!" % str(v))
#                 debug("cluster: %s" % str(current_cluster_membership_hashset))
#                 debug("visited by DFS: %s" % str(visted))
#                 sleep_debug(.25)
#             else:
#                 is_a_cut = False
#
#             if not is_a_cut:
#                 # TODO: check that this makes sense
#                 numerator = current_cluster_weight_in - remove_candidates[v]._in
#                 denominator = current_cluster_weight_in + current_cluster_weight_out - \
#                               remove_candidates[v]._out + self.penalty_value_per_node * (
#                                       len(current_cluster) - 1)
#                 proposed_score = numerator / denominator
#                 if proposed_score > best_change_score:
#                     best_change = v
#                     best_change_score = proposed_score
#                 debug("##################### REMOVE Consideration ########################")
#                 debug("v: %s" % str(v))
#                 debug("proposed_score: %s" % str(proposed_score))
#                 debug("best_change_score: %s" % str(best_change_score))
#                 debug("best_change: %s" % str(best_change))
#                 debug("numerator: %s" % str(numerator))
#                 debug("denominator: %s" % str(denominator))
#                 debug("current_cluster_weight_in: %s" % str(current_cluster_weight_in))
#                 debug(
#                     "remove_candidates[v]._in: %s" % str(remove_candidates[v]._in))
#                 debug("current_cluster_weight_out: %s" % str(current_cluster_weight_out))
#                 debug("remove_candidates[v]._out: %s" % str(
#                     remove_candidates[v]._out))
#                 debug("len(current_cluster): %s" % str(len(current_cluster)))
#                 sleep_debug(1)
#     return best_change, best_change_score
#
#
# def remove(self,
#            remove_candidates,
#            add_candidates,
#            current_cluster,
#            change_vertex,
#            change_vertex_score,
#            cc_weight_in,
#            cc_weight_out):
#     debug("REMOVE: ", change_vertex, "change_vertex_score: ", change_vertex_score)
#     sleep_debug(1)
#     ###############################################################################################
#     # Update the current_cluster 's score, and overall weight into and out of the current cluster #
#     ###############################################################################################
#     change_vertex_in = remove_candidates[change_vertex]._in
#     change_vertex_out = remove_candidates[change_vertex]._out
#     cc_weight_in -= change_vertex_in
#     cc_weight_out = cc_weight_out - change_vertex_out + change_vertex_in
#
#     #######################################################################
#     # Remove the change vertex from remove_candidates and current_cluster #
#     #######################################################################
#     to_remove = remove_candidates[change_vertex].copy()
#     del remove_candidates[change_vertex]
#     del current_cluster[change_vertex]
#
#     ################################################
#     # Also add the change vertex to add_candidates #
#     ################################################
#     add_candidates[change_vertex] = to_remove
#
#     #########################################################################
#     # iterate over neighbors of change_vertex, and update each Relationship #
#     #########################################################################
#     def update_v(v, edge_weight, collection):
#         collection[v]._in -= edge_weight
#         collection[v]._out += edge_weight
#         collection[v].num_edges_to -= 1
#         collection[v].num_edges_from += 1
#         #######################
#         # sanity check
#         #######################
#         thresh = -.001
#         a = collection[v]._in < thresh
#         b = collection[v]._out < thresh
#         c = collection[v].num_edges_to < thresh
#         d = collection[v].num_edges_from < thresh
#         if a or b or c or d:
#             print("oh no %s; %s%s%s%s"%(str(v), a, b, c, d))
#             exit()
#
#     for v in self.graph.hash_graph[change_vertex]:
#         edge_weight = self.graph.hash_graph[change_vertex][v]
#         # note that v may be in both the current_cluster and in remove_candidates
#         if v in remove_candidates:
#             update_v(v, edge_weight, remove_candidates)
#         if v in current_cluster:
#             update_v(v, edge_weight, current_cluster)
#             if current_cluster[v].num_edges_from == 1:
#                 remove_candidates[v] = current_cluster[v].copy()
#         if v in add_candidates:
#             update_v(v, edge_weight, add_candidates)
#             if add_candidates[v].num_edges_to == 0:
#                 del add_candidates[v]
#     return cc_weight_in, cc_weight_out
#
#
# def find_best_suboptimal_add(self, add_candidates, current_cluster, cc_weight_in, cc_weight_out):
#     best_suboptimal_change = None
#     best_suboptimal_change_score = -10000
#     for v in add_candidates:
#         numerator = cc_weight_in + add_candidates[v]._in
#         denominator = cc_weight_in + cc_weight_out + add_candidates[
#             v]._out + self.penalty_value_per_node * (len(current_cluster) + 1)
#         proposed_score = numerator / denominator
#         if proposed_score > best_suboptimal_change_score:
#             best_suboptimal_change = v
#             best_suboptimal_change_score = proposed_score
#         debug("##################### SUBOPTIMAL ADD Consideration ########################")
#         debug("v: %s" % str(v))
#         debug("proposed_score: %s" % str(proposed_score))
#         debug("best_suboptimal_change_score: %s" % str(best_suboptimal_change_score))
#         debug("best_suboptimal_change: %s" % str(best_suboptimal_change))
#         debug("numerator: %s" % str(numerator))
#         debug("denominator: %s" % str(denominator))
#         debug("current_cluster_weight_in: %s" % str(cc_weight_in))
#         debug("add_candidates[v]._in: %s" % str(add_candidates[v]._in))
#         debug("current_cluster_weight_out: %s" % str(cc_weight_out))
#         debug("add_candidates[v]._out: %s" % str(add_candidates[v]._out))
#         debug("len(current_cluster): %s" % str(len(current_cluster)))
#         sleep_debug(.25)
#     return best_suboptimal_change, best_suboptimal_change_score
#
#
# def add_shake(self, add_candidates, remove_candidates, current_cluster, cc_weight_in, cc_weight_out, round_no, last_failed_add_round_no):
#     for i in range(self.number_of_bad_adds):
#         best_suboptimal_change, best_suboptimal_score = find_best_suboptimal_add(self, add_candidates, current_cluster,
#                                                                                  cc_weight_in, cc_weight_out)
#         if best_suboptimal_change:
#             print("adding_suboptimally")
#             sleep_debug(1)
#             round_no += 1
#             last_failed_add_round_no = -5
#             # TODO handle round_no numbers
#             cc_weight_in, cc_weight_out = add(self, add_candidates, current_cluster, remove_candidates,
#                                               best_suboptimal_change, best_suboptimal_score, cc_weight_in,
#                                               cc_weight_out)
#     return best_suboptimal_score, \
#            cc_weight_in, \
#            cc_weight_out, \
#            round_no, \
#            last_failed_add_round_no


# from graph import *
# from src.common.common import *
#
# def find_best_add_list(self, add_candidates, current_cluster, current_score, current_cluster_weight_in, current_cluster_weight_out):
#     best_change_score = current_score
#     best_change_list = []
#     for v in add_candidates:
#         numerator = current_cluster_weight_in + add_candidates[v]._in
#         denominator = current_cluster_weight_in + current_cluster_weight_out + add_candidates[
#             v]._out + self.penalty_value_per_node * (len(current_cluster) + 1)
#         proposed_score = numerator / denominator
#         if proposed_score > best_change_score:
#             best_change = v
#             best_change_score = proposed_score
#             best_change_list.append((proposed_score, best_change))
#
#     def filter_list(best_change_list, length_to_keep):
#         for index, tup in enumerate(best_change_list):
#             new_tuple = (-1 * tup[0], tup[1])
#             best_change_list[index] = new_tuple
#
#         length_to_keep = min(len(best_change_list), length_to_keep)
#         heapq.heapify(best_change_list)
#         update_list = []
#         for i in range(length_to_keep):
#             update_list.append(heapq.heappop(best_change_list))
#         best_change_list = update_list
#
#     print(best_change_list)
#     filter_list(best_change_list, 7)
#     retval = [t[1] for t in best_change_list]
#     return retval
#
#
# def careful_find_best_2neighborhood_add(self, add_candidates, current_cluster, current_score, current_cluster_weight_in, current_cluster_weight_out):
#     # if len(current_cluster)<5:
#     #     return find_best_add(self, add_candidates, current_cluster, current_score, current_cluster_weight_in, current_cluster_weight_out)
#     # best_change_list = find_best_add_list(self, add_candidates, current_cluster, current_score, current_cluster_weight_in, current_cluster_weight_out)
#
#     best_change = None
#     best_change_score = current_score
#     best_proposed_score = current_score
#     neighborhood_gain = {v: {"in": 0,
#                             "out": 0}
#                                         for v in add_candidates}
#     good_neighbors = dict()
#     for v in add_candidates:
#     # for v in best_change_list:
#         numerator = current_cluster_weight_in + \
#                     add_candidates[v]._in
#         denominator = current_cluster_weight_in + \
#                       current_cluster_weight_out + \
#                       self.penalty_value_per_node * (len(current_cluster) + 1)+\
#                       add_candidates[v]._out
#
#         actual_score = numerator/denominator
#         #################################################################
#         #  CONSIDER WHAT COULD BE GAINED ON THE NEWLY EXPOSED BOUNDARY  #
#         #################################################################
#         distance_2_neighbors = self.graph.hash_graph[v]
#         in_and_out_for_distance_2_neighbors = dict()
#         for d2n in distance_2_neighbors:
#             if d2n in current_cluster:
#                 #########################################
#                 #  THIS HAS ALREADY BEEN ACCOUNTED FOR  #
#                 #########################################
#                 '''
#                 already captured in the term current_cluster_weight_in
#                 '''
#                 continue
#             else:
#                 in_and_out_for_distance_2_neighbors[d2n]={"to_current_cluster":0,
#                                   "to_v":0,
#                                   "to_N(v)": 0,
#                                   "out":0}
#                 distance_3_neighbors = self.graph.hash_graph[d2n]
#                 for d3n in distance_3_neighbors:
#                     edge_weight = self.graph.hash_graph[d2n][d3n]
#                     if d3n is v:
#                         in_and_out_for_distance_2_neighbors[d2n]["to_v"] += 3*edge_weight
#                     if d3n in current_cluster:
#                         in_and_out_for_distance_2_neighbors[d2n]["to_current_cluster"] += 3 * edge_weight
#                     elif d3n in distance_2_neighbors:
#                         in_and_out_for_distance_2_neighbors[d2n]["to_N(v)"] = edge_weight
#                     else:
#                         in_and_out_for_distance_2_neighbors[d2n]["out"] = edge_weight
#         d2n_to_grab = set()
#         for d2n in in_and_out_for_distance_2_neighbors:
#             # TODO: figure out what to do with "to_N(v)"
#             inbound = in_and_out_for_distance_2_neighbors[d2n]["to_current_cluster"] + \
#                       in_and_out_for_distance_2_neighbors[d2n]["to_v"]
#             outbound = in_and_out_for_distance_2_neighbors[d2n]["out"]
#             if inbound>=outbound:
#             # if inbound/(inbound+outbound)>actual_score:
#                 d2n_to_grab.add(d2n)
#                 neighborhood_gain[v]["in"]+=inbound
#                 neighborhood_gain[v]["out"]+=outbound
#         good_neighbors[v] = d2n_to_grab
#         # in_N_v_grabbed = 0
#         # for a in d2n_to_grab:
#         #     for b in d2n_to_grab:
#         #         if b in self.graph.hash_graph[a]:
#         #             in_N_v_grabbed += self.graph.hash_graph[a][b]
#         # dict_gain[v]["in"]+=in_N_v_grabbed
#         # print(v, d2n_to_grab)
#         # factor = len(current_cluster)*.09
#         # factor = logistic.cdf(len(current_cluster)/10)
#         # factor = .2* logistic.cdf(len(current_cluster))
#         factor = len(current_cluster)/10
#         numerator+= factor*neighborhood_gain[v]["in"]
#         denominator+= factor*\
#                       (neighborhood_gain[v]["in"] + neighborhood_gain[v]["out"])
#
#         proposed_score = numerator / denominator
#         if proposed_score > best_proposed_score:
#             best_change = v
#             best_proposed_score = proposed_score
#             best_change_score = actual_score
#         sleep_debug(.25)
#     # print("-------------",best_change)
#     return best_change, best_change_score, #good_neighbors[best_change]
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
# def find_best_2neighborhood_add(self, add_candidates, current_cluster, current_score, current_cluster_weight_in, current_cluster_weight_out):
#     best_change = None
#     best_change_score = current_score
#     factor = .3
#     for v in add_candidates:
#         numerator = current_cluster_weight_in + add_candidates[v]._in
#         denominator = current_cluster_weight_in + current_cluster_weight_out + add_candidates[
#             v]._out + self.penalty_value_per_node * (len(current_cluster) + 1)
#         #####################
#         # CONSIDER WHAT COULD BE GAINED ON THE NEWLY EXPOSED BOUNDARY
#         #####################
#         # find the neighbors of v
#         distance_2_neighbors = self.graph.hash_graph[v]
#         for d2n in distance_2_neighbors:
#             if d2n in current_cluster:
#                 numerator += factor * self.graph.hash_graph[v][d2n]
#             else:
#                 distance_3_neighbors = self.graph.hash_graph[d2n]
#                 for d3n in distance_3_neighbors:
#                     if d3n in current_cluster:
#                         numerator += factor * self.graph.hash_graph[d2n][d3n]
#                     # elif d3n in add_candidates:
#                     #     numerator += self.graph.hash_graph[d2n][d3n]
#                     else:
#                         denominator += factor * self.graph.hash_graph[d2n][d3n]
#
#         proposed_score = numerator / denominator
#         if proposed_score > best_change_score:
#             best_change = v
#             best_change_score = proposed_score
#         debug("##################### ADD Consideration ########################")
#         debug("v: %s" % str(v))
#         debug("proposed_score: %s" % str(proposed_score))
#         debug("best_change_score: %s" % str(best_change_score))
#         debug("best_change: %s" % str(best_change))
#         debug("numerator: %s" % str(numerator))
#         debug("denominator: %s" % str(denominator))
#         debug("current_cluster_weight_in: %s" % str(current_cluster_weight_in))
#         debug("add_candidates[v]._in: %s" % str(add_candidates[v]._in))
#         debug("current_cluster_weight_out: %s" % str(current_cluster_weight_out))
#         debug("add_candidates[v]._out: %s" % str(add_candidates[v]._out))
#         debug("len(current_cluster): %s" % str(len(current_cluster)))
#         sleep_debug(.25)
#     return best_change, best_change_score
#
#
#
#
#
#
#
#
#
# def initialize_complex(self, current_seed):
#     ###########################################
#     # initalize the current cluster
#     ###########################################
#     current_cluster = dict()
#     seed_weight_to = 0
#     seed_num_edges_to = 0
#     seed_weight_from = sum([self.graph.hash_graph[current_seed][tar] for tar in self.graph.hash_graph[current_seed]])
#     seed_num_edges_from = len(self.graph.hash_graph[current_seed])
#     current_cluster[current_seed] = Relationship(seed_weight_to,
#                                                  seed_num_edges_to,
#                                                  seed_weight_from,
#                                                  seed_num_edges_from)
#     current_score = 0
#     current_cluster_weight_in = 0
#     current_cluster_weight_out = seed_weight_from
#
#     ###########################################
#     # initalize the candidates for removal
#     ###########################################
#     remove_candidates = dict()
#     remove_candidates[current_seed] = current_cluster[current_seed].copy()
#
#     ###########################################
#     # initialize the candidates for addition
#     ###########################################
#     add_candidates = dict()
#     for target in self.graph.hash_graph[current_seed]:
#         target_weight_to = self.graph.hash_graph[current_seed][target]
#         target_num_edges = len(self.graph.hash_graph[target])
#         target_num_edges_to = 1
#         target_num_edges_from = target_num_edges - 1
#         target_weight_from = sum([self.graph.hash_graph[target][tar] for tar in self.graph.hash_graph[target] if
#                            tar != current_seed])
#         add_candidates[target] = Relationship(target_weight_to,
#                                               target_num_edges_to,
#                                               target_weight_from,
#                                               target_num_edges_from)
#     return current_cluster,   \
#            remove_candidates, \
#            add_candidates,    \
#            current_score,     \
#            current_cluster_weight_in, \
#            current_cluster_weight_out
#
#
# # maybe if we add vertex that is best in terms of 2 neighborhood
# #       but may or may not be the best in terms of 1 neighborhood
# # after we add this vertex, we add the helpful neighbor vertices, without considering the whole boundary
# def find_best_2neighborhood_add(self, add_candidates, current_cluster, current_score, current_cluster_weight_in, current_cluster_weight_out):
#     best_change = None
#     best_change_score = current_score
#     factor = .3
#     for v in add_candidates:
#         numerator = current_cluster_weight_in + add_candidates[v]._in
#         denominator = current_cluster_weight_in + current_cluster_weight_out + add_candidates[
#             v]._out + self.penalty_value_per_node * (len(current_cluster) + 1)
#         #####################
#         # CONSIDER WHAT COULD BE GAINED ON THE NEWLY EXPOSED BOUNDARY
#         #####################
#         # find the neighbors of v
#         distance_2_neighbors = self.graph.hash_graph[v]
#         for d2n in distance_2_neighbors:
#             if d2n in current_cluster:
#                 numerator += factor * self.graph.hash_graph[v][d2n]
#             else:
#                 distance_3_neighbors = self.graph.hash_graph[d2n]
#                 for d3n in distance_3_neighbors:
#                     if d3n in current_cluster:
#                         numerator += factor * self.graph.hash_graph[d2n][d3n]
#                     # elif d3n in add_candidates:
#                     #     numerator += self.graph.hash_graph[d2n][d3n]
#                     else:
#                         denominator += factor * self.graph.hash_graph[d2n][d3n]
#
#         proposed_score = numerator / denominator
#         if proposed_score > best_change_score:
#             best_change = v
#             best_change_score = proposed_score
#         debug("##################### ADD Consideration ########################")
#         debug("v: %s" % str(v))
#         debug("proposed_score: %s" % str(proposed_score))
#         debug("best_change_score: %s" % str(best_change_score))
#         debug("best_change: %s" % str(best_change))
#         debug("numerator: %s" % str(numerator))
#         debug("denominator: %s" % str(denominator))
#         debug("current_cluster_weight_in: %s" % str(current_cluster_weight_in))
#         debug("add_candidates[v]._in: %s" % str(add_candidates[v]._in))
#         debug("current_cluster_weight_out: %s" % str(current_cluster_weight_out))
#         debug("add_candidates[v]._out: %s" % str(add_candidates[v]._out))
#         debug("len(current_cluster): %s" % str(len(current_cluster)))
#         sleep_debug(.25)
#     return best_change, best_change_score
#
#
# def find_best_add(self, add_candidates, current_cluster, current_score, current_cluster_weight_in, current_cluster_weight_out):
#     best_change = None
#     best_change_score = current_score
#     # best_change_list = []
#     for v in add_candidates:
#         numerator = current_cluster_weight_in +\
#                     add_candidates[v]._in
#         denominator = current_cluster_weight_in+\
#                       current_cluster_weight_out +\
#                       add_candidates[v]._out +\
#                       self.penalty_value_per_node * (len(current_cluster) + 1)
#         proposed_score = numerator / denominator
#         if proposed_score > best_change_score:
#             best_change = v
#             best_change_score = proposed_score
#         # UNFINISHED BEST_CHANGE_LIST
#         #     best_change_list = [v]
#         # elif proposed_score == best_change_score:
#         #     best_change_list.append(v)
#         debug("##################### ADD Consideration ########################")
#         debug("v: %s" % str(v))
#         debug("proposed_score: %s" % str(proposed_score))
#         debug("best_change_score: %s" % str(best_change_score))
#         debug("best_change: %s" % str(best_change))
#         debug("numerator: %s" % str(numerator))
#         debug("denominator: %s" % str(denominator))
#         debug("current_cluster_weight_in: %s" % str(current_cluster_weight_in))
#         debug("add_candidates[v]._in: %s" % str(add_candidates[v]._in))
#         debug("current_cluster_weight_out: %s" % str(current_cluster_weight_out))
#         debug("add_candidates[v]._out: %s" % str(add_candidates[v]._out))
#         debug("len(current_cluster): %s" % str(len(current_cluster)))
#         sleep_debug(.25)
#     # UNFINISHED BEST_CHANGE_LIST
#     # if len( best_change_list)>1:
#     #     print(best_change_list)
#     #     for tie in best_change_list:
#     #         print(add_candidates[tie].stringify())
#     return best_change, best_change_score
#
#
# def add(self, add_candidates, current_cluster, remove_candidates,
#         change_vertex, change_vertex_score, cc_weight_in, cc_weight_out):
#     debug("\n", "ADD: %s" % str(change_vertex), "change_vertex_score: %s" % str(change_vertex_score), "\n")
#     #################################################################
#     # update the overall weight into and out of the current_cluster #
#     #################################################################
#     change_vertex_in = add_candidates[change_vertex]._in
#     change_vertex_out = add_candidates[change_vertex]._out
#     cc_weight_in += change_vertex_in
#     cc_weight_out = cc_weight_out - change_vertex_in + change_vertex_out
#
#     ###################################
#     # Move the change vertex from add_candidates to current_cluster
#     ###################################
#     to_add = add_candidates[change_vertex].copy()
#     del add_candidates[change_vertex]
#     current_cluster[change_vertex] = to_add.copy()
#
#     ###################################
#     # Change vertex to remove_candidates if applicable
#     ###################################
#     if to_add.num_edges_from:
#         remove_candidates[change_vertex] = to_add.copy()
#
#     def update_v(v, edge_weight, collection):
#         collection[v]._in += edge_weight
#         collection[v]._out -= edge_weight
#         collection[v].num_edges_to += 1
#         collection[v].num_edges_from -= 1
#         #######################
#         # sanity check
#         #######################
#         thresh = -.001
#         a = collection[v]._in < thresh
#         b = collection[v]._out < thresh
#         c = collection[v].num_edges_to < thresh
#         d = collection[v].num_edges_from < thresh
#         if a or b or c or d:
#             print("oh no %s; %s%s%s%s"%(str(v), a, b, c, d))
#             exit()
#
#     def initialize_new_add_candidate_for_V(add_vertex, current_cluster, add_candidates):
#         num_edges_to = 0
#         weight_to = 0
#         num_edges_from = 0
#         weight_from = 0
#         ###################################
#         # iterate over the neighbors of v
#         ###################################
#         for neighbor in self.graph.hash_graph[add_vertex]:
#             weight_prime = self.graph.hash_graph[add_vertex][neighbor]
#             if neighbor in current_cluster:
#                 num_edges_to += 1
#                 weight_to += weight_prime
#             else:
#                 num_edges_from += 1
#                 weight_from += weight_prime
#
#         add_candidates[add_vertex] = Relationship(weight_to,
#                                          num_edges_to,
#                                          weight_from,
#                                          num_edges_from)
#
#     #######################################################################
#     # iterate over neighbors of change_vertex, and update each Relationship
#     #######################################################################
#     for v in self.graph.hash_graph[change_vertex]:
#         edge_weight = self.graph.hash_graph[v][change_vertex]
#         if v in add_candidates:
#             update_v(v, edge_weight, add_candidates)
#         if v in current_cluster:
#             update_v(v, edge_weight, current_cluster)
#         # note that v may be in both the current_cluster and in remove_candidates
#         # remove_candidates is a subset of current_cluster
#         if v in remove_candidates:
#             update_v(v, edge_weight, remove_candidates)
#             # Check that a candidate for removal is still on the boundary
#             if remove_candidates[v].num_edges_from == 0:
#                 del remove_candidates[v]
#         # handle the case that v is on the new boundary
#         # add v to add_candidates
#         if v not in add_candidates and v not in current_cluster:
#             initialize_new_add_candidate_for_V(v, current_cluster, add_candidates)
#
#     return cc_weight_in, cc_weight_out
#
#
# def find_best_remove(self, remove_candidates, current_cluster, current_cluster_weight_in, current_cluster_weight_out, current_score):
#     best_change = None
#     best_change_score = current_score
#     # check that
#     #   (1) the cluster has more than one element
#     if len(current_cluster) > 1:
#         current_cluster_membership_hashset = [vertex for vertex in current_cluster]
#         for v in remove_candidates:
#             if self.care_about_cuts:
#                 # TODO: check if there is a cut.
#                 #   Implement more efficiently using a Dynamic Connectivity algorithm
#                 is_a_cut = True
#                 visted = set()
#                 start_point = None
#                 for potential_start_point in current_cluster_membership_hashset:
#                     if potential_start_point != v:
#                         start_point = potential_start_point
#                         # consider break statement here
#                 # check that
#                 #   (2) removal of vertex under consideration will not disconnect cluster
#                 debug("DFS starting vertex: %s" % str(start_point), "DFS ignore vertex: %s" % str(v))
#                 dfs(self, start_point, v, current_cluster_membership_hashset, visted)
#                 debug("=============================")
#                 if len(visted) == -1 + len(current_cluster_membership_hashset):
#                     is_a_cut = False
#                     debug("%s is NOT a CUT" % str(v))
#                 if is_a_cut:
#                     debug("%s is a CUT!" % str(v))
#                 debug("cluster: %s" % str(current_cluster_membership_hashset))
#                 debug("visited by DFS: %s" % str(visted))
#                 sleep_debug(.25)
#             else:
#                 is_a_cut = False
#
#             if not is_a_cut:
#                 # TODO: check that this makes sense
#                 numerator = current_cluster_weight_in - remove_candidates[v]._in
#                 denominator = current_cluster_weight_in + current_cluster_weight_out - \
#                               remove_candidates[v]._out + self.penalty_value_per_node * (
#                                       len(current_cluster) - 1)
#                 proposed_score = numerator / denominator
#                 if proposed_score > best_change_score:
#                     best_change = v
#                     best_change_score = proposed_score
#                 debug("##################### REMOVE Consideration ########################")
#                 debug("v: %s" % str(v))
#                 debug("proposed_score: %s" % str(proposed_score))
#                 debug("best_change_score: %s" % str(best_change_score))
#                 debug("best_change: %s" % str(best_change))
#                 debug("numerator: %s" % str(numerator))
#                 debug("denominator: %s" % str(denominator))
#                 debug("current_cluster_weight_in: %s" % str(current_cluster_weight_in))
#                 debug(
#                     "remove_candidates[v]._in: %s" % str(remove_candidates[v]._in))
#                 debug("current_cluster_weight_out: %s" % str(current_cluster_weight_out))
#                 debug("remove_candidates[v]._out: %s" % str(
#                     remove_candidates[v]._out))
#                 debug("len(current_cluster): %s" % str(len(current_cluster)))
#                 sleep_debug(1)
#     return best_change, best_change_score
#
#
# def remove(self,
#            remove_candidates,
#            add_candidates,
#            current_cluster,
#            change_vertex,
#            change_vertex_score,
#            cc_weight_in,
#            cc_weight_out):
#     debug("REMOVE: ", change_vertex, "change_vertex_score: ", change_vertex_score)
#     sleep_debug(1)
#     ###############################################################################################
#     # Update the current_cluster 's score, and overall weight into and out of the current cluster #
#     ###############################################################################################
#     change_vertex_in = remove_candidates[change_vertex]._in
#     change_vertex_out = remove_candidates[change_vertex]._out
#     cc_weight_in -= change_vertex_in
#     cc_weight_out = cc_weight_out - change_vertex_out + change_vertex_in
#
#     #######################################################################
#     # Remove the change vertex from remove_candidates and current_cluster #
#     #######################################################################
#     to_remove = remove_candidates[change_vertex].copy()
#     del remove_candidates[change_vertex]
#     del current_cluster[change_vertex]
#
#     ################################################
#     # Also add the change vertex to add_candidates #
#     ################################################
#     add_candidates[change_vertex] = to_remove
#
#     #########################################################################
#     # iterate over neighbors of change_vertex, and update each Relationship #
#     #########################################################################
#     def update_v(v, edge_weight, collection):
#         collection[v]._in -= edge_weight
#         collection[v]._out += edge_weight
#         collection[v].num_edges_to -= 1
#         collection[v].num_edges_from += 1
#         #######################
#         # sanity check
#         #######################
#         thresh = -.001
#         a = collection[v]._in < thresh
#         b = collection[v]._out < thresh
#         c = collection[v].num_edges_to < thresh
#         d = collection[v].num_edges_from < thresh
#         if a or b or c or d:
#             print("oh no %s; %s%s%s%s"%(str(v), a, b, c, d))
#             exit()
#
#     for v in self.graph.hash_graph[change_vertex]:
#         edge_weight = self.graph.hash_graph[change_vertex][v]
#         # note that v may be in both the current_cluster and in remove_candidates
#         if v in remove_candidates:
#             update_v(v, edge_weight, remove_candidates)
#         if v in current_cluster:
#             update_v(v, edge_weight, current_cluster)
#             if current_cluster[v].num_edges_from == 1:
#                 remove_candidates[v] = current_cluster[v].copy()
#         if v in add_candidates:
#             update_v(v, edge_weight, add_candidates)
#             if add_candidates[v].num_edges_to == 0:
#                 del add_candidates[v]
#     return cc_weight_in, cc_weight_out
#
#
# def find_best_suboptimal_add(self, add_candidates, current_cluster, cc_weight_in, cc_weight_out):
#     best_suboptimal_change = None
#     best_suboptimal_change_score = -10000
#     for v in add_candidates:
#         numerator = cc_weight_in + add_candidates[v]._in
#         denominator = cc_weight_in + cc_weight_out + add_candidates[
#             v]._out + self.penalty_value_per_node * (len(current_cluster) + 1)
#         proposed_score = numerator / denominator
#         if proposed_score > best_suboptimal_change_score:
#             best_suboptimal_change = v
#             best_suboptimal_change_score = proposed_score
#         debug("##################### SUBOPTIMAL ADD Consideration ########################")
#         debug("v: %s" % str(v))
#         debug("proposed_score: %s" % str(proposed_score))
#         debug("best_suboptimal_change_score: %s" % str(best_suboptimal_change_score))
#         debug("best_suboptimal_change: %s" % str(best_suboptimal_change))
#         debug("numerator: %s" % str(numerator))
#         debug("denominator: %s" % str(denominator))
#         debug("current_cluster_weight_in: %s" % str(cc_weight_in))
#         debug("add_candidates[v]._in: %s" % str(add_candidates[v]._in))
#         debug("current_cluster_weight_out: %s" % str(cc_weight_out))
#         debug("add_candidates[v]._out: %s" % str(add_candidates[v]._out))
#         debug("len(current_cluster): %s" % str(len(current_cluster)))
#         sleep_debug(.25)
#     return best_suboptimal_change, best_suboptimal_change_score
#
#
# def add_shake(self, add_candidates, remove_candidates, current_cluster, cc_weight_in, cc_weight_out, round_no, last_failed_add_round_no):
#     for i in range(self.number_of_bad_adds):
#         best_suboptimal_change, best_suboptimal_score = find_best_suboptimal_add(self, add_candidates, current_cluster,
#                                                                                  cc_weight_in, cc_weight_out)
#         if best_suboptimal_change:
#             print("adding_suboptimally")
#             sleep_debug(1)
#             round_no += 1
#             last_failed_add_round_no = -5
#             # TODO handle round_no numbers
#             cc_weight_in, cc_weight_out = add(self, add_candidates, current_cluster, remove_candidates,
#                                               best_suboptimal_change, best_suboptimal_score, cc_weight_in,
#                                               cc_weight_out)
#     return best_suboptimal_score, \
#            cc_weight_in, \
#            cc_weight_out, \
#            round_no, \
#            last_failed_add_round_no
#
#
# # from graph import *
# # from src.common.common import *
# #
# # def find_best_add_list(self, add_candidates, current_cluster, current_score, current_cluster_weight_in, current_cluster_weight_out):
# #     best_change_score = current_score
# #     best_change_list = []
# #     for v in add_candidates:
# #         numerator = current_cluster_weight_in + add_candidates[v]._in
# #         denominator = current_cluster_weight_in + current_cluster_weight_out + add_candidates[
# #             v]._out + self.penalty_value_per_node * (len(current_cluster) + 1)
# #         proposed_score = numerator / denominator
# #         if proposed_score > best_change_score:
# #             best_change = v
# #             best_change_score = proposed_score
# #             best_change_list.append((proposed_score, best_change))
# #
# #     def filter_list(best_change_list, length_to_keep):
# #         for index, tup in enumerate(best_change_list):
# #             new_tuple = (-1 * tup[0], tup[1])
# #             best_change_list[index] = new_tuple
# #
# #         length_to_keep = min(len(best_change_list), length_to_keep)
# #         heapq.heapify(best_change_list)
# #         update_list = []
# #         for i in range(length_to_keep):
# #             update_list.append(heapq.heappop(best_change_list))
# #         best_change_list = update_list
# #
# #     print(best_change_list)
# #     filter_list(best_change_list, 7)
# #     retval = [t[1] for t in best_change_list]
# #     return retval
# #
# #
# # def careful_find_best_2neighborhood_add(self, add_candidates, current_cluster, current_score, current_cluster_weight_in, current_cluster_weight_out):
# #     # if len(current_cluster)<5:
# #     #     return find_best_add(self, add_candidates, current_cluster, current_score, current_cluster_weight_in, current_cluster_weight_out)
# #     # best_change_list = find_best_add_list(self, add_candidates, current_cluster, current_score, current_cluster_weight_in, current_cluster_weight_out)
# #
# #     best_change = None
# #     best_change_score = current_score
# #     best_proposed_score = current_score
# #     neighborhood_gain = {v: {"in": 0,
# #                             "out": 0}
# #                                         for v in add_candidates}
# #     good_neighbors = dict()
# #     for v in add_candidates:
# #     # for v in best_change_list:
# #         numerator = current_cluster_weight_in + \
# #                     add_candidates[v]._in
# #         denominator = current_cluster_weight_in + \
# #                       current_cluster_weight_out + \
# #                       self.penalty_value_per_node * (len(current_cluster) + 1)+\
# #                       add_candidates[v]._out
# #
# #         actual_score = numerator/denominator
# #         #################################################################
# #         #  CONSIDER WHAT COULD BE GAINED ON THE NEWLY EXPOSED BOUNDARY  #
# #         #################################################################
# #         distance_2_neighbors = self.graph.hash_graph[v]
# #         in_and_out_for_distance_2_neighbors = dict()
# #         for d2n in distance_2_neighbors:
# #             if d2n in current_cluster:
# #                 #########################################
# #                 #  THIS HAS ALREADY BEEN ACCOUNTED FOR  #
# #                 #########################################
# #                 '''
# #                 already captured in the term current_cluster_weight_in
# #                 '''
# #                 continue
# #             else:
# #                 in_and_out_for_distance_2_neighbors[d2n]={"to_current_cluster":0,
# #                                   "to_v":0,
# #                                   "to_N(v)": 0,
# #                                   "out":0}
# #                 distance_3_neighbors = self.graph.hash_graph[d2n]
# #                 for d3n in distance_3_neighbors:
# #                     edge_weight = self.graph.hash_graph[d2n][d3n]
# #                     if d3n is v:
# #                         in_and_out_for_distance_2_neighbors[d2n]["to_v"] += 3*edge_weight
# #                     if d3n in current_cluster:
# #                         in_and_out_for_distance_2_neighbors[d2n]["to_current_cluster"] += 3 * edge_weight
# #                     elif d3n in distance_2_neighbors:
# #                         in_and_out_for_distance_2_neighbors[d2n]["to_N(v)"] = edge_weight
# #                     else:
# #                         in_and_out_for_distance_2_neighbors[d2n]["out"] = edge_weight
# #         d2n_to_grab = set()
# #         for d2n in in_and_out_for_distance_2_neighbors:
# #             # TODO: figure out what to do with "to_N(v)"
# #             inbound = in_and_out_for_distance_2_neighbors[d2n]["to_current_cluster"] + \
# #                       in_and_out_for_distance_2_neighbors[d2n]["to_v"]
# #             outbound = in_and_out_for_distance_2_neighbors[d2n]["out"]
# #             if inbound>=outbound:
# #             # if inbound/(inbound+outbound)>actual_score:
# #                 d2n_to_grab.add(d2n)
# #                 neighborhood_gain[v]["in"]+=inbound
# #                 neighborhood_gain[v]["out"]+=outbound
# #         good_neighbors[v] = d2n_to_grab
# #         # in_N_v_grabbed = 0
# #         # for a in d2n_to_grab:
# #         #     for b in d2n_to_grab:
# #         #         if b in self.graph.hash_graph[a]:
# #         #             in_N_v_grabbed += self.graph.hash_graph[a][b]
# #         # dict_gain[v]["in"]+=in_N_v_grabbed
# #         # print(v, d2n_to_grab)
# #         # factor = len(current_cluster)*.09
# #         # factor = logistic.cdf(len(current_cluster)/10)
# #         # factor = .2* logistic.cdf(len(current_cluster))
# #         factor = len(current_cluster)/10
# #         numerator+= factor*neighborhood_gain[v]["in"]
# #         denominator+= factor*\
# #                       (neighborhood_gain[v]["in"] + neighborhood_gain[v]["out"])
# #
# #         proposed_score = numerator / denominator
# #         if proposed_score > best_proposed_score:
# #             best_change = v
# #             best_proposed_score = proposed_score
# #             best_change_score = actual_score
# #         sleep_debug(.25)
# #     # print("-------------",best_change)
# #     return best_change, best_change_score, #good_neighbors[best_change]
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# # def find_best_2neighborhood_add(self, add_candidates, current_cluster, current_score, current_cluster_weight_in, current_cluster_weight_out):
# #     best_change = None
# #     best_change_score = current_score
# #     factor = .3
# #     for v in add_candidates:
# #         numerator = current_cluster_weight_in + add_candidates[v]._in
# #         denominator = current_cluster_weight_in + current_cluster_weight_out + add_candidates[
# #             v]._out + self.penalty_value_per_node * (len(current_cluster) + 1)
# #         #####################
# #         # CONSIDER WHAT COULD BE GAINED ON THE NEWLY EXPOSED BOUNDARY
# #         #####################
# #         # find the neighbors of v
# #         distance_2_neighbors = self.graph.hash_graph[v]
# #         for d2n in distance_2_neighbors:
# #             if d2n in current_cluster:
# #                 numerator += factor * self.graph.hash_graph[v][d2n]
# #             else:
# #                 distance_3_neighbors = self.graph.hash_graph[d2n]
# #                 for d3n in distance_3_neighbors:
# #                     if d3n in current_cluster:
# #                         numerator += factor * self.graph.hash_graph[d2n][d3n]
# #                     # elif d3n in add_candidates:
# #                     #     numerator += self.graph.hash_graph[d2n][d3n]
# #                     else:
# #                         denominator += factor * self.graph.hash_graph[d2n][d3n]
# #
# #         proposed_score = numerator / denominator
# #         if proposed_score > best_change_score:
# #             best_change = v
# #             best_change_score = proposed_score
# #         debug("##################### ADD Consideration ########################")
# #         debug("v: %s" % str(v))
# #         debug("proposed_score: %s" % str(proposed_score))
# #         debug("best_change_score: %s" % str(best_change_score))
# #         debug("best_change: %s" % str(best_change))
# #         debug("numerator: %s" % str(numerator))
# #         debug("denominator: %s" % str(denominator))
# #         debug("current_cluster_weight_in: %s" % str(current_cluster_weight_in))
# #         debug("add_candidates[v]._in: %s" % str(add_candidates[v]._in))
# #         debug("current_cluster_weight_out: %s" % str(current_cluster_weight_out))
# #         debug("add_candidates[v]._out: %s" % str(add_candidates[v]._out))
# #         debug("len(current_cluster): %s" % str(len(current_cluster)))
# #         sleep_debug(.25)
# #     return best_change, best_change_score
# #
# #
# #
# #
# #
# #
# #
# #
# #
# # def initialize_complex(self, current_seed):
# #     ###########################################
# #     # initalize the current cluster
# #     ###########################################
# #     current_cluster = dict()
# #     seed_weight_to = 0
# #     seed_num_edges_to = 0
# #     seed_weight_from = sum([self.graph.hash_graph[current_seed][tar] for tar in self.graph.hash_graph[current_seed]])
# #     seed_num_edges_from = len(self.graph.hash_graph[current_seed])
# #     current_cluster[current_seed] = Relationship(seed_weight_to,
# #                                                  seed_num_edges_to,
# #                                                  seed_weight_from,
# #                                                  seed_num_edges_from)
# #     current_score = 0
# #     current_cluster_weight_in = 0
# #     current_cluster_weight_out = seed_weight_from
# #
# #     ###########################################
# #     # initalize the candidates for removal
# #     ###########################################
# #     remove_candidates = dict()
# #     remove_candidates[current_seed] = current_cluster[current_seed].copy()
# #
# #     ###########################################
# #     # initialize the candidates for addition
# #     ###########################################
# #     add_candidates = dict()
# #     for target in self.graph.hash_graph[current_seed]:
# #         target_weight_to = self.graph.hash_graph[current_seed][target]
# #         target_num_edges = len(self.graph.hash_graph[target])
# #         target_num_edges_to = 1
# #         target_num_edges_from = target_num_edges - 1
# #         target_weight_from = sum([self.graph.hash_graph[target][tar] for tar in self.graph.hash_graph[target] if
# #                            tar != current_seed])
# #         add_candidates[target] = Relationship(target_weight_to,
# #                                               target_num_edges_to,
# #                                               target_weight_from,
# #                                               target_num_edges_from)
# #     return current_cluster,   \
# #            remove_candidates, \
# #            add_candidates,    \
# #            current_score,     \
# #            current_cluster_weight_in, \
# #            current_cluster_weight_out
# #
# #
# # # maybe if we add vertex that is best in terms of 2 neighborhood
# # #       but may or may not be the best in terms of 1 neighborhood
# # # after we add this vertex, we add the helpful neighbor vertices, without considering the whole boundary
# # def find_best_2neighborhood_add(self, add_candidates, current_cluster, current_score, current_cluster_weight_in, current_cluster_weight_out):
# #     best_change = None
# #     best_change_score = current_score
# #     factor = .3
# #     for v in add_candidates:
# #         numerator = current_cluster_weight_in + add_candidates[v]._in
# #         denominator = current_cluster_weight_in + current_cluster_weight_out + add_candidates[
# #             v]._out + self.penalty_value_per_node * (len(current_cluster) + 1)
# #         #####################
# #         # CONSIDER WHAT COULD BE GAINED ON THE NEWLY EXPOSED BOUNDARY
# #         #####################
# #         # find the neighbors of v
# #         distance_2_neighbors = self.graph.hash_graph[v]
# #         for d2n in distance_2_neighbors:
# #             if d2n in current_cluster:
# #                 numerator += factor * self.graph.hash_graph[v][d2n]
# #             else:
# #                 distance_3_neighbors = self.graph.hash_graph[d2n]
# #                 for d3n in distance_3_neighbors:
# #                     if d3n in current_cluster:
# #                         numerator += factor * self.graph.hash_graph[d2n][d3n]
# #                     # elif d3n in add_candidates:
# #                     #     numerator += self.graph.hash_graph[d2n][d3n]
# #                     else:
# #                         denominator += factor * self.graph.hash_graph[d2n][d3n]
# #
# #         proposed_score = numerator / denominator
# #         if proposed_score > best_change_score:
# #             best_change = v
# #             best_change_score = proposed_score
# #         debug("##################### ADD Consideration ########################")
# #         debug("v: %s" % str(v))
# #         debug("proposed_score: %s" % str(proposed_score))
# #         debug("best_change_score: %s" % str(best_change_score))
# #         debug("best_change: %s" % str(best_change))
# #         debug("numerator: %s" % str(numerator))
# #         debug("denominator: %s" % str(denominator))
# #         debug("current_cluster_weight_in: %s" % str(current_cluster_weight_in))
# #         debug("add_candidates[v]._in: %s" % str(add_candidates[v]._in))
# #         debug("current_cluster_weight_out: %s" % str(current_cluster_weight_out))
# #         debug("add_candidates[v]._out: %s" % str(add_candidates[v]._out))
# #         debug("len(current_cluster): %s" % str(len(current_cluster)))
# #         sleep_debug(.25)
# #     return best_change, best_change_score
# #
# #
# # def find_best_add(self, add_candidates, current_cluster, current_score, current_cluster_weight_in, current_cluster_weight_out):
# #     best_change = None
# #     best_change_score = current_score
# #     # best_change_list = []
# #     for v in add_candidates:
# #         numerator = current_cluster_weight_in +\
# #                     add_candidates[v]._in
# #         denominator = current_cluster_weight_in+\
# #                       current_cluster_weight_out +\
# #                       add_candidates[v]._out +\
# #                       self.penalty_value_per_node * (len(current_cluster) + 1)
# #         proposed_score = numerator / denominator
# #         if proposed_score > best_change_score:
# #             best_change = v
# #             best_change_score = proposed_score
# #         # UNFINISHED BEST_CHANGE_LIST
# #         #     best_change_list = [v]
# #         # elif proposed_score == best_change_score:
# #         #     best_change_list.append(v)
# #         debug("##################### ADD Consideration ########################")
# #         debug("v: %s" % str(v))
# #         debug("proposed_score: %s" % str(proposed_score))
# #         debug("best_change_score: %s" % str(best_change_score))
# #         debug("best_change: %s" % str(best_change))
# #         debug("numerator: %s" % str(numerator))
# #         debug("denominator: %s" % str(denominator))
# #         debug("current_cluster_weight_in: %s" % str(current_cluster_weight_in))
# #         debug("add_candidates[v]._in: %s" % str(add_candidates[v]._in))
# #         debug("current_cluster_weight_out: %s" % str(current_cluster_weight_out))
# #         debug("add_candidates[v]._out: %s" % str(add_candidates[v]._out))
# #         debug("len(current_cluster): %s" % str(len(current_cluster)))
# #         sleep_debug(.25)
# #     # UNFINISHED BEST_CHANGE_LIST
# #     # if len( best_change_list)>1:
# #     #     print(best_change_list)
# #     #     for tie in best_change_list:
# #     #         print(add_candidates[tie].stringify())
# #     return best_change, best_change_score
# #
# #
# # def add(self, add_candidates, current_cluster, remove_candidates,
# #         change_vertex, change_vertex_score, cc_weight_in, cc_weight_out):
# #     debug("\n", "ADD: %s" % str(change_vertex), "change_vertex_score: %s" % str(change_vertex_score), "\n")
# #     #################################################################
# #     # update the overall weight into and out of the current_cluster #
# #     #################################################################
# #     change_vertex_in = add_candidates[change_vertex]._in
# #     change_vertex_out = add_candidates[change_vertex]._out
# #     cc_weight_in += change_vertex_in
# #     cc_weight_out = cc_weight_out - change_vertex_in + change_vertex_out
# #
# #     ###################################
# #     # Move the change vertex from add_candidates to current_cluster
# #     ###################################
# #     to_add = add_candidates[change_vertex].copy()
# #     del add_candidates[change_vertex]
# #     current_cluster[change_vertex] = to_add.copy()
# #
# #     ###################################
# #     # Change vertex to remove_candidates if applicable
# #     ###################################
# #     if to_add.num_edges_from:
# #         remove_candidates[change_vertex] = to_add.copy()
# #
# #     def update_v(v, edge_weight, collection):
# #         collection[v]._in += edge_weight
# #         collection[v]._out -= edge_weight
# #         collection[v].num_edges_to += 1
# #         collection[v].num_edges_from -= 1
# #         #######################
# #         # sanity check
# #         #######################
# #         thresh = -.001
# #         a = collection[v]._in < thresh
# #         b = collection[v]._out < thresh
# #         c = collection[v].num_edges_to < thresh
# #         d = collection[v].num_edges_from < thresh
# #         if a or b or c or d:
# #             print("oh no %s; %s%s%s%s"%(str(v), a, b, c, d))
# #             exit()
# #
# #     def initialize_new_add_candidate_for_V(add_vertex, current_cluster, add_candidates):
# #         num_edges_to = 0
# #         weight_to = 0
# #         num_edges_from = 0
# #         weight_from = 0
# #         ###################################
# #         # iterate over the neighbors of v
# #         ###################################
# #         for neighbor in self.graph.hash_graph[add_vertex]:
# #             weight_prime = self.graph.hash_graph[add_vertex][neighbor]
# #             if neighbor in current_cluster:
# #                 num_edges_to += 1
# #                 weight_to += weight_prime
# #             else:
# #                 num_edges_from += 1
# #                 weight_from += weight_prime
# #
# #         add_candidates[add_vertex] = Relationship(weight_to,
# #                                          num_edges_to,
# #                                          weight_from,
# #                                          num_edges_from)
# #
# #     #######################################################################
# #     # iterate over neighbors of change_vertex, and update each Relationship
# #     #######################################################################
# #     for v in self.graph.hash_graph[change_vertex]:
# #         edge_weight = self.graph.hash_graph[v][change_vertex]
# #         if v in add_candidates:
# #             update_v(v, edge_weight, add_candidates)
# #         if v in current_cluster:
# #             update_v(v, edge_weight, current_cluster)
# #         # note that v may be in both the current_cluster and in remove_candidates
# #         # remove_candidates is a subset of current_cluster
# #         if v in remove_candidates:
# #             update_v(v, edge_weight, remove_candidates)
# #             # Check that a candidate for removal is still on the boundary
# #             if remove_candidates[v].num_edges_from == 0:
# #                 del remove_candidates[v]
# #         # handle the case that v is on the new boundary
# #         # add v to add_candidates
# #         if v not in add_candidates and v not in current_cluster:
# #             initialize_new_add_candidate_for_V(v, current_cluster, add_candidates)
# #
# #     return cc_weight_in, cc_weight_out
# #
# #
# # def find_best_remove(self, remove_candidates, current_cluster, current_cluster_weight_in, current_cluster_weight_out, current_score):
# #     best_change = None
# #     best_change_score = current_score
# #     # check that
# #     #   (1) the cluster has more than one element
# #     if len(current_cluster) > 1:
# #         current_cluster_membership_hashset = [vertex for vertex in current_cluster]
# #         for v in remove_candidates:
# #             if self.care_about_cuts:
# #                 # TODO: check if there is a cut.
# #                 #   Implement more efficiently using a Dynamic Connectivity algorithm
# #                 is_a_cut = True
# #                 visted = set()
# #                 start_point = None
# #                 for potential_start_point in current_cluster_membership_hashset:
# #                     if potential_start_point != v:
# #                         start_point = potential_start_point
# #                         # consider break statement here
# #                 # check that
# #                 #   (2) removal of vertex under consideration will not disconnect cluster
# #                 debug("DFS starting vertex: %s" % str(start_point), "DFS ignore vertex: %s" % str(v))
# #                 dfs(self, start_point, v, current_cluster_membership_hashset, visted)
# #                 debug("=============================")
# #                 if len(visted) == -1 + len(current_cluster_membership_hashset):
# #                     is_a_cut = False
# #                     debug("%s is NOT a CUT" % str(v))
# #                 if is_a_cut:
# #                     debug("%s is a CUT!" % str(v))
# #                 debug("cluster: %s" % str(current_cluster_membership_hashset))
# #                 debug("visited by DFS: %s" % str(visted))
# #                 sleep_debug(.25)
# #             else:
# #                 is_a_cut = False
# #
# #             if not is_a_cut:
# #                 # TODO: check that this makes sense
# #                 numerator = current_cluster_weight_in - remove_candidates[v]._in
# #                 denominator = current_cluster_weight_in + current_cluster_weight_out - \
# #                               remove_candidates[v]._out + self.penalty_value_per_node * (
# #                                       len(current_cluster) - 1)
# #                 proposed_score = numerator / denominator
# #                 if proposed_score > best_change_score:
# #                     best_change = v
# #                     best_change_score = proposed_score
# #                 debug("##################### REMOVE Consideration ########################")
# #                 debug("v: %s" % str(v))
# #                 debug("proposed_score: %s" % str(proposed_score))
# #                 debug("best_change_score: %s" % str(best_change_score))
# #                 debug("best_change: %s" % str(best_change))
# #                 debug("numerator: %s" % str(numerator))
# #                 debug("denominator: %s" % str(denominator))
# #                 debug("current_cluster_weight_in: %s" % str(current_cluster_weight_in))
# #                 debug(
# #                     "remove_candidates[v]._in: %s" % str(remove_candidates[v]._in))
# #                 debug("current_cluster_weight_out: %s" % str(current_cluster_weight_out))
# #                 debug("remove_candidates[v]._out: %s" % str(
# #                     remove_candidates[v]._out))
# #                 debug("len(current_cluster): %s" % str(len(current_cluster)))
# #                 sleep_debug(1)
# #     return best_change, best_change_score
# #
# #
# # def remove(self,
# #            remove_candidates,
# #            add_candidates,
# #            current_cluster,
# #            change_vertex,
# #            change_vertex_score,
# #            cc_weight_in,
# #            cc_weight_out):
# #     debug("REMOVE: ", change_vertex, "change_vertex_score: ", change_vertex_score)
# #     sleep_debug(1)
# #     ###############################################################################################
# #     # Update the current_cluster 's score, and overall weight into and out of the current cluster #
# #     ###############################################################################################
# #     change_vertex_in = remove_candidates[change_vertex]._in
# #     change_vertex_out = remove_candidates[change_vertex]._out
# #     cc_weight_in -= change_vertex_in
# #     cc_weight_out = cc_weight_out - change_vertex_out + change_vertex_in
# #
# #     #######################################################################
# #     # Remove the change vertex from remove_candidates and current_cluster #
# #     #######################################################################
# #     to_remove = remove_candidates[change_vertex].copy()
# #     del remove_candidates[change_vertex]
# #     del current_cluster[change_vertex]
# #
# #     ################################################
# #     # Also add the change vertex to add_candidates #
# #     ################################################
# #     add_candidates[change_vertex] = to_remove
# #
# #     #########################################################################
# #     # iterate over neighbors of change_vertex, and update each Relationship #
# #     #########################################################################
# #     def update_v(v, edge_weight, collection):
# #         collection[v]._in -= edge_weight
# #         collection[v]._out += edge_weight
# #         collection[v].num_edges_to -= 1
# #         collection[v].num_edges_from += 1
# #         #######################
# #         # sanity check
# #         #######################
# #         thresh = -.001
# #         a = collection[v]._in < thresh
# #         b = collection[v]._out < thresh
# #         c = collection[v].num_edges_to < thresh
# #         d = collection[v].num_edges_from < thresh
# #         if a or b or c or d:
# #             print("oh no %s; %s%s%s%s"%(str(v), a, b, c, d))
# #             exit()
# #
# #     for v in self.graph.hash_graph[change_vertex]:
# #         edge_weight = self.graph.hash_graph[change_vertex][v]
# #         # note that v may be in both the current_cluster and in remove_candidates
# #         if v in remove_candidates:
# #             update_v(v, edge_weight, remove_candidates)
# #         if v in current_cluster:
# #             update_v(v, edge_weight, current_cluster)
# #             if current_cluster[v].num_edges_from == 1:
# #                 remove_candidates[v] = current_cluster[v].copy()
# #         if v in add_candidates:
# #             update_v(v, edge_weight, add_candidates)
# #             if add_candidates[v].num_edges_to == 0:
# #                 del add_candidates[v]
# #     return cc_weight_in, cc_weight_out
# #
# #
# # def find_best_suboptimal_add(self, add_candidates, current_cluster, cc_weight_in, cc_weight_out):
# #     best_suboptimal_change = None
# #     best_suboptimal_change_score = -10000
# #     for v in add_candidates:
# #         numerator = cc_weight_in + add_candidates[v]._in
# #         denominator = cc_weight_in + cc_weight_out + add_candidates[
# #             v]._out + self.penalty_value_per_node * (len(current_cluster) + 1)
# #         proposed_score = numerator / denominator
# #         if proposed_score > best_suboptimal_change_score:
# #             best_suboptimal_change = v
# #             best_suboptimal_change_score = proposed_score
# #         debug("##################### SUBOPTIMAL ADD Consideration ########################")
# #         debug("v: %s" % str(v))
# #         debug("proposed_score: %s" % str(proposed_score))
# #         debug("best_suboptimal_change_score: %s" % str(best_suboptimal_change_score))
# #         debug("best_suboptimal_change: %s" % str(best_suboptimal_change))
# #         debug("numerator: %s" % str(numerator))
# #         debug("denominator: %s" % str(denominator))
# #         debug("current_cluster_weight_in: %s" % str(cc_weight_in))
# #         debug("add_candidates[v]._in: %s" % str(add_candidates[v]._in))
# #         debug("current_cluster_weight_out: %s" % str(cc_weight_out))
# #         debug("add_candidates[v]._out: %s" % str(add_candidates[v]._out))
# #         debug("len(current_cluster): %s" % str(len(current_cluster)))
# #         sleep_debug(.25)
# #     return best_suboptimal_change, best_suboptimal_change_score
# #
# #
# # def add_shake(self, add_candidates, remove_candidates, current_cluster, cc_weight_in, cc_weight_out, round_no, last_failed_add_round_no):
# #     for i in range(self.number_of_bad_adds):
# #         best_suboptimal_change, best_suboptimal_score = find_best_suboptimal_add(self, add_candidates, current_cluster,
# #                                                                                  cc_weight_in, cc_weight_out)
# #         if best_suboptimal_change:
# #             print("adding_suboptimally")
# #             sleep_debug(1)
# #             round_no += 1
# #             last_failed_add_round_no = -5
# #             # TODO handle round_no numbers
# #             cc_weight_in, cc_weight_out = add(self, add_candidates, current_cluster, remove_candidates,
# #                                               best_suboptimal_change, best_suboptimal_score, cc_weight_in,
# #                                               cc_weight_out)
# #     return best_suboptimal_score, \
# #            cc_weight_in, \
# #            cc_weight_out, \
# #            round_no, \
# #            last_failed_add_round_no