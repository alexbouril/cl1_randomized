from common import *
from graph import *

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


# maybe if we add vertex that is best in terms of 2 neighborhood
#       but may or may not be the best in terms of 1 neighborhood
# after we add this vertex, we add the helpful neighbor vertices, without considering the whole boundary
def find_best_2neighborhood_add(self, add_candidates, current_cluster, current_score, current_cluster_weight_in, current_cluster_weight_out):
    best_change = None
    best_change_score = current_score
    for v in add_candidates:
        numerator = current_cluster_weight_in + add_candidates[v].sum_weight_to
        denominator = current_cluster_weight_in + current_cluster_weight_out + add_candidates[
            v].sum_weight_from + self.penalty_value_per_node * (len(current_cluster) + 1)
        #####################
        # CONSIDER WHAT COULD BE GAINED ON THE NEWLY EXPOSED BOUNDARY
        #####################
        # find the neighbors of v
        neighbors = self.graph.hash_graph[v]
        for neighbor in neighbors:
            # find out who the neighbor is connected to
            distance_two_neighbors = self.graph.hash_graph[neighbor]
            for d2n in distance_two_neighbors:
                distance_three_neighbors = self.graph.hash_graph[d2n]
                for d3n in distance_three_neighbors:
                    # weight to {current cluster, weight to v, weight to neighbors}

                    numerator+=0

                    #weight from {current cluster, weight to v, weight to neighbors}
                    denominator+=0

        proposed_score = numerator / denominator
        if proposed_score > best_change_score:
            best_change = v
            best_change_score = proposed_score
        debug("##################### ADD Consideration ########################")
        debug("v: %s" % str(v))
        debug("proposed_score: %s" % str(proposed_score))
        debug("best_change_score: %s" % str(best_change_score))
        debug("best_change: %s" % str(best_change))
        debug("numerator: %s" % str(numerator))
        debug("denominator: %s" % str(denominator))
        debug("current_cluster_weight_in: %s" % str(current_cluster_weight_in))
        debug("add_candidates[v].sum_weight_to: %s" % str(add_candidates[v].sum_weight_to))
        debug("current_cluster_weight_out: %s" % str(current_cluster_weight_out))
        debug("add_candidates[v].sum_weight_from: %s" % str(add_candidates[v].sum_weight_from))
        debug("len(current_cluster): %s" % str(len(current_cluster)))
        sleep_debug(.25)
    return best_change, best_change_score


def find_best_add(self, add_candidates, current_cluster, current_score, current_cluster_weight_in, current_cluster_weight_out):
    best_change = None
    best_change_score = current_score
    for v in add_candidates:
        numerator = current_cluster_weight_in + add_candidates[v].sum_weight_to
        denominator = current_cluster_weight_in + \
                      current_cluster_weight_out + \
                      add_candidates[v].sum_weight_from + \
                      self.penalty_value_per_node * (len(current_cluster) + 1)
        proposed_score = numerator / denominator
        if proposed_score == best_change_score:
            print("oh wait, there was a removal tie")
            time.sleep(5)
        if proposed_score > best_change_score:
            best_change = v
            best_change_score = proposed_score
        debug("##################### ADD Consideration ########################")
        debug("v: %s" % str(v))
        debug("proposed_score: %s" % str(proposed_score))
        debug("best_change_score: %s" % str(best_change_score))
        debug("best_change: %s" % str(best_change))
        debug("numerator: %s" % str(numerator))
        debug("denominator: %s" % str(denominator))
        debug("current_cluster_weight_in: %s" % str(current_cluster_weight_in))
        debug("add_candidates[v].sum_weight_to: %s" % str(add_candidates[v].sum_weight_to))
        debug("current_cluster_weight_out: %s" % str(current_cluster_weight_out))
        debug("add_candidates[v].sum_weight_from: %s" % str(add_candidates[v].sum_weight_from))
        debug("len(current_cluster): %s" % str(len(current_cluster)))
        sleep_debug(.25)
    return best_change, best_change_score


def add(self, add_candidates, current_cluster, remove_candidates,
        change_vertex, change_vertex_score, cc_weight_in, cc_weight_out):
    debug("\n", "ADD: %s" % str(change_vertex), "change_vertex_score: %s" % str(change_vertex_score), "\n")
    # update the overall weight into and out of the current_cluster
    cc_weight_in += add_candidates[change_vertex].sum_weight_to
    cc_weight_out += add_candidates[change_vertex].sum_weight_from

    ###################################
    # Move the change vertex from add_candidates to current_cluster
    ###################################
    to_add = add_candidates[change_vertex].copy()
    del add_candidates[change_vertex]
    current_cluster[change_vertex] = to_add.copy()

    ###################################
    # Change vertex to remove_candidates if applicable
    ###################################
    if to_add.num_edges_from:
        remove_candidates[change_vertex] = to_add.copy()

    s = current_cluster[change_vertex].stringify()
    def update_v(v, edge_weight, collection):
        # a = collection[v].sum_weight_to < 0
        # b = collection[v].sum_weight_from < 0
        # c = collection[v].num_edges_to < 0
        # d = collection[v].num_edges_from < 0
        #
        # if a or b or c or d:
        #     print("darn %s; %s%s%s%s"%(str(v), a, b, c, d))
        #     # exit()

        collection[v].sum_weight_to += edge_weight
        collection[v].sum_weight_from -= edge_weight
        collection[v].num_edges_to += 1
        collection[v].num_edges_from -= 1

        thresh = -.001
        a = collection[v].sum_weight_to < thresh
        b = collection[v].sum_weight_from < thresh
        c = collection[v].num_edges_to < thresh
        d = collection[v].num_edges_from < thresh
        if a or b or c or d:
            print("oh no %s; %s%s%s%s"%(str(v), a, b, c, d))
            exit()


    #######################################################################
    # iterate over neighbors of change_vertex, and update each Relationship
    #######################################################################
    for v in self.graph.hash_graph[change_vertex]:
        edge_weight = self.graph.hash_graph[v][change_vertex]
        if v in add_candidates:
            update_v(v, edge_weight, add_candidates)
        if v in current_cluster:
            update_v(v, edge_weight, current_cluster)
        # note that v may be in both the current_cluster and in remove_candidates
        # remove_candidates is a subset of current_cluster
        if v in remove_candidates:
            update_v(v, edge_weight, remove_candidates)
            # Check that a candidate for removal+ is still on the boundary
            if remove_candidates[v].num_edges_from == 0:
                del remove_candidates[v]
        # handle the case that v is on the new boundary
        # add v to add_candidates
        if v not in add_candidates and v not in current_cluster:
            num_edges_to = 0
            weight_to = 0
            num_edges_from = 0
            weight_from = 0
            ###################################
            # iterate over the neighbors of v
            ###################################
            for v_prime in self.graph.hash_graph[v]:
                weight_prime = self.graph.hash_graph[v][v_prime]
                if v_prime in current_cluster:
                    num_edges_to += 1
                    weight_to += weight_prime
                else:
                    num_edges_from += 1
                    weight_from += weight_prime

            add_candidates[v] = Relationship(weight_to,
                                             num_edges_to,
                                             weight_from,
                                             num_edges_from)
            # a = num_edges_to < 0
            # b = weight_to < 0
            # c = num_edges_from < 0
            # d = weight_from < 0
            # if a or b or c or d:
            #     print("shoot")
            #     # exit(0)


    return cc_weight_in, cc_weight_out


def find_best_remove(self, remove_candidates, current_cluster, current_cluster_weight_in, current_cluster_weight_out, current_score):
    best_change = None
    best_change_score = current_score
    # check that
    #   (1) the cluster has more than one element
    if len(current_cluster) > 1:
        current_cluster_membership_hashset = [vertex for vertex in current_cluster]
        for v in remove_candidates:
            if self.care_about_cuts:
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
                debug("DFS starting vertex: %s" % str(start_point), "DFS ignore vertex: %s" % str(v))
                dfs(self, start_point, v, current_cluster_membership_hashset, visted)
                debug("=============================")
                if len(visted) == -1 + len(current_cluster_membership_hashset):
                    is_a_cut = False
                    debug("%s is NOT a CUT" % str(v))
                if is_a_cut:
                    debug("%s is a CUT!" % str(v))
                debug("cluster: %s" % str(current_cluster_membership_hashset))
                debug("visited by DFS: %s" % str(visted))
                sleep_debug(.25)
            else:
                is_a_cut = False

            if not is_a_cut:
                # TODO: check that this makes sense
                numerator = current_cluster_weight_in - \
                            remove_candidates[v].sum_weight_to
                denominator = current_cluster_weight_in + \
                              current_cluster_weight_out - \
                              remove_candidates[v].sum_weight_from + \
                              self.penalty_value_per_node * (len(current_cluster) - 1)
                proposed_score = numerator / denominator
                if proposed_score==best_change_score:
                    print("oh wait, there was a removal tie")
                    time.sleep(5)
                if proposed_score > best_change_score:
                    best_change = v
                    best_change_score = proposed_score
                debug("##################### REMOVE Consideration ########################")
                debug("v: %s" % str(v))
                debug("proposed_score: %s" % str(proposed_score))
                debug("best_change_score: %s" % str(best_change_score))
                debug("best_change: %s" % str(best_change))
                debug("numerator: %s" % str(numerator))
                debug("denominator: %s" % str(denominator))
                debug("current_cluster_weight_in: %s" % str(current_cluster_weight_in))
                debug(
                    "remove_candidates[v].sum_weight_to: %s" % str(remove_candidates[v].sum_weight_to))
                debug("current_cluster_weight_out: %s" % str(current_cluster_weight_out))
                debug("remove_candidates[v].sum_weight_from: %s" % str(
                    remove_candidates[v].sum_weight_from))
                debug("len(current_cluster): %s" % str(len(current_cluster)))
                sleep_debug(1)
    return best_change, best_change_score


def remove(self, remove_candidates, add_candidates, current_cluster, change_vertex, change_vertex_score, cc_weight_in, cc_weight_out):
    debug("REMOVE: ", change_vertex, "change_vertex_score: ", change_vertex_score)
    sleep_debug(1)
    # Update the current_cluster 's score, and overall weight into and out of the current cluster
    current_score = change_vertex_score
    cc_weight_in -= remove_candidates[change_vertex].sum_weight_to
    cc_weight_out -= remove_candidates[change_vertex].sum_weight_from

    # Remove the change vertex from remove_candidates and current_cluster
    to_remove = remove_candidates[change_vertex].copy()
    del remove_candidates[change_vertex]
    del current_cluster[change_vertex]

    # Also add the change vertex to add_candidates
    add_candidates[change_vertex] = to_remove.copy()

    def update_v(v, edge_weight, collection):
        collection[v].sum_weight_to -= edge_weight
        collection[v].sum_weight_from += edge_weight
        collection[v].num_edges_to -= 1
        collection[v].num_edges_from += 1

    # AFTER this is done, THEN do the following
    for v in self.graph.hash_graph[
        change_vertex]:  # iterate over neighbors of change_vertex, and update each Relationship
        edge_weight = self.graph.hash_graph[change_vertex][v]
        # note that v may be in both the current_cluster and in remove_candidates
        if v in remove_candidates:
            update_v(v, edge_weight, remove_candidates)
        if v in current_cluster:
            update_v(v, edge_weight, current_cluster)
            # TODO: if v is in the current cluster, and is newly also on the boundary, add v to remove_candidates
            if current_cluster[v].num_edges_from == 1:
                remove_candidates[v] = current_cluster[v].copy()
        if v in add_candidates:
            update_v(v, edge_weight, add_candidates)
            # TODO: if v in no longer on the boundary after change_vertex is removed from current_cluster, remove v from add_candidates
            if add_candidates[v].num_edges_to == 0:
                del add_candidates[v]
    return cc_weight_in, cc_weight_out


def find_best_suboptimal_add(self, add_candidates, current_cluster, cc_weight_in, cc_weight_out):
    best_suboptimal_change = None
    best_suboptimal_change_score = -10000
    for v in add_candidates:
        numerator = cc_weight_in + add_candidates[v].sum_weight_to
        denominator = cc_weight_in + cc_weight_out + add_candidates[
            v].sum_weight_from + self.penalty_value_per_node * (len(current_cluster) + 1)
        proposed_score = numerator / denominator
        if proposed_score > best_suboptimal_change_score:
            best_suboptimal_change = v
            best_suboptimal_change_score = proposed_score
        debug("##################### SUBOPTIMAL ADD Consideration ########################")
        debug("v: %s" % str(v))
        debug("proposed_score: %s" % str(proposed_score))
        debug("best_suboptimal_change_score: %s" % str(best_suboptimal_change_score))
        debug("best_suboptimal_change: %s" % str(best_suboptimal_change))
        debug("numerator: %s" % str(numerator))
        debug("denominator: %s" % str(denominator))
        debug("current_cluster_weight_in: %s" % str(cc_weight_in))
        debug("add_candidates[v].sum_weight_to: %s" % str(add_candidates[v].sum_weight_to))
        debug("current_cluster_weight_out: %s" % str(cc_weight_out))
        debug("add_candidates[v].sum_weight_from: %s" % str(add_candidates[v].sum_weight_from))
        debug("len(current_cluster): %s" % str(len(current_cluster)))
        sleep_debug(.25)
    return best_suboptimal_change, best_suboptimal_change_score


def add_shake(self, add_candidates, remove_candidates, current_cluster, cc_weight_in, cc_weight_out, round_no, last_failed_add_round_no):
    for i in range(self.number_of_bad_adds):
        best_suboptimal_change, best_suboptimal_score = find_best_suboptimal_add(self, add_candidates, current_cluster,
                                                                                 cc_weight_in, cc_weight_out)
        if best_suboptimal_change:
            print("adding_suboptimally")
            sleep_debug(1)
            round_no += 1
            last_failed_add_round_no = -5
            # TODO handle round_no numbers
            cc_weight_in, cc_weight_out = add(self, add_candidates, current_cluster, remove_candidates,
                                              best_suboptimal_change, best_suboptimal_score, cc_weight_in,
                                              cc_weight_out)
    return best_suboptimal_score, \
           cc_weight_in, \
           cc_weight_out, \
           round_no, \
           last_failed_add_round_no