from common import *

def original_construction(self):
    # TODO: move dfs outside of original_construction()
    def dfs(current_vertex, ignore_vertex, current_cluster_membership_hashset, visited):
        visited.add(current_vertex)
        for neighbor in self.graph.hash_graph[current_vertex]:
            if neighbor not in current_cluster_membership_hashset:
                continue
            elif neighbor in visited:
                continue
            elif neighbor == ignore_vertex:
                continue
            else:
                dfs(neighbor, ignore_vertex, current_cluster_membership_hashset, visited)

    considered_vertices = set()
    index = 0

    while index < len(self.vertices_by_degree):
        current_seed = self.vertices_by_degree[index][0]
        current_seed_degree = self.vertices_by_degree[index][1]
        debug("current_seed: ", current_seed)
        debug("current_seed_degree: ", current_seed_degree)
        if current_seed in considered_vertices:
            debug("SKIP %s" % str(current_seed))
            index += 1
            continue
        else:
            debug("Starting cluster #%s" % str(len(self.initial_clustering)))
            # initialize the current cluster starting with the selected seed
            current_cluster = dict()
            weight_to = 0
            num_edges_to = 0
            weight_from = sum([self.graph.hash_graph[current_seed][tar] for tar in self.graph.hash_graph[current_seed]])
            num_edges_from = len(self.graph.hash_graph[current_seed])
            current_cluster[current_seed] = Relationship(weight_to, num_edges_to, weight_from, num_edges_from)
            current_score = 0
            current_cluster_weight_in = 0
            current_cluster_weight_out = weight_from
            # initalize the candidates for removal
            remove_candidates = dict()
            remove_candidates[current_seed] = Relationship(weight_to, num_edges_to, weight_from, num_edges_from)
            # initialize the candidates for addition
            add_candidates = dict()
            for target in self.graph.hash_graph[current_seed]:
                weight_to = self.graph.hash_graph[current_seed][target]
                num_edges = len(self.graph.hash_graph[target])
                num_edges_to = 1
                num_edges_from = num_edges - 1
                weight_from = sum([self.graph.hash_graph[target][tar] for tar in self.graph.hash_graph[target] if
                                   tar != current_seed])
                add_candidates[target] = Relationship(weight_to, num_edges_to, weight_from, num_edges_from)

            def find_best_add():
                best_change = None
                best_change_score = current_score
                for v in add_candidates:
                    numerator = current_cluster_weight_in + add_candidates[v].sum_weight_to
                    denominator = current_cluster_weight_in + current_cluster_weight_out + add_candidates[
                        v].sum_weight_from + self.penalty_value_per_node * (len(current_cluster) + 1)
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

            def add(change_vertex, change_vertex_score, cc_weight_in, cc_weight_out):
                debug("\n", "ADD: %s" % str(change_vertex), "change_vertex_score: %s" % str(change_vertex_score), "\n")
                # update the overall weight into and out of the current_cluster
                cc_weight_in += add_candidates[change_vertex].sum_weight_to
                cc_weight_out += add_candidates[change_vertex].sum_weight_from

                # Move the change vertex from add_candidates to current_cluster
                to_add = add_candidates[change_vertex].copy()
                current_cluster[change_vertex] = to_add
                del add_candidates[change_vertex]

                # Also add the change vertex to remove_candidates if applicable
                if to_add.num_edges_from:
                    remove_candidates[change_vertex] = to_add

                for v in self.graph.hash_graph[
                    change_vertex]:  # iterate over neighbors of change_vertex, and update each Relationship
                    edge_weight = self.graph.hash_graph[change_vertex][v]
                    if v in add_candidates:
                        add_candidates[v].sum_weight_to = edge_weight + add_candidates[v].sum_weight_to
                        add_candidates[v].sum_weight_from = -1 * edge_weight + add_candidates[v].sum_weight_from
                        add_candidates[v].num_edges_to = 1 + add_candidates[v].num_edges_to
                        add_candidates[v].num_edges_from = -1 + add_candidates[v].num_edges_from
                    if v in current_cluster:
                        current_cluster[v].sum_weight_to = edge_weight + current_cluster[v].sum_weight_to
                        current_cluster[v].sum_weight_from = -1 * edge_weight + current_cluster[v].sum_weight_from
                        current_cluster[v].num_edges_to = 1 + current_cluster[v].num_edges_to
                        current_cluster[v].num_edges_from = -1 + current_cluster[v].num_edges_from
                    # note that v may be in both the current_cluster and in remove_candidates
                    # remove_candidates is a subset of current_cluster
                    if v in remove_candidates:
                        remove_candidates[v].sum_weight_to = edge_weight + remove_candidates[v].sum_weight_to
                        remove_candidates[v].sum_weight_from = -1 * edge_weight + remove_candidates[
                            v].sum_weight_from
                        remove_candidates[v].num_edges_to = 1 + remove_candidates[v].num_edges_to
                        remove_candidates[v].num_edges_from = -1 + remove_candidates[v].num_edges_from
                        # Check that a candidate for removal is still on the boundary
                        if remove_candidates[v].num_edges_from == 0:
                            del remove_candidates[v]
                    # handle the case that v is on the new boundary
                    # add v to add_candidates
                    if v not in add_candidates and v not in current_cluster:
                        num_edges_to = 0
                        weight_to = 0
                        num_edges_from = 0
                        weight_from = 0
                        for v_prime in self.graph.hash_graph[v]:
                            weight_prime = self.graph.hash_graph[v][v_prime]
                            if v_prime in current_cluster:
                                num_edges_to += 1
                                weight_to += weight_prime
                            else:
                                num_edges_from += 1
                                weight_from += weight_prime
                        add_candidates[v] = Relationship(weight_to, num_edges_to, weight_from, num_edges_from)
                return cc_weight_in, cc_weight_out

            def find_best_remove(current_score):
                best_change = None
                best_change_score = current_score
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
                            dfs(start_point, v, current_cluster_membership_hashset, visted)
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
                            numerator = current_cluster_weight_in - remove_candidates[v].sum_weight_to
                            denominator = current_cluster_weight_in + current_cluster_weight_out - \
                                          remove_candidates[v].sum_weight_from + self.penalty_value_per_node * (
                                                  len(current_cluster) - 1)
                            proposed_score = numerator / denominator
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

            def remove(change_vertex, change_vertex_score, cc_weight_in, cc_weight_out):
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
                add_candidates[change_vertex] = to_remove

                # AFTER this is done, THEN do the following
                for v in self.graph.hash_graph[
                    change_vertex]:  # iterate over neighbors of change_vertex, and update each Relationship
                    edge_weight = self.graph.hash_graph[change_vertex][v]
                    # note that v may be in both the current_cluster and in remove_candidates
                    if v in remove_candidates:
                        remove_candidates[v].sum_weight_to = -1 * edge_weight + remove_candidates[v].sum_weight_to
                        remove_candidates[v].sum_weight_from = edge_weight + remove_candidates[v].sum_weight_from
                        remove_candidates[v].num_edges_to = -1 + remove_candidates[v].num_edges_to
                        remove_candidates[v].num_edges_from = 1 + remove_candidates[v].num_edges_from
                    if v in current_cluster:
                        current_cluster[v].sum_weight_to = -1 * edge_weight + current_cluster[v].sum_weight_to
                        current_cluster[v].sum_weight_from = edge_weight + current_cluster[v].sum_weight_from
                        current_cluster[v].num_edges_to = -1 + current_cluster[v].num_edges_to
                        current_cluster[v].num_edges_from = 1 + current_cluster[v].num_edges_from
                        # TODO: if v is in the current cluster, and is newly also on the boundary, add v to remove_candidates
                        if current_cluster[v].num_edges_from == 1:
                            remove_candidates[v] = current_cluster[v].copy()
                    if v in add_candidates:
                        # update the relationship
                        add_candidates[v].sum_weight_to = -1 * edge_weight + add_candidates[v].sum_weight_to
                        add_candidates[v].sum_weight_from = edge_weight + add_candidates[v].sum_weight_from
                        add_candidates[v].num_edges_to = - 1 + add_candidates[v].num_edges_to
                        add_candidates[v].num_edges_from = 1 + add_candidates[v].num_edges_from
                        # TODO: if v in no longer on the boundary after change_vertex is removed from current_cluster, remove v from add_candidates
                        if add_candidates[v].num_edges_to == 0:
                            del add_candidates[v]
                return cc_weight_in, cc_weight_out

            improvement_flag = True
            while improvement_flag and (add_candidates or remove_candidates):
                improvement_flag = False
                ############################################################
                # CONSIDER ADDING A VERTEX ON THE BOUNDARY
                ############################################################
                best_add_change, best_add_change_score = find_best_add()

                ############################################################
                # CONSIDER REMOVING A VERTEX ON THE BOUNDARY
                ############################################################
                best_remove_change, best_remove_change_score = find_best_remove(current_score)

                ############################################################
                # IF THERE ARE VERTICES THAT CAN BE ADDED OR REMOVED TO IMPROVE COHESIVENESS, TAKE THE ACTION THAT RESULTS IN THE GREATEST IMPROVEMENT
                ############################################################
                should_add = False
                should_remove = False
                if best_add_change and best_remove_change:
                    if best_add_change_score > best_remove_change_score:
                        should_add = True
                    else:
                        should_remove = True
                elif best_add_change:
                    should_add = True
                elif best_remove_change:
                    should_remove = True

                if should_add:
                    improvement_flag = True
                    current_score = best_add_change_score
                    add(best_add_change, best_add_change_score, current_cluster_weight_in, current_cluster_weight_out)
                elif should_remove:
                    improvement_flag = True
                    current_score = best_remove_change_score
                    remove(best_remove_change, best_remove_change_score, current_cluster_weight_in,
                           current_cluster_weight_out)

            # add current_cluster to the list of clusters
            self.initial_clustering.append(current_cluster)
            index += 1
            if not self.seed_from_all:
                for v in current_cluster:
                    considered_vertices.add(v)
            print("CLUSTER #%s: %s" % (str(len(self.initial_clustering)), str([vertex for vertex in current_cluster])))
            # time.sleep(2)


def original_construction_2(self):
    # TODO: move dfs outside of original_construction()
    def dfs(current_vertex, ignore_vertex, current_cluster_membership_hashset, visited):
        visited.add(current_vertex)
        for neighbor in self.graph.hash_graph[current_vertex]:
            if neighbor not in current_cluster_membership_hashset:
                continue
            elif neighbor in visited:
                continue
            elif neighbor == ignore_vertex:
                continue
            else:
                dfs(neighbor, ignore_vertex, current_cluster_membership_hashset, visited)

    considered_vertices = set()
    index = 0

    while index < len(self.vertices_by_degree):
        current_seed = self.vertices_by_degree[index][0]
        current_seed_degree = self.vertices_by_degree[index][1]
        debug("current_seed: ", current_seed)
        debug("current_seed_degree: ", current_seed_degree)
        if current_seed in considered_vertices:
            debug("SKIP %s" % str(current_seed))
            index += 1
            continue
        else:
            debug("Starting cluster #%s" % str(len(self.initial_clustering)))
            # initialize the current cluster starting with the selected seed
            current_cluster = dict()
            weight_to = 0
            num_edges_to = 0
            weight_from = sum([self.graph.hash_graph[current_seed][tar] for tar in self.graph.hash_graph[current_seed]])
            num_edges_from = len(self.graph.hash_graph[current_seed])
            current_cluster[current_seed] = Relationship(weight_to, num_edges_to, weight_from, num_edges_from)
            current_score = 0
            current_cluster_weight_in = 0
            current_cluster_weight_out = weight_from
            # initalize the candidates for removal
            remove_candidates = dict()
            remove_candidates[current_seed] = Relationship(weight_to, num_edges_to, weight_from, num_edges_from)
            # initialize the candidates for addition
            add_candidates = dict()
            for target in self.graph.hash_graph[current_seed]:
                weight_to = self.graph.hash_graph[current_seed][target]
                num_edges = len(self.graph.hash_graph[target])
                num_edges_to = 1
                num_edges_from = num_edges - 1
                weight_from = sum([self.graph.hash_graph[target][tar] for tar in self.graph.hash_graph[target] if
                                   tar != current_seed])
                add_candidates[target] = Relationship(weight_to, num_edges_to, weight_from, num_edges_from)

            improvement_flag = True

            while improvement_flag and (add_candidates or remove_candidates):
                improvement_flag = False
                ############################################################
                # CONSIDER ADDING A VERTEX ON THE BOUNDARY
                ############################################################
                best_change = None
                best_change_score = current_score
                for v in add_candidates:
                    numerator = current_cluster_weight_in + add_candidates[v].sum_weight_to
                    # TODO: check denominator
                    denominator = current_cluster_weight_in + current_cluster_weight_out + add_candidates[
                        v].sum_weight_from + self.penalty_value_per_node * (len(current_cluster) + 1)

                    proposed_score = numerator / denominator

                    if proposed_score > best_change_score:
                        best_change = v
                        best_change_score = proposed_score
                    debug("##################### ADD Consideration ########################")
                    debug("v: %s" % str(v), "proposed_score: %s" % str(proposed_score))
                    debug("best_change_score: %s" % str(best_change_score), "best_change: %s" % str(best_change))
                    debug("numerator: %s" % str(numerator), "denominator: %s" % str(denominator))
                    debug("current_cluster_weight_in: %s" % str(current_cluster_weight_in))
                    debug("add_candidates[v].sum_weight_to: %s" % str(add_candidates[v].sum_weight_to))
                    debug("current_cluster_weight_out: %s" % str(current_cluster_weight_out))
                    debug("add_candidates[v].sum_weight_from: %s" % str(add_candidates[v].sum_weight_from))
                    sleep_debug(.25)

                ############################################################
                # IF THERE ARE VERTICES THAT CAN BE ADDED TO IMPROVE COHESIVENESS, ADD THE VERTEX THAT RESULTS IN THE GREATEST IMPROVEMENT
                ############################################################
                if best_change:
                    debug("\n", "ADD: %s" % str(best_change), "best_change_score: %s" % str(best_change_score), "\n")
                    improvement_flag = True

                    # update the current_cluster 's score
                    current_score = best_change_score
                    # update the overall weight into and out of the current_cluster
                    current_cluster_weight_in += add_candidates[best_change].sum_weight_to
                    current_cluster_weight_out += add_candidates[best_change].sum_weight_from
                    # copy the relationship, and add the change vertex to current_cluster
                    to_add = add_candidates[best_change].copy()
                    current_cluster[best_change] = to_add
                    # remove the change vertex from add candidates
                    del add_candidates[best_change]
                    # Also add the change vertex to remove_candidates if the change vertex is on the boundary
                    if to_add.num_edges_from:
                        remove_candidates[best_change] = to_add
                    ############################################################
                    # iterate over neighbors of best_change, and update each Relationship
                    ############################################################
                    for v in self.graph.hash_graph[best_change]:
                        edge_weight = self.graph.hash_graph[best_change][v]
                        if v in add_candidates:
                            add_candidates[v].sum_weight_to = edge_weight + add_candidates[v].sum_weight_to
                            add_candidates[v].sum_weight_from = -1 * edge_weight + add_candidates[v].sum_weight_from
                            add_candidates[v].num_edges_to = 1 + add_candidates[v].num_edges_to
                            add_candidates[v].num_edges_from = -1 + add_candidates[v].num_edges_from
                        if v in current_cluster:
                            current_cluster[v].sum_weight_to = edge_weight + current_cluster[v].sum_weight_to
                            current_cluster[v].sum_weight_from = -1 * edge_weight + current_cluster[v].sum_weight_from
                            current_cluster[v].num_edges_to = 1 + current_cluster[v].num_edges_to
                            current_cluster[v].num_edges_from = -1 + current_cluster[v].num_edges_from
                        # NOTE: v may be in both the current_cluster and in remove_candidates
                        # remove_candidates is a subset of current_cluster
                        if v in remove_candidates:
                            remove_candidates[v].sum_weight_to = edge_weight + remove_candidates[v].sum_weight_to
                            remove_candidates[v].sum_weight_from = -1 * edge_weight + remove_candidates[
                                v].sum_weight_from
                            remove_candidates[v].num_edges_to = 1 + remove_candidates[v].num_edges_to
                            remove_candidates[v].num_edges_from = -1 + remove_candidates[v].num_edges_from
                            # Check that a candidate for removal is still on the boundary
                            if remove_candidates[v].num_edges_from == 0:
                                del remove_candidates[v]
                        # handle the case that v is on the new boundary
                        if v not in add_candidates and v not in current_cluster:
                            # add v to add_candidates
                            num_edges_to = 0
                            weight_to = 0
                            num_edges_from = 0
                            weight_from = 0
                            for v_prime in self.graph.hash_graph[v]:
                                weight_prime = self.graph.hash_graph[v][v_prime]
                                if v_prime in current_cluster:
                                    num_edges_to += 1
                                    weight_to += weight_prime
                                else:
                                    num_edges_from += 1
                                    weight_from += weight_prime
                            add_candidates[v] = Relationship(weight_to, num_edges_to, weight_from, num_edges_from)

                else:
                    debug("\n", "No improvement by ADDING", "\n")

                ############################################################
                # CONSIDER REMOVING A VERTEX ON THE BOUNDARY
                ############################################################
                # check that the cluster has more than one element
                if len(current_cluster) > 1:
                    best_change = None
                    best_change_score = current_score
                    current_cluster_membership_hashset = [vertex for vertex in current_cluster]
                    for v in remove_candidates:
                        # TODO: check if there is a cut.
                        #   Implement more efficiently using a Dynamic Connectivity algorithm
                        is_a_cut = True
                        visted = set()
                        start_point = None
                        # find a start point for dfs.  the start point should not be 'v'.
                        for potential_start_point in current_cluster_membership_hashset:
                            if potential_start_point != v:
                                start_point = potential_start_point
                                # TODO: consider break statement here
                        # check that removal of vertex under consideration will not disconnect cluster
                        debug("DFS starting vertex: %s" % str(start_point), "DFS ignore vertex: %s" % str(v),
                              "=============================")
                        dfs(start_point, v, current_cluster_membership_hashset, visted)
                        if len(visted) == -1 + len(current_cluster_membership_hashset):
                            is_a_cut = False
                        debug("vertex %s is a cut? %s" % (str(v), str(is_a_cut)),
                              "cluster: %s" % str(current_cluster_membership_hashset),
                              "visited by DFS: %s" % str(visted))
                        sleep_debug(.25)
                        if not is_a_cut:
                            # TODO: check that this makes sense
                            numerator = current_cluster_weight_in - remove_candidates[v].sum_weight_to
                            denominator = current_cluster_weight_in + current_cluster_weight_out - remove_candidates[
                                v].sum_weight_from + self.penalty_value_per_node * (len(current_cluster) - 1)
                            proposed_score = numerator / denominator
                            if proposed_score > best_change_score:
                                best_change = v
                                best_change_score = proposed_score
                            debug("##################### REMOVE Consideration ########################")
                            debug("v: %s" % str(v), "proposed_score: %s" % str(proposed_score))
                            debug("best_change_score: %s" % str(best_change_score),
                                  "best_change: %s" % str(best_change))
                            debug("numerator: %s" % str(numerator), "denominator: %s" % str(denominator))
                            debug("current_cluster_weight_in: %s" % str(current_cluster_weight_in))
                            debug("remove_candidates[v].sum_weight_to: %s" % str(remove_candidates[v].sum_weight_to))
                            debug("current_cluster_weight_out: %s" % str(current_cluster_weight_out))
                            debug(
                                "remove_candidates[v].sum_weight_from: %s" % str(remove_candidates[v].sum_weight_from))
                            sleep_debug(.5)
                    ############################################################
                    # IF THERE ARE VERTICES THAT CAN BE REMOVED TO IMPROVE COHESIVENESS, REMOVE THE VERTEX THAT RESULTS IN THE GREATEST IMPROVEMENT
                    ############################################################
                    if best_change:
                        improvement_flag = True
                        debug("REMOVE: ", best_change, "best_change_score: ", best_change_score)
                        sleep_debug(2)
                        # Update the current_cluster 's score, and overall weight into and out of the current cluster
                        current_score = best_change_score
                        current_cluster_weight_in -= remove_candidates[best_change].sum_weight_to
                        current_cluster_weight_out -= remove_candidates[best_change].sum_weight_from

                        # Remove the change vertex from remove_candidates and current_cluster
                        to_remove = remove_candidates[best_change].copy()
                        del remove_candidates[best_change]
                        del current_cluster[best_change]

                        # Also add the change vertex to add_candidates
                        # TODO: figure out if we should put a removed vertex back in add_candidates
                        #   I think it makes sense to do.  Check with original authors.
                        add_candidates[best_change] = to_remove

                        ############################################################
                        # iterate over neighbors of best_change, and update each Relationship
                        ############################################################
                        for v in self.graph.hash_graph[best_change]:
                            edge_weight = self.graph.hash_graph[best_change][v]
                            # note that v may be in both the current_cluster and in remove_candidates
                            if v in remove_candidates:
                                remove_candidates[v].sum_weight_to = -1 * edge_weight + remove_candidates[
                                    v].sum_weight_to
                                remove_candidates[v].sum_weight_from = edge_weight + remove_candidates[
                                    v].sum_weight_from
                                remove_candidates[v].num_edges_to = -1 + remove_candidates[v].num_edges_to
                                remove_candidates[v].num_edges_from = 1 + remove_candidates[v].num_edges_from
                            if v in current_cluster:
                                current_cluster[v].sum_weight_to = -1 * edge_weight + current_cluster[v].sum_weight_to
                                current_cluster[v].sum_weight_from = edge_weight + current_cluster[v].sum_weight_from
                                current_cluster[v].num_edges_to = -1 + current_cluster[v].num_edges_to
                                current_cluster[v].num_edges_from = 1 + current_cluster[v].num_edges_from
                                # TODO: if v is in the current cluster, and is newly also on the boundary, add v to remove_candidates
                                if current_cluster[v].num_edges_from == 1:
                                    remove_candidates[v] = current_cluster[v].copy()
                            if v in add_candidates:
                                # update the relationship
                                add_candidates[v].sum_weight_to = -1 * edge_weight + add_candidates[v].sum_weight_to
                                add_candidates[v].sum_weight_from = edge_weight + add_candidates[v].sum_weight_from
                                add_candidates[v].num_edges_to = - 1 + add_candidates[v].num_edges_to
                                add_candidates[v].num_edges_from = 1 + add_candidates[v].num_edges_from
                                # TODO: if v in no longer on the boundary after best_change is removed from current_cluster, remove v from add_candidates
                                if add_candidates[v].num_edges_to == 0:
                                    del add_candidates[v]
                    else:
                        debug("\n", "No improvement by REMOVING", "\n")

                sleep_debug(1)

            # add current_cluster to the list of clusters
            self.initial_clustering.append(current_cluster)
            index += 1
            if not self.seed_from_all:
                for v in current_cluster:
                    considered_vertices.add(v)
            print("CLUSTER #%s: %s" % (str(len(self.initial_clustering)), str([vertex for vertex in current_cluster])))
            # time.sleep(2)
