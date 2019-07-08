exec(open("/home/alonzo/PycharmProjects/relinking/common.py").read())
import numpy
import sys
import time
DEBUG = False
SLEEP_DEBUG = False
def debug(*argv):
    if DEBUG:
        for arg in argv:
            print(arg)

def sleep_debug(t):
    if SLEEP_DEBUG:
        time.sleep(t)

class Relationship:
    def __init__(self, sum_weight_to:float, num_edges_to:int, sum_weight_from:float, num_edges_from:int):
        self.sum_weight_to = sum_weight_to
        self.num_edges_to = num_edges_to
        self.sum_weight_from = sum_weight_from
        self.num_edges_from = num_edges_from

    def copy(self):
        return Relationship(self.sum_weight_to, self.num_edges_to, self.sum_weight_from, self.num_edges_from)

class CL1_Randomized:

    def __init__(self, base_file_path, original_graph_filename, quality_function_name):
        self.graph = Graph(base_file_path+"/"+original_graph_filename)
        self.vertices_by_degree = self.sort_vertices_by_degree()
        self.quality_function_name = quality_function_name
        self.clustering = list()

    def sort_vertices_by_degree(self):
        retval = [[k, len(self.graph.hash_graph[k])] for k in self.graph.hash_graph]
        retval.sort(key=lambda x: x[1], reverse=True)
        return retval

    def get_clusters(self):
        clusters = []
        for cluster in self.clustering:
            clusters.append([vertex for vertex in cluster])
        return clusters

    def print_translated_clusters(self):
        clusters = self.get_clusters()
        for cluster in clusters:
            s = ""
            for id in cluster:
                s+=self.graph.id_to_name[id]+"\t"
            print(s)


    def original_construction(self):

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
                debug("SKIP %s"%str(current_seed))
                index+=1
                continue
            else:
                debug("Starting cluster #%s"%str(len(self.clustering)))
                # time.sleep(3)
                # TODO: ignore vertices that have been removed before during construction of current cluster
                ignore_vertices = set()

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
                    num_edges_from = num_edges - num_edges_from
                    weight_from = sum([self.graph.hash_graph[target][tar] for tar in self.graph.hash_graph[target] if tar != current_seed])
                    add_candidates[target] = Relationship(weight_to, num_edges_to, weight_from, num_edges_from)


                improvement_flag=True

                while improvement_flag and (add_candidates or remove_candidates):
                    improvement_flag = False

                    # Consider adding a vertex on the boundary
                    #
                    #
                    #
                    best_change = None
                    best_change_score = current_score
                    for v in add_candidates:
                        numerator = current_cluster_weight_in + add_candidates[v].sum_weight_to
                        denominator = current_cluster_weight_in + current_cluster_weight_out + add_candidates[v].sum_weight_to + add_candidates[v].sum_weight_from
                        proposed_score = numerator / denominator

                        if proposed_score>best_change_score:
                            best_change = v
                            best_change_score = proposed_score
                        debug("##################### ADD Consideration ########################")
                        debug("v: %s"%str(v))
                        debug("proposed_score: %s"% str(proposed_score))
                        debug("best_change_score: %s"%str(best_change_score))
                        debug("best_change: %s"%str(best_change))
                        debug("numerator: %s" %str(numerator))
                        debug("denominator: %s" %str(denominator))
                        debug("current_cluster_weight_in: %s"%str(current_cluster_weight_in))
                        debug("add_candidates[v].sum_weight_to: %s"%str(add_candidates[v].sum_weight_to))
                        debug("current_cluster_weight_out: %s"%str(current_cluster_weight_out))
                        debug("add_candidates[v].sum_weight_from: %s"%str(add_candidates[v].sum_weight_from))
                        sleep_debug(.25)

                    if best_change:
                        debug("\n", "ADD: %s" %str(best_change), "best_change_score: %s"%str(best_change_score), "\n")
                        improvement_flag = True

                        # update the current_cluster 's score
                        current_score = best_change_score
                        # update the overall weight into and out of the current_cluster
                        current_cluster_weight_in += add_candidates[best_change].sum_weight_to
                        current_cluster_weight_out += add_candidates[best_change].sum_weight_from

                        # Move the change vertex from add_candidates to current_cluster
                        to_add = add_candidates[best_change].copy()
                        current_cluster[best_change] = to_add
                        del add_candidates[best_change]

                        # Also add the change vertex to remove_candidates if applicable
                        if to_add.num_edges_from:
                            remove_candidates[best_change] = to_add

                        for v in self.graph.hash_graph[best_change]:     # iterate over neighbors of best_change, and update each Relationship
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
                            # note that v may be in both the current_cluster and in remove_candidates
                            # remove_candidates is a subset of current_cluster
                            if v in remove_candidates:
                                remove_candidates[v].sum_weight_to = edge_weight + remove_candidates[v].sum_weight_to
                                remove_candidates[v].sum_weight_from = -1 * edge_weight + remove_candidates[v].sum_weight_from
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
                                        num_edges_to+=1
                                        weight_to+=weight_prime
                                    else:
                                        num_edges_from+=1
                                        weight_from+=weight_prime
                                add_candidates[v] = Relationship(weight_to, num_edges_to, weight_from, num_edges_from)

                    else:
                        debug("\n","No improvement by ADDING", "\n")

                    # sleep_debug(1)

                    # Consider removing a vertex on the boundary
                    #
                    #
                    #
                    # check that
                    #   (1) the cluster has more than one element
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
                            for potential_start_point in current_cluster_membership_hashset:
                                if potential_start_point != v:
                                    start_point = potential_start_point
                                    # consider break statement here
                            # check that
                            #   (2) removal of vertex under consideration will not disconnect cluster
                            debug("DFS starting vertex: %s"%str(start_point), "DFS ignore vertex: %s"%str(v))
                            dfs(start_point, v,  current_cluster_membership_hashset, visted)
                            debug("=============================")
                            if len(visted) == -1 + len(current_cluster_membership_hashset):
                                is_a_cut = False
                                debug("%s is NOT a CUT"%str(v))
                            if is_a_cut:
                                debug("%s is a CUT!" % str(v))
                            debug("cluster: %s"% str(current_cluster_membership_hashset))
                            debug("visited by DFS: %s"%str(visted))
                            sleep_debug(.25)
                            if not is_a_cut:
                                # TODO: check that this makes sense
                                numerator = current_cluster_weight_in - remove_candidates[v].sum_weight_to
                                denominator = current_cluster_weight_in + current_cluster_weight_out - remove_candidates[v].sum_weight_to - remove_candidates[v].sum_weight_from
                                proposed_score = numerator / denominator
                                if proposed_score>best_change_score:
                                    best_change = v
                                    best_change_score = proposed_score
                                debug("##################### REMOVE Consideration ########################")
                                debug("v: %s"%str(v))
                                debug("proposed_score: %s" % str(proposed_score))
                                debug("best_change_score: %s" % str(best_change_score))
                                debug("best_change: %s" % str(best_change))
                                debug("numerator: %s" % str(numerator))
                                debug("denominator: %s" % str(denominator))
                                debug("current_cluster_weight_in: %s" % str(current_cluster_weight_in))
                                debug("remove_candidates[v].sum_weight_to: %s" % str(remove_candidates[v].sum_weight_to))
                                debug("current_cluster_weight_out: %s" % str(current_cluster_weight_out))
                                debug("remove_candidates[v].sum_weight_from: %s" % str(remove_candidates[v].sum_weight_from))
                                sleep_debug(1)
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
                            add_candidates[best_change] = to_remove

                            # AFTER this is done, THEN do the following
                            for v in self.graph.hash_graph[best_change]:  # iterate over neighbors of best_change, and update each Relationship
                                edge_weight = self.graph.hash_graph[best_change][v]
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
                                    # TODO: if v in no longer on the boundary after best_change is removed from current_cluster, remove v from add_candidates
                                    if add_candidates[v].num_edges_to == 0:
                                        del add_candidates[v]
                        else:
                            debug("\n", "No improvement by REMOVING", "\n")

                    sleep_debug(1)

                #add current_cluster to the list of clusters
                self.clustering.append(current_cluster)
                index += 1
                for v in current_cluster:
                    considered_vertices.add(v)
                print("CLUSTER #%s: %s"%(str(len(self.clustering)), str([vertex for vertex in current_cluster])))
                # time.sleep(2)


    def merger(self, clusters):# takes a list of clusters
        num_clusters = len(clusters)
        for i in range(num_clusters):
            for j in range(num_clusters):
                print("Dummy merger")








    # def quality_function(self, current_cluster:dict, change_vertex = None, change_relationship=None, action=""):
    #     if action:
    #         if action is 'add':
    #             dummy = -1
    #         if action is 'remove':
    #             dummy = -2
    #     else:
    #         dummy = -3
    #     return dummyPycharmProjects


    # def change_add(self,  current_cluster:dict, change_vertex = None, change_relationship=None, action =""):
    #     return "dummy"


    # def quality_function_naive(self, current_cluster:dict, change_vertex=None, change_relationship=None, action =""):
    #     weight_in = 0
    #     weight_out = 0
    #     for source in current_cluster:
    #         for target in self.graph.hash_graph[source]:
    #             edge_weight = self.graph.hash_graph[source][target]
    #             if target in current_cluster:
    #                 weight_in+=edge_weight
    #             else:
    #                 weight_out+=edge_weight
    #     weight_in /= 2
    #     cohesiveness = weight_in / (weight_in + weight_out)


