from common import *
import numpy
import sys
import math
#TODO: make the random proportional bad_adds option work for randomized_construction()
#TODO: check if Cl1 are always connected, in original implementation, and in my implementation
#TODO: check if Cl1 found are similar to my found
#TODO: figure out why found and not_found have added cardinality larger than gsc_appearing_in_dataset when the threshold is less than 1
#TODO: figure out how sensitive the clusters are to the choice of seeds
#TODO: figure out which complexes are missed by my algo



def jaccard_similarity(l1:list, l2:list)->float:
    set1 = set(l1)
    set2 = set(l2)
    numerator = len(set1.intersection(set2))
    denominator = len(set1.union(set2))
    return numerator/denominator

class CL1_Randomized:
    def __init__(self,
                 base_file_path,
                 original_graph_filename,
                 quality_function_name,
                 output_filename="unnamed.txt",
                 density_threshold = .3,
                 merge_threshold = .8,
                 penalty_value_per_node = 2,
                 randomized_construction_bool= False,
                 number_of_shakes = 0,
                 number_of_bad_adds = 0,
                 bad_add_probability = 0,
                 add_with_proportional_probability = False,
                 sort_seeds_by="degree",
                 care_about_cuts=True,
                 seed_from_all = False,
                 gsc_appearance_ratio_threshold=1,
                 found_gsc_jaccard_threshold = 1,
                 gold_standard_filename = ""):
        self.base_file_path = base_file_path
        self.graph = Graph(base_file_path+"/"+original_graph_filename)
        self.vertices_by_degree = self.sort_vertices_by_degree()
        self.vertices_by_weight = self.sort_vertices_by_weight()
        self.quality_function_name = quality_function_name
        self.output_filename = output_filename
        self.density_threshold = density_threshold
        self.merge_threshold = merge_threshold
        self.penalty_value_per_node = penalty_value_per_node
        self.randomized_construction_bool = randomized_construction_bool
        self.number_of_shakes = number_of_shakes
        self.number_of_bad_adds = number_of_bad_adds
        self.bad_add_probability = bad_add_probability
        # TODO implement the feature of adding with proportional probability
        self.add_with_proportional_probability = add_with_proportional_probability
        self.sort_seeds_by = sort_seeds_by
        self.care_about_cuts = care_about_cuts
        self.seed_from_all = seed_from_all
        self.gsc_appearance_ratio_threshold = gsc_appearance_ratio_threshold
        self.found_gsc_jaccard_threshold = found_gsc_jaccard_threshold
        self.gold_standard_filename = gold_standard_filename

        # the gold standard complexes whose entire set of proteins appear in the dataset
        # TODO: prune gold_standard_complexes_appearing_in_dataset if it is unlikely that is "SHOULD" be detected by algorithm
        self.gold_standard_complexes_appearing_in_dataset = list()
        self.gold_standard_complexes = list()
        self.clustering = list()
        self.cluster_list = list()
        self.merged_cluster_list = list()
        self.sizeThreshold_merged_cluster_list = list()
        self.densityThreshold_sizeThreshold_merged_cluster_list = list()
        # the gold standard complexes appearing in the dataset that were  found by the algorithm
        self.found = []
        # the gold standard complexes appearing in the dataset that were not found by the algorithm
        self.not_found = []

    def sort_vertices_by_degree(self):
        retval = [[k, len(self.graph.hash_graph[k])] for k in self.graph.hash_graph]
        retval.sort(key=lambda x: x[1], reverse=True)
        return retval

    def sort_vertices_by_weight(self):
        retval = [[k, sum([self.graph.hash_graph[k][target] for target in self.graph.hash_graph[k]])] for k in self.graph.hash_graph]
        retval.sort(key=lambda x: x[1], reverse=True)
        return retval

    def get_clusters(self):
        clusters = []
        for cluster in self.clustering:
            clusters.append([vertex for vertex in cluster])
        return clusters

    def make_cluster_list(self):
        self.cluster_list = [set(cluster) for cluster in self.get_clusters()]

    def write_final_clusters(self):
        f = open(self.base_file_path+"/"+ self.output_filename, "w+")
        counter = 1
        for cluster in self.densityThreshold_sizeThreshold_merged_cluster_list:
            s = ""
            for id in cluster:
                s+=self.graph.id_to_name[id]+"\t"
            print("Cluster #%s: "%str(counter),len(cluster), " proteins")
            print(s)
            f.write(s + "\n")
            counter+=1
        f.close()

    def get_quality(self):
        import subprocess
        res = subprocess.check_output(["python2",
                                       "cl1_reproducibility/reproducibility/scripts/match_standalone.py",
                                       self.gold_standard_filename,
                                       self.base_file_path +"/"+self.output_filename])
        for line in res.splitlines():
            print(line)

    def sizeThreshold(self):
        self.sizeThreshold_merged_cluster_list = [cluster for cluster in self.merged_cluster_list if len(cluster) > 2]

    def densityThreshold(self):

        def checkDensity(cluster_set):
            in_weight=0
            for source in cluster_set:
                for target in self.graph.hash_graph[source]:
                    if target in cluster_set:
                        in_weight+=self.graph.hash_graph[source][target]
            n = len(cluster_set)
            # TODO: check that authors didn't mean n* (n+1)/2
            denominator = (n * (n-1)) /2
            return in_weight/denominator

        for cluster_set in self.sizeThreshold_merged_cluster_list:
            density = checkDensity(cluster_set)
            if density>self.density_threshold:
                self.densityThreshold_sizeThreshold_merged_cluster_list.append(cluster_set)
        return "Dummy"

    def process(self):
        if self.randomized_construction_bool:
            self.randomized_construction()
        else:
            self.original_construction()
        self.make_cluster_list()
        self.merger()
        self.sizeThreshold()
        self.densityThreshold()
        self.write_final_clusters()
        self.gold_standard_complex_appearance()
        self.found_gsc()
        self.not_found_gsc()
        print("################ FOUND AND UNFOUND GSC DETAILS #######################")
        self.found_and_unfound_details()
        print("################ QUALITY #######################")
        self.get_quality()

    def gold_standard_complex_appearance(self,gold_standard_filename=""):
        if gold_standard_filename:
            gsf = gold_standard_filename
        else:
            gsf = self.gold_standard_filename

        f = open(gsf, "r")
        li = []
        counter = 0
        for line in f.readlines():
            complex = line.split()
            li.append(complex)
            counter+=1
        self.gold_standard_complexes = li
        intermediate = []
        for complex in li:
            complex_appears = True
            proteins_from_current_gsc_found_in_dataset = []
            for protein in complex:
                if protein not in self.graph.name_to_id:
                    complex_appears=False
                else:
                    proteins_from_current_gsc_found_in_dataset.append(protein)
            if complex_appears:
                intermediate.append(complex)
            elif len(proteins_from_current_gsc_found_in_dataset)/float(len(complex)) >= self.gsc_appearance_ratio_threshold and len(proteins_from_current_gsc_found_in_dataset)>=3:
                intermediate.append(proteins_from_current_gsc_found_in_dataset)


        final=[]
        for complex in intermediate:
            current = []
            for protein in complex:
                current.append(self.graph.name_to_id[protein])
            final.append(current)
        self.gold_standard_complexes_appearing_in_dataset = final
        return final

    def found_gsc(self):
        retval = []
        for A in self.gold_standard_complexes_appearing_in_dataset:
            for B in self.densityThreshold_sizeThreshold_merged_cluster_list:
                if jaccard_similarity(A,B) >= self.found_gsc_jaccard_threshold:
                    retval.append(A)
        self.found = retval
        return retval

    def not_found_gsc(self):
        retval = []
        for A in self.gold_standard_complexes_appearing_in_dataset:
            found = False
            for B in self.densityThreshold_sizeThreshold_merged_cluster_list:
                if jaccard_similarity(A,B) >= self.found_gsc_jaccard_threshold:
                    found = True
            if not found:
                retval.append(A)
        self.not_found = retval
        return retval

    def found_and_unfound_details(self):
        print("FOUND")
        cohesiveness_tot = 0
        density_tot = 0
        length_tot = 0
        for c in self.found:
            c_cohesiveness = self.cohesiveness(c)
            cohesiveness_tot += c_cohesiveness
            c_density = self.density(c)
            density_tot += c_density
            length_tot += len(c)
            print(c)
            print(len(c), c_cohesiveness, c_density, len(c))
        c1 = cohesiveness_tot / float(len(self.found))
        d1 = density_tot / float(len(self.found))
        l1 = length_tot / float(len(self.found))
        print("--------------------------------------")

        print("NOT FOUND")
        cohesiveness_tot = 0
        density_tot = 0
        length_tot = 0
        for c in self.not_found:
            c_cohesiveness = self.cohesiveness(c)
            cohesiveness_tot += c_cohesiveness
            c_density = self.density(c)
            density_tot += c_density
            length_tot += len(c)
            print(c)
            print(len(c), c_cohesiveness, c_density, len(c))
        c2 = cohesiveness_tot / float(len(self.not_found))
        d2 = density_tot / float(len(self.not_found))
        l2 = length_tot / float(len(self.not_found))
        print("--------------------------------------")
        print("GSC appearing in dataset: ", len(self.gold_standard_complexes_appearing_in_dataset))
        print("found: ", len(self.found))
        print("not found: ", len(self.not_found))
        print("average cohesiveness of found= ", c1)
        print("average density of found= ", d1)
        print("average length of found= ", l1)
        print("average cohesiveness of NOT found= ", c2)
        print("average density of NOT found= ", d2)
        print("average length of NOT found= ", l2)
        print(len(self.gold_standard_complexes), " reference complexes")
        print(len(self.gold_standard_complexes_appearing_in_dataset), " appear in the dataset")

    def cohesiveness(self, list_of_proteins):
        weight_in = 0
        weight_out = 0
        for source in list_of_proteins:
            for target in self.graph.hash_graph[source]:
                if target in list_of_proteins:
                    weight_in += self.graph.hash_graph[source][target] /2.0
                else:
                    weight_out += self.graph.hash_graph[source][target]
        return weight_in/(weight_in + weight_out)

    def density(self, list_of_proteins):
        in_weight = 0
        for source in list_of_proteins:
            for target in self.graph.hash_graph[source]:
                if target in list_of_proteins:
                    in_weight += self.graph.hash_graph[source][target]
        n = len(list_of_proteins)
        # TODO: check that authors didn't mean n* (n+1)/2
        denominator = (n * (n - 1)) / 2
        return in_weight / denominator

    def modularity(self, list_of_proteins):
        return 1

    def randomized_construction(self):

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

        if self.sort_seeds_by == 'degree':
            sorted_seeds = self.vertices_by_degree
        else:
            sorted_seeds = self.vertices_by_weight

        while index < len(sorted_seeds):
            current_seed = sorted_seeds[index][0]
            current_seed_degree = sorted_seeds[index][1]
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


                last_failed_add_round = -777
                last_failed_remove_round = -666
                round = 0

                def find_best_suboptimal_add(cc_weight_in, cc_weight_out):
                    best_suboptimal_change = None
                    best_suboptimal_change_score = -10000
                    for v in add_candidates:
                        numerator = cc_weight_in + add_candidates[v].sum_weight_to
                        denominator = cc_weight_in + cc_weight_out + add_candidates[v].sum_weight_from + self.penalty_value_per_node * (len(current_cluster)+1)
                        proposed_score = numerator / denominator
                        if proposed_score>best_change_score:
                            best_suboptimal_change = v
                            best_suboptimal_change_score = proposed_score
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
                        debug("len(current_cluster): %s"%str(len(current_cluster)))
                        sleep_debug(.25)
                    return best_suboptimal_change, best_suboptimal_change_score

                def find_best_add():
                    best_change = None
                    best_change_score = current_score
                    for v in add_candidates:
                        numerator = current_cluster_weight_in + add_candidates[v].sum_weight_to
                        denominator = current_cluster_weight_in + current_cluster_weight_out + add_candidates[v].sum_weight_from + self.penalty_value_per_node * (len(current_cluster)+1)
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
                        debug("len(current_cluster): %s"%str(len(current_cluster)))
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
                    # TODO: figure out if we should put a removed vertex back in add_candidates
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

                def add_shake(cc_weight_in, cc_weight_out, round, last_failed_add_round):
                    for i in range(self.number_of_bad_adds):
                        best_suboptimal_change, best_suboptimal_score = find_best_suboptimal_add(cc_weight_in, cc_weight_out)
                        if best_suboptimal_change:
                            debug("suboptimal add")
                            sleep_debug(1)
                            round += 1
                            last_failed_add_round = -5
                            # TODO handle round numbers
                            cc_weight_in, cc_weight_out = add(best_suboptimal_change,
                                                                                        best_suboptimal_score,
                                                                                        cc_weight_in,
                                                                                        cc_weight_out)
                    return best_suboptimal_score, cc_weight_in, cc_weight_out, round, last_failed_add_round


                local_number_of_shakes_remaining = self.number_of_shakes
                while (add_candidates or remove_candidates) and abs(last_failed_remove_round-last_failed_add_round) != 1:
                    debug("Current cluster #%s" % str(len(self.clustering)))
                    decider = numpy.random.rand()
                    # Consider ADDING a vertex on the boundary
                    #
                    if (decider <= .5 or last_failed_remove_round == round) and last_failed_add_round!=round:
                        round+=1
                        best_change, best_change_score = find_best_add()
                        if best_change:
                            current_score = best_change_score
                            current_cluster_weight_in, current_cluster_weight_out = add(best_change, best_change_score, current_cluster_weight_in, current_cluster_weight_out)
                        else:
                            debug("\n","No improvement by ADDING", "\n")
                            last_failed_add_round = round

                    # Consider REMOVING a vertex on the boundary
                    #
                    if (decider>.5 or last_failed_add_round == round) and  last_failed_remove_round!=round:
                        round+=1
                        # check that
                        #   (1) the cluster has more than one element
                        debug("length of current cluster: ", len(current_cluster))
                        best_change, best_change_score = find_best_remove(current_score)
                        if best_change:
                            current_cluster_weight_in, current_cluster_weight_out = remove(best_change, best_change_score, current_cluster_weight_in, current_cluster_weight_out)
                        else:
                            debug("\n", "No improvement by REMOVING, len(current_cluster) = 1", "\n")
                            last_failed_remove_round = round

                    # If stuck in local optimum, consider taking a 'bad' add step
                    #
                    if local_number_of_shakes_remaining and add_candidates and abs(last_failed_remove_round - last_failed_add_round) == 1:
                        local_number_of_shakes_remaining -= 1
                        current_score, current_cluster_weight_in, current_cluster_weight_out, round, last_failed_add_round = \
                            add_shake(current_cluster_weight_in, current_cluster_weight_out, round, last_failed_add_round)

                    debug("$$$$$$$$$", last_failed_add_round, last_failed_remove_round, decider)



                #add current_cluster to the list of clusters
                self.clustering.append(current_cluster)
                index += 1
                if not self.seed_from_all:
                    for v in current_cluster:
                        considered_vertices.add(v)
                print("CLUSTER #%s: %s"%(str(len(self.clustering)), str([vertex for vertex in current_cluster])))
                print(last_failed_add_round, last_failed_remove_round)
                # time.sleep(2)

    def merger(self):# takes a list of clusters
        threshold = self.merge_threshold
        indices_of_considered_clusters = {}
        hash_graph = dict()

        def similarity(A, B):
            numerator = len(A.intersection(B))**2
            denominator = len(A) * len(B)
            return numerator / denominator

        def dfs(index, local_visited):
            local_visited.add(index)
            for neigbor in hash_graph[index]:
                if neigbor not in local_visited:
                    dfs(neigbor, local_visited)

        for i in range(len(self.cluster_list)):
            if i not in hash_graph:
                hash_graph[i]=set()
            for j in range(i+1, len(self.cluster_list)):
                if similarity(self.cluster_list[i], self.cluster_list[j])>threshold:
                    if i in hash_graph:
                        hash_graph[i].add(j)
                    else:
                        hash_graph[i] = {j}
                    if j in hash_graph:
                        hash_graph[j].add(i)
                    else:
                        hash_graph[j] = {i}

        global_visited = set()
        merge_indices = list()
        for index in range(len(self.cluster_list)):
            if index not in global_visited:
                local_visited = set()
                dfs(index, local_visited)
                # TODO: is the copy necessary?
                merge_indices.append(local_visited.copy())
                for reached in local_visited:
                    global_visited.add(reached)


        new_clusters = list()
        for group in merge_indices:
            new_cluster = set()
            for cluster_index in group:
                new_cluster = new_cluster.union(self.cluster_list[cluster_index])
            new_clusters.append(new_cluster)

        self.merged_cluster_list = new_clusters

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

                    # Consider ADDING a vertex on the boundary
                    #
                    #
                    #
                    best_change = None
                    best_change_score = current_score
                    for v in add_candidates:
                        numerator = current_cluster_weight_in + add_candidates[v].sum_weight_to
                        denominator = current_cluster_weight_in + current_cluster_weight_out + add_candidates[v].sum_weight_from + self.penalty_value_per_node * (len(current_cluster)+1)

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

                    # Consider REMOVING a vertex on the boundary
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
                                denominator = current_cluster_weight_in + current_cluster_weight_out - remove_candidates[v].sum_weight_from + self.penalty_value_per_node * (len(current_cluster)-1)
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
                if not self.seed_from_all:
                    for v in current_cluster:
                        considered_vertices.add(v)
                print("CLUSTER #%s: %s"%(str(len(self.clustering)), str([vertex for vertex in current_cluster])))
                # time.sleep(2)



