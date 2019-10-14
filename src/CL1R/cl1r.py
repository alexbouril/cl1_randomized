from src.COMMON.cmn import *
from src.GRAPH.graph import *
from src.CONSTRUCTION.randomized_construction import randomized_construction
from src.CONSTRUCTION.original_construction import original_construction
from src.QUALITY.quality import density, cohesiveness
from src.MERGING.merge import merge
from src.THRESHOLDING.thresholding_functions import sizeThreshold, densityThreshold

# TODO: package/syspath,1-1 comparison, Cython, Threading, Command-line args, currying, class variable, hidden class variable, numpy, inheritance, tracemalloc



#TODO implement find 1, 2, 3 neighborhood of current cluster
#TODO implement add based on cohesiveness of 2 neighborhood
#TODO GRAPH cluster with their 1,2,3 neighborhoods in the background
#TODO Fix object saving

#TODO try MERGING a run of add_shake enabled with a run of add_shake disabled


#TODO: make the random proportional bad_adds option work for randomized_construction()
#TODO: check if Cl1 are always connected, in original implementation, and in my implementation
#TODO: check if Cl1 found are similar to my found
#TODO: figure out why found and not_found have added cardinality larger than gsc_appearing_in_dataset when the threshold is less than 1
#TODO: figure out how sensitive the clusters are to the choice of seeds
#TODO: figure out which complexes are missed by my algo
#TODO: multi-pass merger


class CL1_Randomized:
    def __init__(self,
                 base_file_path,
                 original_graph_filename,
                 quality_function_name,
                 density_threshold = .3,
                 merge_threshold = .8,
                 penalty_value_per_node = 2,
                 randomized_construction_bool= False,
                 rng_seed = None,
                 number_of_shakes = 0,
                 number_of_bad_adds = 0,
                 bad_add_probability = 0,
                 add_with_proportional_probability = False,
                 sort_seeds_by="degree",
                 care_about_cuts=True,
                 seed_from_all = False,
                 gsc_appearance_ratio_threshold=1,
                 gsc_appearance_density_threshold = 0,
                 found_gsc_jaccard_threshold = 1,
                 gold_standard_filename = ""):
        self.time_of_run = str(datetime.datetime.now()).replace(" ","_").replace(".",":")
        self.run_title = original_graph_filename.replace(".txt", "")+"+"+self.time_of_run
        self.argument_dict = locals()

        ############################################################
        # STORE ARGUMENTS
        ############################################################
        self.base_file_path = base_file_path
        self.graph = Graph(base_file_path+"/"+original_graph_filename)

        print("size of GRAPH in bytes: %s"%str(sys.getsizeof(self.graph.hash_graph)))
        print("number of nodes: %s"%str(self.graph.num_proteins))
        print("number of edges: %s"%str(self.graph.num_edges))

        self.vertices_by_degree = sort_vertices_by_degree(self)
        self.vertices_by_weight = sort_vertices_by_weight(self)
        self.quality_function_name = quality_function_name
        self.output_filename = "complexes+" + self.run_title+".txt"
        self.density_threshold = density_threshold
        self.merge_threshold = merge_threshold
        self.penalty_value_per_node = penalty_value_per_node
        self.randomized_construction_bool = randomized_construction_bool
        self.rng_seed = rng_seed
        numpy.random.seed(rng_seed)
        self.rng_initial_state = numpy.random.get_state()
        self.number_of_shakes = number_of_shakes
        self.number_of_bad_adds = number_of_bad_adds
        self.bad_add_probability = bad_add_probability
        # TODO implement the feature of adding with proportional probability
        self.add_with_proportional_probability = add_with_proportional_probability
        self.sort_seeds_by = sort_seeds_by
        self.care_about_cuts = care_about_cuts
        self.seed_from_all = seed_from_all
        self.gsc_appearance_ratio_threshold = gsc_appearance_ratio_threshold
        # TODO implement the density threshold for gsc appearance
        self.gsc_appearance_density_threshold = gsc_appearance_density_threshold
        self.found_gsc_jaccard_threshold = found_gsc_jaccard_threshold
        self.gold_standard_filename = gold_standard_filename

        ############################################################
        # SET UP THE LOGGER
        ############################################################
        # self.logger = setup_custom_logger('CL1R')
        # self.logger.info(f'(RUN:{self.run_title}) Starting CL1R v{__version__} ({__status__}) [{__website__}]')

        ############################################################
        # OUTPUTS OF CALCULATIONS
        ############################################################
        # TODO: prune gold_standard_complexes_appearing_in_dataset if it is unlikely that is "SHOULD" be detected by algorithm
        # the gold standard complexes whose entire set of proteins appear in the dataset
        self.gold_standard_complexes_appearing_in_dataset = list()
        self.gold_standard_complexes = list()
        self.initial_clustering = list()
        self.initial_clustering_seeds = list()
        self.initial_cluster_list = list()
        self.zipper_merged_cluster_list = list()
        self.merged_cluster_list = list()
        self.sizeThreshold_merged_cluster_list = list()
        self.densityThreshold_sizeThreshold_merged_cluster_list = list()
        # the gold standard complexes appearing in the dataset that were  found by the algorithm
        self.found = []
        # the gold standard complexes appearing in the dataset that were not found by the algorithm
        self.not_found = []
        self.quality_report = {}
        self.final_clusters_stats = dict()
        self.gsc_appearing_stats = dict()
        self.gsc_appearing_found_stats = dict()
        self.gsc_appearing_notFound_stats = dict()

        ############################################################
        # STORE INFORMATION ABOUT STATES OF CLUSTER DURING CONSTRUCTION
        ############################################################
        self.construction_log = dict()

        ############################################################
        # RUN THE ALGORITHM
        ############################################################
        self.process()


    def process(self):
        ############################################################
        # INITIAL CONSTRUCTION
        ############################################################
        if self.randomized_construction_bool:
            randomized_construction(self)
        # TODO reimplement the oc function using new classes and package framework
        else:
            original_construction(self)

        ############################################################
        # GET THE INITIAL CLUSTER LIST
        ############################################################
        #TODO document better what the point of this step is
        def make_initial_cluster_list():
            # TODO: rid redundant work
            def get_initial_clusters():
                clusters = []
                for cluster in self.initial_clustering:
                    clusters.append([vertex for vertex in cluster])
                return clusters

            self.initial_cluster_list = [set(cluster) for cluster in get_initial_clusters()]
        make_initial_cluster_list()

        ############################################################
        # MERGE THE INITIAL CLUSTERS
        ############################################################
        self.merged_cluster_list+= merge(merge_threshold=self.merge_threshold,
                                         source=self.initial_cluster_list)
        ############################################################
        # ZIPPER MERGE COMPLEXES
        ############################################################
        def zipper_merger():  # takes a list of clusters
            threshold = self.merge_threshold
            hash_graph = dict()

            def find_one_neighborhood_of_complex(complex):
                complex_set = set([v for v in complex])
                neighbor_set = set()
                for v in complex_set:
                    for u in self.graph.hash_graph[v]:
                        if u not in complex_set:
                            neighbor_set.add(u)
                return complex_set, neighbor_set

            def similarity(A, B):
                A_complex_set, A_neighbor_set = find_one_neighborhood_of_complex(A)
                B_complex_set, B_neighbor_set = find_one_neighborhood_of_complex(B)



                """Implements the overlap score described in the paper

                :param A: a set
                :param B: a set
                :return: the overlap score
                """
                numerator = len(A.intersection(B)) ** 2
                denominator = len(A) * len(B)
                return numerator / denominator

            def dfs(index, local_visited):
                local_visited.add(index)
                for neigbor in hash_graph[index]:
                    if neigbor not in local_visited:
                        dfs(neigbor, local_visited)

            for i in range(len(self.initial_cluster_list)):
                if i not in hash_graph:
                    hash_graph[i] = set()
                for j in range(i + 1, len(self.initial_cluster_list)):
                    if similarity(self.initial_cluster_list[i], self.initial_cluster_list[j]) > threshold:
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
            for index in range(len(self.initial_cluster_list)):
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
                    new_cluster = new_cluster.union(self.initial_cluster_list[cluster_index])
                new_clusters.append(new_cluster)

            self.merged_cluster_list = new_clusters

        # merger()



        ############################################################
        # THRESHOLD BASED ON SIZE
        ############################################################
        self.sizeThreshold_merged_cluster_list+=sizeThreshold(source=self.merged_cluster_list,
                                                              size_thresh=2)
        ############################################################
        # THRESHOLD BASED ON DENSITY
        ############################################################
        self.densityThreshold_sizeThreshold_merged_cluster_list +=densityThreshold(hash_graph=self.graph.hash_graph,
                                                                                   li=self.sizeThreshold_merged_cluster_list,
                                                                                   dens_thresh=self.density_threshold)
        ############################################################
        # WRITE THE FINAL RESULT INTO A TEXT FILE
        ############################################################
        def write_final_clusters():
            f = open("../complexes/" + self.output_filename, "w+")
            counter = 1
            for cluster in self.densityThreshold_sizeThreshold_merged_cluster_list:
                s = ""
                for id in cluster:
                    s += self.graph.id_to_name[id] + "\t"
                print("Cluster #%s: " % str(counter), len(cluster), " proteins")
                print(s)
                f.write(s + "\n")
                counter += 1
            f.close()

        write_final_clusters()


        ############################################################
        # DETERMINE WHICH GOLD STANDARD COMPLEXES APPEAR IN CURRENT DATASET
        ############################################################
        def gold_standard_complex_appearance(gold_standard_filename=""):
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
                counter += 1
            self.gold_standard_complexes = li
            intermediate = []
            for complex in li:
                complex_appears = True
                proteins_from_current_gsc_found_in_dataset = []
                for protein in complex:
                    if protein not in self.graph.name_to_id:
                        complex_appears = False
                    else:
                        proteins_from_current_gsc_found_in_dataset.append(protein)
                if complex_appears:
                    intermediate.append(complex)
                elif len(proteins_from_current_gsc_found_in_dataset) / float(
                        len(complex)) >= self.gsc_appearance_ratio_threshold and len(
                        proteins_from_current_gsc_found_in_dataset) >= 3:
                    intermediate.append(proteins_from_current_gsc_found_in_dataset)
            final = []
            for complex in intermediate:
                current = []
                for protein in complex:
                    current.append(self.graph.name_to_id[protein])
                final.append(current)
            self.gold_standard_complexes_appearing_in_dataset = final
            return final

        gold_standard_complex_appearance()


        ############################################################
        # DETERMINE WHICH APPEARING GSC ARE FOUND BY ALGORITHM
        ############################################################
        def found_gsc():
            retval = []
            for A in self.gold_standard_complexes_appearing_in_dataset:
                for B in self.densityThreshold_sizeThreshold_merged_cluster_list:
                    if jaccard_similarity(A, B) >= self.found_gsc_jaccard_threshold:
                        retval.append(A)
            self.found = retval
            return retval

        found_gsc()



        ############################################################
        # DETERMINE WHICH GSC APPEARING IN DATASET ARE NOT FOUND BY ALGORITHM
        ############################################################
        def not_found_gsc():
            retval = []
            for A in self.gold_standard_complexes_appearing_in_dataset:
                found = False
                for B in self.densityThreshold_sizeThreshold_merged_cluster_list:
                    if jaccard_similarity(A, B) >= self.found_gsc_jaccard_threshold:
                        found = True
                if not found:
                    retval.append(A)
            self.not_found = retval
            return retval

        not_found_gsc()

        ############################################################
        # FUNCTION TO CALCULATE STATS ABOUT A LIST OF CLUSTERS
        ############################################################
        def calculate_clusters_stats(cluster_list, output_map):
            output_map['clusters'] = dict()
            total_density = 0
            total_cohesiveness = 0
            total_size = 0
            for cluster in cluster_list:
                key = tuple([self.graph.id_to_name[id] for id in cluster])
                output_map['clusters'][key] = dict()
                _cohesiveness = cohesiveness(self, list(cluster))
                _size = len(cluster)
                _density = density(self, list(cluster))
                output_map['clusters'][key]['cohesiveness'] = _cohesiveness
                output_map['clusters'][key]['size'] = _size
                output_map['clusters'][key]['density'] = _density
                total_cohesiveness+=_cohesiveness
                total_density += _density
                total_size +=_size
            if len(cluster_list)==0:
                print("the current cluster list has length 0")
                output_map['average_cohesiveness'] = 0
                output_map['average_density'] = 0
                output_map['average_size'] = 0
                return
            output_map['average_cohesiveness'] = total_cohesiveness /float(len(cluster_list))
            output_map['average_density'] = total_density /float(len(cluster_list))
            output_map['average_size'] = total_size /float(len(cluster_list))


        ############################################################
        # CALCULATE STATS ABOUT THE FINAL RESULT
        ############################################################
        calculate_clusters_stats(
            self.densityThreshold_sizeThreshold_merged_cluster_list,
            self.final_clusters_stats)

        ############################################################
        # CALCULATE STATS ABOUT THE GSC APPEARING IN DATASET
        ############################################################
        calculate_clusters_stats(
            self.gold_standard_complexes_appearing_in_dataset,
            self.gsc_appearing_stats)

        ############################################################
        # CALCULATE STATS ABOUT THE GSC APPEARING IN DATASET FOUND BY ALGORITHM
        ############################################################
        print(len(self.found))
        calculate_clusters_stats(
            self.found,
            self.gsc_appearing_found_stats)

        ############################################################
        # CALCULATE STATS ABOUT THE GSC APPEARING IN DATASET NOT FOUND BY ALGORITHM
        ############################################################
        calculate_clusters_stats(
            self.not_found,
            self.gsc_appearing_notFound_stats)

        ############################################################
        # DISPLAY DETAILS OF THE GSC THAT APPEARED IN THE DATASET BASED ON WHETHER OR NOT THE ALGORITHM FOUND THEM
        ############################################################
        print("################ FOUND AND UNFOUND GSC DETAILS #######################")

        def found_and_unfound_details():
            print("FOUND")
            cohesiveness_tot = 0
            density_tot = 0
            length_tot = 0
            for c in self.found:
                c_cohesiveness = cohesiveness(self, c)
                cohesiveness_tot += c_cohesiveness
                c_density = density(self, c)
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
                c_cohesiveness = cohesiveness(self, c)
                cohesiveness_tot += c_cohesiveness
                c_density = density(self, c)
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

        if len(self.found):
            found_and_unfound_details()

        ##########################################################################
        #  DETERMINE THE QUALITY OF THE RESULT USING ORIGINAL AUTHORS' MEASURES  #
        ##########################################################################
        print("################ QUALITY #######################")

        def get_quality():
            import subprocess
            res = subprocess.check_output(["python2",
                                           "../cl1_reproducibility/reproducibility/scripts/match_standalone.py",
                                           self.gold_standard_filename,
                                           "../complexes/" + self.output_filename])
            for line in res.splitlines():
                print(line)
                a = str(line)
                a = a.replace("b", "").replace("=", "").replace("\'", "").split()
                self.quality_report[a[0]] = float(a[1])

            print(self.run_title)

        get_quality()

        ############################################################
        # LOG RUN INFO
        ############################################################
        f = open("../../run_log","a+")
        f.write(str(self.run_title)+"\n")
        f.write(str(self.argument_dict)+"\n")
        f.write(str(self.quality_report)+"\n\n")
        f.close()

        ############################################################
        # STORE THE CURRENT OBJECT USING PICKLE
        ############################################################
        def store_self():
            f_name = "../pickles/pickle+" + self.run_title
            f = open(f_name, 'ab')
            pickle.dump(self, f)
            f.close()

        store_self()

        f_name = "../pickles/most_recent"
        f = open(f_name, 'wb')
        title = {'title': "pickles/pickle+"+self.run_title}
        pickle.dump(title, f)







