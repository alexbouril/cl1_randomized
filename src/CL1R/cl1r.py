from src.COMMON.cmn import *
from src.GRAPH.graph import *
from src.CONSTRUCTION.randomized_construction import randomized_construction
from src.CONSTRUCTION.original_construction import original_construction
from src.QUALITY.quality import density, cohesiveness
from src.MERGING.merge import merge
from src.MERGING.zipper_merge import zipper_merge
from src.THRESHOLDING.thresholding_functions import sizeThreshold, densityThreshold


# TODO: package/syspath,1-1 comparison, Cython, Threading, Command-line args, currying, class variable, hidden class variable, numpy, inheritance, tracemalloc

# TODO: check the effect of using weights for density, not just number edges
# TODO: check about setting upper bound for lengths!  should we allow complexes with length greater than 300?


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
                 use_original_penalty=True,
                 randomized_construction_bool= False,
                 rng_seed = None,
                 use_mixed_measure_find=False,
                 number_consecutive_mixed_measure_finds=1,
                 use_suboptimal_adds=False,
                 number_consecutive_suboptimal_adds=1,
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
                 gold_standard_filename = "",
                 save_self=False,
                 process=True):
        self.time_of_run = str(datetime.datetime.now()).replace(" ","_").replace(".",":")
        self.run_title = original_graph_filename.replace(".txt", "")+"+"+self.time_of_run
        self.argument_dict = locals()
        ############################################################
        # STORE ARGUMENTS
        ############################################################
        self.base_file_path = base_file_path
        self.original_graph_filename = original_graph_filename
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
        if use_original_penalty:
            self.penalty_value_per_node = penalty_value_per_node
        else:
            self.penalty_value_per_node = .5 * sum([len(self.graph.hash_graph[v]) for v in self.graph.hash_graph])/len(self.graph.hash_graph)
            self.penalty_value_per_node = sum([sum([self.graph.hash_graph[v][u] for u in self.graph.hash_graph[v]]) for v in self.graph.hash_graph])/len(self.graph.hash_graph)
        self.randomized_construction_bool = randomized_construction_bool
        self.rng_seed = rng_seed
        np.random.seed(rng_seed)
        self.use_mixed_measure_find = use_mixed_measure_find
        self.number_consecutive_mixed_measure_finds = number_consecutive_mixed_measure_finds
        self.use_suboptimal_adds=use_suboptimal_adds
        self.number_consecutive_suboptimal_adds=number_consecutive_suboptimal_adds
        self.rng_initial_state = np.random.get_state()
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
        self.save_self = save_self
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
        if process:
            self.process()

    def process(self):
        ############################################################
        # INITIAL CONSTRUCTION
        ############################################################
        # if self.randomized_construction_bool:
        #     randomized_construction(self)
        # # TODO reimplement the oc function using new classes and package framework
        # else:
        #     original_construction(self)
        randomized_construction(self)
        # self.seed_from_all=False
        # randomized_construction(self)
        # original_construction(self)
        # bsdfa = len(self.initial_clustering)
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
        print("len(self.merged_cluster_list)",len(self.merged_cluster_list))
        ############################################################
        # ZIPPER MERGE COMPLEXES
        ############################################################
        # print("____________________________zipper_merge_________________________________")
        # time.sleep(1)
        # zipper_merged_clusters = zipper_merge(self, self.initial_cluster_list)
        # # add the zipper merged clusters to the merged cluster list
        # self.merged_cluster_list+=zipper_merged_clusters
        # print("len(zipper_merged_clusters)", len(zipper_merged_clusters))
        # time.sleep(1)
        ############################################################
        # REMERGE
        ############################################################
        # self.merged_cluster_list= merge(merge_threshold=self.merge_threshold,
        #                                  source=self.merged_cluster_list)
        # print("len(self.merged_cluster_list)",len(self.merged_cluster_list))
        ############################################################
        # THRESHOLD BASED ON COHESIVENESS???
        ############################################################
        # todo
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
        # self.densityThreshold_sizeThreshold_merged_cluster_list  = self.sizeThreshold_merged_cluster_list
        ############################################################
        # count the number of duplicates
        ############################################################
        result = self.densityThreshold_sizeThreshold_merged_cluster_list
        to_delete = set()
        for idx, cluster in enumerate(result):
            for idx2, cluster2 in enumerate(result):
                if idx == idx2 or idx2 in to_delete:
                    continue
                if cluster == cluster2:
                    to_delete.add(idx2)
        print("number of duplicates: ", len(to_delete))
        time.sleep(1)
        ############################################################
        # WRITE THE FINAL RESULT INTO A TEXT FILE
        ############################################################

        def write_final_clusters(collection, filename):
            f = open(filename, "w+")
            counter = 1
            for cluster in collection:
                s = ""
                for id in cluster:
                    s += self.graph.id_to_name[id] + "\t"
                print("Cluster #%s: " % str(counter), len(cluster), " proteins")
                print(s)
                f.write(s + "\n")
                counter += 1
            f.close()
        write_final_clusters(self.densityThreshold_sizeThreshold_merged_cluster_list,
                             "../complexes/" + self.output_filename)
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

        def get_quality(reference_file, complexes_file):
            import subprocess
            res = subprocess.check_output(["python2",
                                           "../cl1_reproducibility/reproducibility/scripts/match_standalone.py",
                                           reference_file,
                                           complexes_file])
            composite_score = 0
            acc = -1
            cws = -1
            frac= -1
            mmr = -1
            for line in res.splitlines():
                print(line)
                a = str(line)
                a = a.replace("b", "").replace("=", "").replace("\'", "").split()
                self.quality_report[a[0]] = float(a[1])
                if a[0] in ['acc', "frac","mmr"]:# "cws"]:
                    val =float(a[1])
                    composite_score+=val
                    attribute = a[0]
                    if attribute=='acc':
                        acc = val
                    elif attribute=="cws":
                        cws= val
                    elif attribute=='frac':
                        frac=val
                    elif attribute =='mmr':
                        mmr=val
            print(self.run_title)
            return composite_score, acc, cws, frac, mmr

        scores = {}
        dataset_to_cyto_default_result = {"gavin2006_socioaffinities_rescaled.txt": "cyto_gavin_paper_defaults.txt",
                                          "collins2007.txt":"cyto_collins_paper_defaults.txt",
                                          "krogan2006_extended.txt":"cyto_krogan_paper_defaults.txt",
                                          "biogrid_yeast_physical_unweighted+naively_weighted.txt":"cyto_biog+nW_cl1_fromUnused_overlapPoint8.txt"
        }
        if self.original_graph_filename in dataset_to_cyto_default_result:
            scores["mips"]={"mine":{}, "theirs":{}}
            scores["sgd"]={"mine":{}, "theirs":{}}
            print("using mips_3_100.txt")
            print("Mine")
            scores["mips"]['mine']["total"],\
            scores["mips"]['mine']["acc"],\
            scores["mips"]['mine']["cws"], \
            scores["mips"]['mine']["frac"],\
            scores["mips"]['mine']["mmr"]\
                = get_quality(self.gold_standard_filename, "../complexes/" + self.output_filename)
            print("Theirs")
            scores["mips"]['theirs']["total"],\
            scores["mips"]['theirs']["acc"],\
            scores["mips"]['theirs']["cws"],\
            scores["mips"]['theirs']["frac"],\
            scores["mips"]['theirs']["mmr"]\
                =get_quality(self.gold_standard_filename, "../cl1_datasets/datasets/%s"%dataset_to_cyto_default_result[self.original_graph_filename])
            self.gold_standard_filename = "../cl1_gold_standard/gold_standard/sgd.txt"
            print("using sgd.txt")
            print("Mine")
            scores["sgd"]['mine']["total"], \
            scores["sgd"]['mine']["acc"], \
            scores["sgd"]['mine']["cws"], \
            scores["sgd"]['mine']["frac"], \
            scores["sgd"]['mine']["mmr"] \
                =get_quality(self.gold_standard_filename, "../complexes/" + self.output_filename)
            print("Theirs")
            scores["sgd"]['theirs']["total"], \
            scores["sgd"]['theirs']["acc"], \
            scores["sgd"]['theirs']["cws"], \
            scores["sgd"]['theirs']["frac"], \
            scores["sgd"]['theirs']["mmr"] \
                =get_quality(self.gold_standard_filename, "../cl1_datasets/datasets/%s"%dataset_to_cyto_default_result[self.original_graph_filename])


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
        if self.save_self:
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
        print(scores)
        self.scores = scores







