from src.CL1R.cl1r import *
def convert_to_new_names(reference, graph):
     return [[graph.name_to_id.get(protein,-1) for protein in line.split()] for line in reference.readlines()]
def get_edges(complex, graph):
    list_of_edges=[]
    for i in range(len(complex)):
        for j in range(len(complex)):
            try:
                if graph.hash_graph[i][j]:
                    list_of_edges.append((i,j, graph.hash_graph[i][j]))
            except KeyError:
                pass
    return list_of_edges

datasets = []
datasets.append("gavin2006_socioaffinities_rescaled.txt")
datasets.append("krogan2006_extended.txt")
datasets.append("collins2007.txt")
datasets.append("biogrid_yeast_physical_unweighted+naively_weighted.txt")

references = []
references.append("mips_3_100.txt")
references.append("sgd.txt")

dataset_base_file_path = "../../cl1_datasets/datasets/"
reference_base_file_path= "../../cl1_gold_standard/gold_standard/"
i=0
j=0
original_graph_filename=datasets[i]
reference_filename=references[j]
original_graph_filename = original_graph_filename
cl1 = CL1_Randomized(dataset_base_file_path,
                           datasets[i],
                           'Dummy_quality',
                           density_threshold=.15,
                           merge_threshold=.9,
                           penalty_value_per_node=2,
                           use_original_penalty=False,
                           randomized_construction_bool=True,
                           rng_seed=None,
                           use_mixed_measure_find=True,
                           number_consecutive_mixed_measure_finds=1,
                           use_suboptimal_adds = True,
                           number_consecutive_suboptimal_adds=4,
                           number_of_shakes=0,
                           number_of_bad_adds=2,
                           sort_seeds_by="weight",
                           care_about_cuts=False,
                           seed_from_all=True,
                           gsc_appearance_ratio_threshold=.9,
                           found_gsc_jaccard_threshold=.8,
                           gold_standard_filename=reference_filename,
                           save_self=False,
                            process=False)
reference = open(reference_base_file_path+reference_filename,"r")
complex_list = convert_to_new_names(reference,cl1.graph)
edge_list = get_edges(complex_list[0], cl1.graph)

protein_involvement = [(sum([protein in cl1.graph.hash_graph for protein in complex])/ len(complex)) for complex in complex_list]
cohesiveness_list=[cohesiveness(cl1, complex, only_allow_good_queries=False) for complex in complex_list ]
density_list = [density(cl1, complex, only_allow_good_queries=False) for complex in complex_list]

for i in range(len(complex_list)):
    if protein_involvement[i]>.70:
        print("----------- {} -----------".format(i))
        print(complex_list[i])
        print("protein_involvement[i]", protein_involvement[i])
        print("density_list[i]", density_list[i])
        print("cohesiveness_list[i]", cohesiveness_list[i])

