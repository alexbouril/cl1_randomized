#exec(open("commands.py").read())
#exec(open("commands.py").read())
#python commands.py
from cl1_randomized import *
a = CL1_Randomized("cl1_datasets/datasets", "gavin2006_socioaffinities_rescaled.txt", 'Dummy_quality',
                   density_threshold=.3,
                   merge_threshold=.8,
                   penalty_value_per_node=2,
                   randomized_construction_bool=True,
                   number_of_shakes=1,
                   number_of_bad_adds=2,
                   sort_seeds_by="degree",
                   care_about_cuts=True,
                   seed_from_all=False,
                   gsc_appearance_ratio_threshold=.8,
                   found_gsc_jaccard_threshold=.8,
                   gold_standard_filename="cl1_gold_standard/gold_standard/mips_3_100.txt")
a.process()



print(a.run_title)
print(a.quality_report)
print(a.initial_clustering)
print(a.final_clusters_stats)
print(a.found)
print(a.gsc_appearing_found_stats)
print(a.gsc_appearing_notFound_stats)
print(a.gsc_appearing_stats)
for d in [a.final_clusters_stats, a.gsc_appearing_stats, a.gsc_appearing_found_stats, a.gsc_appearing_notFound_stats]:
    print(d['average_cohesiveness'], d['average_density'], d['average_size'])





