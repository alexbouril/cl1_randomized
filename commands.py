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
                   sort_seeds_by="weight",
                   care_about_cuts=True,
                   seed_from_all=False,
                   gsc_appearance_ratio_threshold=.8,
                   found_gsc_jaccard_threshold=.8,
                   gold_standard_filename="cl1_gold_standard/gold_standard/mips_3_100.txt")
a.process()










