from src.CL1R.cl1r import CL1_Randomized
from src.COMMON.cmn import *

if __name__=="__main__":
    pr = cProfile.Profile()
    pr.enable()
    # ... do something ...
    a = CL1_Randomized("../cl1_datasets/datasets",
                       "gavin2006_socioaffinities_rescaled.txt",
                       'Dummy_quality',
                       density_threshold=.3,
                       merge_threshold=.9,
                       penalty_value_per_node=2,
                       randomized_construction_bool=True,
                       rng_seed=None,
                       number_of_shakes=0,
                       number_of_bad_adds=2,
                       sort_seeds_by="weight",
                       care_about_cuts=False,
                       seed_from_all=False,
                       gsc_appearance_ratio_threshold=.9,
                       found_gsc_jaccard_threshold=.8,
                       gold_standard_filename="../cl1_gold_standard/gold_standard/mips_3_100.txt") 

    pr.disable()
    s = io.StringIO()
    sortby = SortKey.CUMULATIVE
    ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    ps.print_stats()
    print(s.getvalue())






    # pr = cProfile.Profile()
    # pr.enable()
    # # ... do something ...
    # b = CL1_Randomized("../cl1_datasets/datasets",
    #                    "gavin2006_socioaffinities_rescaled.txt",
    #                    'Dummy_quality',
    #                    density_threshold=.3,
    #                    merge_threshold=.9,
    #                    penalty_value_per_node=2,
    #                    randomized_construction_bool=False,
    #                    rng_seed=None,
    #                    number_of_shakes=0,
    #                    number_of_bad_adds=2,
    #                    sort_seeds_by="weight",
    #                    care_about_cuts=False,
    #                    seed_from_all=True,
    #                    gsc_appearance_ratio_threshold=.9,
    #                    found_gsc_jaccard_threshold=.8,
    #                    gold_standard_filename="../cl1_gold_standard/gold_standard/mips_3_100.txt")
    #
    # pr.disable()
    # s = io.StringIO()
    # sortby = SortKey.CUMULATIVE
    # ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    # ps.print_stats()
    # print(s.getvalue())

# To beat: my implementation of the original
# ORIGINAL
# 189 reference complexes, 254 predicted complexes
# b'acc = 0.3698'
# b'cws = 0.3388'
# b'frac = 0.4074'
# b'mmr = 0.2116'
# b'ppv = 0.4035'
# b'sep = 0.2429'
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# a = CL1_Randomized("cl1_datasets/datasets", "gavin2006_socioaffinities_rescaled.txt", 'Dummy_quality',
#                    density_threshold=.3,
#                    merge_threshold=.9,
#                    penalty_value_per_node=2,
#                    randomized_construction_bool=True,
#                    rng_seed=None,
#                    number_of_shakes=1,
#                    number_of_bad_adds=1,
#                    sort_seeds_by="weight",
#                    care_about_cuts=True,
#                    seed_from_all=True,
#                    gsc_appearance_ratio_threshold=.9,
#                    found_gsc_jaccard_threshold=.8,
#                    gold_standard_filename="cl1_gold_standard/gold_standard/mips_3_100.txt")
# above, without reverting to best seen
# 189
# reference
# complexes, 301
# predicted
# complexes
# b'acc = 0.3741'
# b'cws = 0.3511'
# b'frac = 0.4074'
# b'mmr = 0.2135'
# b'ppv = 0.3985'
# b'sep = 0.2182'
# 189
# reference
# complexes, 300
# predicted
# complexes
# b'acc = 0.3759'
# b'cws = 0.3545'
# b'frac = 0.4180'
# b'mmr = 0.2175'
# b'ppv = 0.3985'
# b'sep = 0.2206'
# 189
# reference
# complexes, 301
# predicted
# complexes
# b'acc = 0.3756'
# b'cws = 0.3545'
# b'frac = 0.4180'
# b'mmr = 0.2173'
# b'ppv = 0.3979'
# b'sep = 0.2211'
# 189
# reference
# complexes, 298
# predicted
# complexes
# b'acc = 0.3750'
# b'cws = 0.3537'
# b'frac = 0.4127'
# b'mmr = 0.2147'
# b'ppv = 0.3977'
# b'sep = 0.2222'
# gavin2006_socioaffinities_rescaled + 2019 - 0
# 9 - 0
# 9_11: 43:34: 255366
# 189
# reference
# complexes, 300
# predicted
# complexes
# b'acc = 0.3760'
# b'cws = 0.3541'
# b'frac = 0.4180'
# b'mmr = 0.2176'
# b'ppv = 0.3992'
# b'sep = 0.2232'
# gavin2006_socioaffinities_rescaled + 2019 - 0
# 9 - 0
# 9_11: 45:03: 468226


# ABOVE, with reverting to best seen
# 189
# reference
# complexes, 378
# predicted
# complexes
# b'acc = 0.3780'
# b'cws = 0.3520'
# b'frac = 0.4127'
# b'mmr = 0.2378'
# b'ppv = 0.4059'
# b'sep = 0.2016'
# 189
# reference
# complexes, 372
# predicted
# complexes
# b'acc = 0.3772'
# b'cws = 0.3499'
# b'frac = 0.4127'
# b'mmr = 0.2377'
# b'ppv = 0.4068'
# b'sep = 0.2034'
# 189
# reference
# complexes, 377
# predicted
# complexes
# b'acc = 0.3786'
# b'cws = 0.3520'
# b'frac = 0.4127'
# b'mmr = 0.2380'
# b'ppv = 0.4073'
# b'sep = 0.2021'












# a = CL1_Randomized("cl1_datasets/datasets", "gavin2006_socioaffinities_rescaled.txt", 'Dummy_quality',
#                    density_threshold=.3,
#                    merge_threshold=.8,
#                    penalty_value_per_node=2,
#                    randomized_construction_bool=True,
#
#                    rng_seed=None,
#                    number_of_shakes=0,
#                    number_of_bad_adds=1,
#                    sort_seeds_by="weight",
#                    care_about_cuts=True,
#                    seed_from_all=False,
#                    gsc_appearance_ratio_threshold=.9,
#                    found_gsc_jaccard_threshold=.8,
#                    gold_standard_filename="cl1_gold_standard/gold_standard/mips_3_100.txt")
# exec("networkx_tests.py")







# a = CL1_Randomized("cl1_datasets/datasets", "krogan2006_extended.txt", 'Dummy_quality',
#                    density_threshold=.3,
#                    merge_threshold=.9,
#                    penalty_value_per_node=2,
#                    randomized_construction_bool=True,
#                    rng_seed=None,
#                    number_of_shakes=1,
#                    number_of_bad_adds=1,
#                    sort_seeds_by="weight",
#                    care_about_cuts=True,
#                    seed_from_all=True,
#                    gsc_appearance_ratio_threshold=.9,
#                    found_gsc_jaccard_threshold=.8,
#                    gold_standard_filename="cl1_gold_standard/gold_standard/mips_3_100.txt")
"""
################ QUALITY #######################
189 reference complexes, 1591 predicted complexes
b'acc = 0.3572'
b'cws = 0.3189'
b'frac = 0.4762'
b'mmr = 0.2599'
b'ppv = 0.4001'
b'sep = 0.1251'
krogan2006_extended+2019-09-19_09:27:51:114908

"""
# a = CL1_Randomized("cl1_datasets/datasets", "krogan2006_extended.txt", 'Dummy_quality',
#                    density_threshold=.3,
#                    merge_threshold=.9,
#                    penalty_value_per_node=2,
#                    randomized_construction_bool=False,
#                    rng_seed=None,
#                    number_of_shakes=0,
#                    number_of_bad_adds=1,
#                    sort_seeds_by="weight",
#                    care_about_cuts=True,
#                    seed_from_all=True,
#                    gsc_appearance_ratio_threshold=.9,
#                    found_gsc_jaccard_threshold=.8,
#                    gold_standard_filename="cl1_gold_standard/gold_standard/mips_3_100.txt")
"""
################ QUALITY #######################
189 reference complexes, 1131 predicted complexes
b'acc = 0.3709'
b'cws = 0.3571'
b'frac = 0.4444'
b'mmr = 0.2570'
b'ppv = 0.3853'
b'sep = 0.1387'
krogan2006_extended+2019-09-19_09:26:25:901378
"""


# a = CL1_Randomized("cl1_datasets/datasets", "biogrid_yeast_physical_unweighted+naively_weighted.txt", 'Dummy_quality',
#                    density_threshold=.3,
#                    merge_threshold=.9,
#                    penalty_value_per_node=2,
#                    randomized_construction_bool=False,
#                    rng_seed=None,
#                    number_of_shakes=1,
#                    number_of_bad_adds=1,
#                    sort_seeds_by="weight",
#                    care_about_cuts=True,
#                    seed_from_all=True,
#                    gsc_appearance_ratio_threshold=.9,
#                    found_gsc_jaccard_threshold=.8,
#                    gold_standard_filename="cl1_gold_standard/gold_standard/mips_3_100.txt")



# a = CL1_Randomized("cl1_datasets/datasets", "krogan2006_extended.txt", 'Dummy_quality',
#                    density_threshold=.3,
#                    merge_threshold=.9,
#                    penalty_value_per_node=2,
#                    randomized_construction_bool=True,
#                    rng_seed=None,
#                    number_of_shakes=1,
#                    number_of_bad_adds=1,
#                    sort_seeds_by="weight",
#                    care_about_cuts=True,
#                    seed_from_all=True,
#                    gsc_appearance_ratio_threshold=.9,
#                    found_gsc_jaccard_threshold=.8,
#                    gold_standard_filename="cl1_gold_standard/gold_standard/mips_3_100.txt")
#