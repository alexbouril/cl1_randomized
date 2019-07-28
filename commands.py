#exec(open("commands.py").read())
#python commands.py
from cl1_randomized import *
a = CL1_Randomized("cl1_datasets/datasets", "gavin2006_socioaffinities_rescaled.txt", 'Dummy_quality',
                   "tarea_regulardensitythreshold18.txt",
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


#
# average cohesiveness of found=  0.7287811118428832
# average density of found=  1.1477001097007222
# average cohesiveness of NOT found=  0.306053788274819
# average density of NOT found=  0.7658721978209649

# average cohesiveness of found=  0.7032572263002919
# average density of found=  1.1277445595704951
# average cohesiveness of NOT found=  0.2990207643943256
# average density of NOT found=  0.7580406234794939









# a = CL1_Randomized("cl1_datasets/datasets", "gavin2006_socioaffinities_rescaled.txt", 'Dummy_quality',
#                    "tarea_regulardensitythreshold18.txt",
#                    density_threshold=.3,
#                    penalty_value_per_node=2,
#                    randomized_construction_bool=True,
#                    sort_seeds_by="degree",
#                    seed_from_all=False,
#                    gold_standard_filename="cl1_gold_standard/gold_standard/mips_3_100.txt")
#
# found:  20
# not found:  45
# average cohesiveness of found=  0.7160846952573257
# average density of found=  1.168026004215686
# average length of found=  5.7
# average cohesiveness of NOT found=  0.3023026995668875
# average density of NOT found=  0.7483534022169867
# average length of NOT found=  7.133333333333334
# #######################################
# 189 reference complexes, 251 predicted complexes
# b'acc = 0.3674'
# b'cws = 0.3355'
# b'frac = 0.4021'
# b'mmr = 0.2112'
# b'ppv = 0.4024'
# b'sep = 0.2407'
# --------------------------------------
# found:  19
# not found:  46
# average cohesiveness of found=  0.7287811118428832
# average density of found=  1.1477001097007222
# average length of found=  5.7894736842105265
# average cohesiveness of NOT found=  0.306053788274819
# average density of NOT found=  0.7658721978209649
# average length of NOT found=  7.065217391304348
# #######################################
# 189 reference complexes, 248 predicted complexes
# b'acc = 0.3655'
# b'cws = 0.3321'
# b'frac = 0.3915'
# b'mmr = 0.2062'
# b'ppv = 0.4024'
# b'sep = 0.2403'
# --------------------------------------
# found:  19
# not found:  46
# average cohesiveness of found=  0.7287811118428832
# average density of found=  1.1477001097007222
# average length of found=  5.7894736842105265
# average cohesiveness of NOT found=  0.306053788274819
# average density of NOT found=  0.7658721978209649
# average length of NOT found=  7.065217391304348
# #######################################
# 189 reference complexes, 254 predicted complexes
# b'acc = 0.3698'
# b'cws = 0.3388'
# b'frac = 0.4074'
# b'mmr = 0.2116'
# b'ppv = 0.4035'
# b'sep = 0.2429'









# a = CL1_Randomized("cl1_datasets/datasets", "gavin2006_socioaffinities_rescaled.txt", 'Dummy_quality',
#                    "tarea_regulardensitythreshold18.txt",
#                    density_threshold=.3,
#                    penalty_value_per_node=2,
#                    randomized_construction_bool=True,
#                    sort_seeds_by="degree",
#                    seed_from_all=True,
#                    gold_standard_filename="cl1_gold_standard/gold_standard/mips_3_100.txt")
#
# --------------------------------------
# found:  23
# not found:  42
# average cohesiveness of found=  0.6826099368166605
# average density of found=  1.1162283312020458
# average length of found=  5.6521739130434785
# average cohesiveness of NOT found=  0.2910782580684111
# average density of NOT found=  0.7467421801531208
# average length of NOT found=  7.261904761904762
# #######################################
# 189 reference complexes, 295 predicted complexes
# b'acc = 0.3757'
# b'cws = 0.3401'
# b'frac = 0.4074'
# b'mmr = 0.2156'
# b'ppv = 0.4151'
# b'sep = 0.2354'

# --------------------------------------
# found:  23
# not found:  42
# average cohesiveness of found=  0.6826099368166605
# average density of found=  1.1162283312020458
# average length of found=  5.6521739130434785
# average cohesiveness of NOT found=  0.2910782580684111
# average density of NOT found=  0.7467421801531208
# average length of NOT found=  7.261904761904762
# #######################################
# 189 reference complexes, 298 predicted complexes
# b'acc = 0.3778'
# b'cws = 0.3431'
# b'frac = 0.4074'
# b'mmr = 0.2158'
# b'ppv = 0.4159'
# b'sep = 0.2370'
# --------------------------------------
# found:  23
# not found:  42
# average cohesiveness of found=  0.6826099368166605
# average density of found=  1.1162283312020458
# average length of found=  5.6521739130434785
# average cohesiveness of NOT found=  0.2910782580684111
# average density of NOT found=  0.7467421801531208
# average length of NOT found=  7.261904761904762
# #######################################
# 189 reference complexes, 296 predicted complexes
# b'acc = 0.3776'
# b'cws = 0.3431'
# b'frac = 0.4074'
# b'mmr = 0.2158'
# b'ppv = 0.4156'
# b'sep = 0.2377'




# from cl1_randomized import *
# a = CL1_Randomized("cl1_datasets/datasets", "gavin2006_socioaffinities_rescaled.txt", 'Dummy_quality',
#                    "tarea_regulardensitythreshold18.txt",
#                    density_threshold=.3,
#                    penalty_value_per_node=2,
#                    randomized_construction_bool=True,
#                    sort_seeds_by="degree",
#                    care_about_cuts=False,
#                    seed_from_all=True,
#                    gold_standard_filename="cl1_gold_standard/gold_standard/mips_3_100.txt")
# --------------------------------------
# found:  23
# not found:  42
# average cohesiveness of found=  0.6826099368166605
# average density of found=  1.1162283312020458
# average length of found=  5.6521739130434785
# average cohesiveness of NOT found=  0.2910782580684111
# average density of NOT found=  0.7467421801531208
# average length of NOT found=  7.261904761904762
# #######################################
# 189 reference complexes, 299 predicted complexes
# b'acc = 0.3762'
# b'cws = 0.3401'
# b'frac = 0.4074'
# b'mmr = 0.2158'
# b'ppv = 0.4162'
# b'sep = 0.2342'
# --------------------------------------
# found:  23
# not found:  42
# average cohesiveness of found=  0.6826099368166605
# average density of found=  1.1162283312020458
# average length of found=  5.6521739130434785
# average cohesiveness of NOT found=  0.2910782580684111
# average density of NOT found=  0.7467421801531208
# average length of NOT found=  7.261904761904762
# #######################################
# 189 reference complexes, 300 predicted complexes
# b'acc = 0.3782'
# b'cws = 0.3439'
# b'frac = 0.4127'
# b'mmr = 0.2175'
# b'ppv = 0.4159'
# b'sep = 0.2350'
# --------------------------------------
# found:  23
# not found:  42
# average cohesiveness of found=  0.6826099368166605
# average density of found=  1.1162283312020458
# average length of found=  5.6521739130434785
# average cohesiveness of NOT found=  0.2910782580684111
# average density of NOT found=  0.7467421801531208
# average length of NOT found=  7.261904761904762
# #######################################
# 189 reference complexes, 296 predicted complexes
# b'acc = 0.3757'
# b'cws = 0.3401'
# b'frac = 0.4074'
# b'mmr = 0.2153'
# b'ppv = 0.4151'
# b'sep = 0.2352'