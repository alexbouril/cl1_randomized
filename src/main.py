from src.CL1R.cl1r import CL1_Randomized
from src.COMMON.cmn import *
import pprint as pp
import collections
import numpy as np
import matplotlib.pyplot as plt
import re
import time
import os
import pprint
import sys
run_id = time.strftime("run_%Y_%m_%d-%H_%M_%S")
dir_path = "../runs/"+run_id
os.mkdir(path=dir_path, mode=0o755)
notes_path = dir_path+"/notes"
notes_file = open(notes_path, 'a')
# sorting by weight vs sorting by degree
# penalty value fi
if __name__=="__main__":
    pr = cProfile.Profile()
    pr.enable()
    scores = dict()
    datasets=[]
    datasets.append("gavin2006_socioaffinities_rescaled.txt")
    datasets.append("krogan2006_extended.txt")
    datasets.append("collins2007.txt")
    datasets.append("biogrid_yeast_physical_unweighted+naively_weighted.txt")
    for i, dataset in enumerate(datasets):
        result = CL1_Randomized("../cl1_datasets/datasets",
                           dataset,
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
                           gold_standard_filename="../cl1_gold_standard/gold_standard/mips_3_100.txt",
                           save_self=False)
        if i==0:
            notes_file.write(pprint.pformat(result.argument_dict, indent=4)+"\n")
        notes_file.write(result.run_title+"\n\n")
        scores[dataset]=result.scores
    pp.pprint(scores)
    improvements = collections.defaultdict(dict)
    for dataset in scores:
        for ref in scores[dataset]:
            improvements[dataset][ref] = -1 + (scores[dataset][ref]['mine']['total'] / scores[dataset][ref]['theirs']['total'])
    notes_file.write(str(datasets)+"\n\n")
    notes_file.write(pprint.pformat(scores, indent=4)+"\n\n")
    notes_file.write(pprint.pformat(improvements, indent=4)+"\n\n")
    pr.disable()
    s = io.StringIO()
    sortby = SortKey.CUMULATIVE
    ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    ps.print_stats()
    print(s.getvalue())
    pp.pprint(scores)
    pp.pprint(improvements)
    #
    # import matplotlib.pyplot as plt
    # for gold_standard in ['mips', 'sgd']:
    #     N = len(datasets)
    #     fracs = tuple([scores[d][gold_standard]["mine"]['frac'] for d in datasets])
    #     accs = tuple([scores[d][gold_standard]["mine"]['acc'] for d in datasets])
    #     mmrs = tuple([scores[d][gold_standard]["mine"]['mmr'] for d in datasets])
    #     ind = np.arange(N)  # the x locations for the groups
    #     width = 0.35  # the width of the bars: can also be len(x) sequence
    #     p1 = plt.bar(ind, fracs, width)
    #     p2 = plt.bar(ind, accs, width,
    #                  bottom=fracs)
    #     p3 = plt.bar(ind, mmrs, width,
    #                  bottom=fracs)
    #     plt.ylabel('Scores')
    #     plt.title('Scores')
    #     plt.xticks(ind, (d for d in datasets))
    #     plt.yticks(np.arange(0, 2, .1))
    #     plt.legend((p1[0], p2[0], p3[0]), ('Fraction Matched', 'Accuracy', 'Maximum Matching Ratio'))
    #     plt.show()


#########################################
    datasets = tuple(datasets)
    segments = 3

    # multi-dimensional data
    for gold_standard in ['mips', 'sgd']:
        elements = 2*len(datasets)
        fracs = np.empty((elements,))
        accs = np.empty((elements,))
        mmrs = np.empty((elements,))
        for measure, result in zip(['frac', 'acc', "mmr"],[fracs, accs, mmrs]):
            mine = [scores[d][gold_standard]["mine"][measure] for d in datasets]
            theirs = [scores[d][gold_standard]["theirs"][measure] for d in datasets]
            result[0::2] = theirs
            result[1::2] = mine

        data =np.array( [fracs,
                         accs,
                         mmrs
                         ])
        section_labels = data.transpose()
        print(data)
        print(section_labels)
        y_pos = np.arange(len(datasets)*2).astype(float)

        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111)
        colors = ['#7bdb7c','#f0c37a','#6face8']
        patch_handles = []
        # left alignment of data starts at zero
        left = np.zeros(len(datasets*2))
        print(y_pos)
        for i in range(len(y_pos)):
            if i%2==1:
                y_pos[i]-=.1
        for i, d in enumerate(data):
            print('y_pos: ', y_pos, d)
            to_append = ax.barh(y_pos, d, color=colors[i % len(colors)],align='center', left=left)
            patch_handles.append(to_append)
            left += d
        print("patch_handles",patch_handles)
        # search all of the bar segments and annotate
        for j in range(len(patch_handles)):
            for i, patch in enumerate(patch_handles[j].get_children()):
                bl = patch.get_xy()
                x = 0.5 * patch.get_width() + bl[0]
                y = 0.5 * patch.get_height() + bl[1]
                print(patch.get_height(), patch.get_width())
                ax.text(x, y, "%s" % (str(section_labels[i, j])), ha='center')
                if j+1 == len(patch_handles):
                    x = patch.get_width() + bl[0]
                    y = 0.5*patch.get_height() + bl[1]
                    ax.text(x,
                            y,
                            "%s" % (np.sum(section_labels[i]).round(decimals=3)),
                            ha='left',
                            bbox={"facecolor":"red", "alpha":.5})

        ax.set_yticks(y_pos)
        y_tick_labels = []
        for d in datasets:
            short_name = re.split(r'[\._]', d)[0]
            y_tick_labels.append(short_name+" THEIRS")
            y_tick_labels.append(short_name+" MINE")
        ax.set_yticklabels(y_tick_labels)
        ax.set_xlabel('Composite Score')
        plt.title('Composite Scores Using %s'%gold_standard.upper())
        plt.legend((patch_handles[0], patch_handles[1], patch_handles[2]),
                   ('Fraction Matched', 'Accuracy', 'Maximum Matching Ratio'),
                   loc='upper left',
                   framealpha=.25,
                   bbox_to_anchor=(1,1))

        plt.savefig(dir_path+"/"+gold_standard+"_bar_chart")

        plt.show()



    exit()




    pr = cProfile.Profile()
    pr.enable()
    # ... do something ...
    scores = CL1_Randomized("../cl1_datasets/datasets",
                       "collins2007.txt",
                       'Dummy_quality',
                       density_threshold=.15,
                       merge_threshold=.9,
                       penalty_value_per_node=10,
                       randomized_construction_bool=True,
                       rng_seed=None,
                       number_of_shakes=0,
                       number_of_bad_adds=2,
                       sort_seeds_by="weight",
                       care_about_cuts=False,
                       seed_from_all=True,
                       gsc_appearance_ratio_threshold=.9,
                       found_gsc_jaccard_threshold=.8,
                       gold_standard_filename="../cl1_gold_standard/gold_standard/mips_3_100.txt",
                       save_self=False)

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
    #                    density_threshold=.15,
    #                    merge_threshold=.9,
    #                    penalty_value_per_node=2,
    #                    randomized_construction_bool=False,
    #                    rng_seed=None,
    #                    number_of_shakes=0,
    #                    number_of_bad_adds=2,
    #                    sort_seeds_by="weight",
    #                    care_about_cuts=False,
    #                    seed_from_all=False,
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
"""
    189
    reference
    complexes, 303
    predicted
    complexes
    b'acc = 0.3883'
    b'cws = 0.3630'
    b'frac = 0.4127'
    b'mmr = 0.2211'
    b'ppv = 0.4155'
    b'sep = 0.2296'
    """

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