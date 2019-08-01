#exec(open("commands.py").read())
#exec(open("commands.py").read())
#python commands.py
from cl1_randomized import *
import matplotlib.pyplot as plt
plt.style.use('seaborn-whitegrid')
import numpy as np


def convert_cluster_list_to_subgraph(cluster_list):
    return 1

def graph_average_over_multiple_runs(number_of_runs):
    def namestr(obj, namespace=locals()):
        return str([name for name in namespace if namespace[name] is obj][0])

    number_of_complexes =[]
    avg_cohesiveness=[]
    avg_density=[]
    avg_acc = []
    avg_cws = []
    avg_frac = []
    avg_mmr = []
    avg_ppv = []
    avg_sep = []

    names_of_runs = []
    arguments = []
    for i in range(number_of_runs):
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
        number_of_complexes.append(len(a.final_clusters_stats['clusters']))
        avg_cohesiveness.append(a.final_clusters_stats["average_cohesiveness"])
        avg_density.append(a.final_clusters_stats["average_density"])
        avg_acc.append(a.quality_report['acc'])
        avg_cws.append(a.quality_report['cws'])
        avg_frac.append(a.quality_report['frac'])
        avg_mmr.append(a.quality_report['mmr'])
        avg_ppv.append(a.quality_report['ppv'])
        avg_sep.append(a.quality_report['sep'])
        names_of_runs.append(a.time_of_run)


    print(number_of_complexes)
    print(avg_cohesiveness)
    print(avg_density)
    print(names_of_runs)



    considerations = [avg_cohesiveness, avg_density,
                      avg_acc, avg_cws, avg_frac, avg_mmr,
                      avg_ppv, avg_sep]
    N = len(considerations)
    data_means=[]
    data_std =[]
    for x in considerations:
        data_means.append(np.mean(x))
        data_std.append(np.std(x))

    ind = np.arange(N)  # the x locations for the groups
    width = 0.35       # the width of the bars

    fig, ax = plt.subplots()
    rects1 = ax.bar(ind, data_means, width, color='r', yerr=data_std)

    women_means = tuple([.1 for i in range(N)])
    women_std = tuple([.05 for i in range(N)])
    rects2 = ax.bar(ind + width, women_means, width, color='b', yerr=women_std)

    # add some text for labels, title and axes ticks
    ax.set_ylabel('Scores')
    ax.set_title('Scores by group and gender')
    ax.set_xticks(ind + width / 2)
    ax.set_xticklabels(tuple([namestr(obj) for obj in considerations]), rotation='vertical')

    ax.legend((rects1[0], rects2[0]), ('Run Type 1', 'Run Type 2'))


    def autolabel(rects):
        """
        Attach a text label above each bar displaying its height
        """
        for rect in rects:
            height = rect.get_height()
            ax.text(rect.get_x() + rect.get_width()/2., 1.05*height,
                    '%d' % int(height),
                    ha='center', va='bottom')

    autolabel(rects1)
    autolabel(rects2)

    plt.show()






graph_average_over_multiple_runs(5)
