from cl1_randomized import *
from common import loadData
name = \
"pickles/pickle+gavin2006_socioaffinities_rescaled+2019-07-31_22:23:49:185488"





x = loadData(name)
print(x.quality_report)
print(x.initial_clustering)
print(x.final_clusters_stats)
print(x.found)
print(x.gsc_appearing_found_stats)
print(x.gsc_appearing_notFound_stats)
print(x.gsc_appearing_stats)
for d in [x.final_clusters_stats, x.gsc_appearing_stats, x.gsc_appearing_found_stats, x.gsc_appearing_notFound_stats]:
    print(d['average_cohesiveness'], d['average_density'], d['average_size'])
    print()