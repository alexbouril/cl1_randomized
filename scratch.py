from src.common.common import *

nice_comment("hello world.  what's going on?")
print(nice_comment("DETERMINE THE QUALITY OF THE RESULT USING ORIGINAL AUTHORS' MEASURES"))

exit()


# get_quality("cl1_gold_standard/gold_standard/mips_3_100.txt", "original_algo_findings.txt")
name = loadData('pickles/most_recent')['title']
print(name)
x = loadData(name)
print(len(x.graph.hash_graph))


def tie_print(g, li):
    for x in li:
        print(g.graph.id_to_name[x])
        print(g.graph.hash_graph[x])
    print("-----------------------------------")
#
# tie_print(x,[1533, 1534])
# tie_print(x,[980, 981])
# tie_print(x,[978, 980, 981])
# tie_print(x,[563, 569])
tie_print(x,[3046, 3048, 3049, 3050, 3051])


to_look_for = [x.graph.id_to_name[e] for e in [3046, 3048, 3049, 3050, 3051]]
f = open("cl1_datasets/datasets/krogan2006_extended.txt", "r")
for line in f.readlines():
    for e in line.split():
        if e in to_look_for:
            print(line)