from common import *

def jaccard_similarity(l1:list, l2:list)->float:
    set1 = set(l1)
    set2 = set(l2)
    numerator = len(set1.intersection(set2))
    denominator = len(set1.union(set2))
    return numerator/denominator