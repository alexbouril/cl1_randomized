import sys
import heapq
import numpy as np
import pickle
import time
import pickle
import datetime
DEBUG = False
SLEEP_DEBUG = False
import numpy
from cluster_quality import cohesiveness, density



def debug(*argv):
    if DEBUG:
        for arg in argv:
            print(arg)

def sleep_debug(t, m = ""):
    print(m)
    if SLEEP_DEBUG:
        time.sleep(t)

def loadData(f_name):
    # for reading also binary mode is important
    f = open(f_name, 'rb')
    obj = pickle.load(f)
    f.close()
    return obj

def jaccard_similarity(l1:list, l2:list)->float:
    set1 = set(l1)
    set2 = set(l2)
    numerator = len(set1.intersection(set2))
    denominator = len(set1.union(set2))
    return numerator/denominator


class Relationship:
    def __init__(self, sum_weight_to:float, num_edges_to:int, sum_weight_from:float, num_edges_from:int):
        self.sum_weight_to = sum_weight_to
        self.num_edges_to = num_edges_to
        self.sum_weight_from = sum_weight_from
        self.num_edges_from = num_edges_from

    def copy(self):
        return Relationship(self.sum_weight_to, self.num_edges_to, self.sum_weight_from, self.num_edges_from)

