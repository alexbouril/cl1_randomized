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
import os
import logging
import pprint as pp

def setup_custom_logger(name):
    """ Setup logger """
    LOGGING_FILE_PATH = './logs/cl1_randomized.log'

    if not os.path.exists(os.path.dirname(LOGGING_FILE_PATH)):
        os.makedirs(os.path.dirname(LOGGING_FILE_PATH), exist_ok=True)

    logging.basicConfig(
        format='[%(asctime)s] %(name)s %(levelname)s (%(funcName)s) %(message)s',
        level=logging.DEBUG,
        handlers=[logging.StreamHandler(),
                  logging.FileHandler(LOGGING_FILE_PATH, encoding='UTF-8')]
    )

    return logging.getLogger(name)


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


def sort_vertices_by_degree(self) -> list:
    """
    :return: a sorted list of tuples corresponding to (vertex id, vertex degree)
    """
    retval = [[k, len(self.graph.hash_graph[k])] for k in self.graph.hash_graph]
    retval.sort(key=lambda x: x[1], reverse=True)
    return retval

def sort_vertices_by_weight(self) -> list:
    """
    :return: a sorted list of tuples corresponding to (vertex id, weight)
    """
    retval = [[k, sum([self.graph.hash_graph[k][target] for target in self.graph.hash_graph[k]])] for k in self.graph.hash_graph]
    retval.sort(key=lambda x: x[1], reverse=True)
    return retval

class Relationship:
    def __init__(self, sum_weight_to:float, num_edges_to:int, sum_weight_from:float, num_edges_from:int):
        self.sum_weight_to = sum_weight_to
        self.num_edges_to = num_edges_to
        self.sum_weight_from = sum_weight_from
        self.num_edges_from = num_edges_from

    def copy(self):
        return Relationship(self.sum_weight_to, self.num_edges_to, self.sum_weight_from, self.num_edges_from)

def jaccard_similarity(a,b):
    a = set(a.copy())
    b = set(b.copy())
    return len(a.intersection(b)) / len(a.union(b))