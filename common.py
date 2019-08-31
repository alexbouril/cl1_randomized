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

def jaccard_similarity(a,b):
    a = set(a.copy())
    b = set(b.copy())
    return len(a.intersection(b)) / len(a.union(b))


class Relationship:
    def __init__(self, sum_weight_to:float, num_edges_to:int, sum_weight_from:float, num_edges_from:int):
        self.sum_weight_to = sum_weight_to
        self.num_edges_to = num_edges_to
        self.sum_weight_from = sum_weight_from
        self.num_edges_from = num_edges_from

    def copy(self):
        return Relationship(self.sum_weight_to, self.num_edges_to, self.sum_weight_from, self.num_edges_from)

    def stringify(self):
        s=""
        s+="weight_to: "+ str(self.sum_weight_to) +"\n"
        s+="weight_from: "+ str(self.sum_weight_from)+"\n"
        s+="num_edges_to: "+ str(self.num_edges_to)+"\n"
        s+="num_edges_from: "+ str(self.num_edges_from)+"\n"
        return s


class ClusterState:
    def __init__(self, current_cluster, add_candidates, remove_candidates, cohesiveness, neighborhood_2=None, neighborhood_3=None):
        self.current_cluster = current_cluster.copy()
        self.add_candidates = add_candidates.copy()
        self.remove_candidates = remove_candidates.copy()
        self.cohesiveness = cohesiveness
        self.neighborhood_2 = neighborhood_2
        self.neighborhood_3 = neighborhood_3

    def stringify_heavy(self):
        s=">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n"
        s+=f"CURRENT_COHESIVENESS: {self.cohesiveness}\n"
        s += "CURRENT CLUSTER:\n"
        for protein in self.current_cluster:
            s+=str(protein)+"\n"
            s+=self.current_cluster[protein].stringify()
        s+="ADD CANDIDATES:\n"
        for protein in self.add_candidates:
            s+=str(protein)+"\n"
            s+=self.add_candidates[protein].stringify()
        s+="REMOVE CANDIDATES:\n"
        for protein in self.remove_candidates:
            s+=str(protein)+"\n"
            s+=self.remove_candidates[protein].stringify()
        s+="<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n"
        return s


    def stringify_lite(self):
        s=">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n"
        s+=f"CURRENT_COHESIVENESS: {self.cohesiveness}\n"
        s += "CURRENT CLUSTER:\n"
        s+= str([protein for protein in self.current_cluster])+"\n"
        s+="ADD CANDIDATES:\n"
        s+= str([protein for protein in self.add_candidates]) +"\n"
        s+="REMOVE CANDIDATES:\n"
        s+= str([protein for protein in self.remove_candidates]) +"\n"
        s+="<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n"
        return s





class Action:
    def __init__(self, action_name, involved_node = None):
        self.action_name = action_name
        self.involved_node  = involved_node

    def stringify(self):
        s=">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n"
        if self.action_name:
            s+=f"ACTION: {str(self.action_name)}\n"
        if self.involved_node:
            s+=f"INVOLVED NODE: {str(self.involved_node)}\n"
        s+="<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n"

        return s



def stringify_construction_log(construction_log):
    s = ""
    for key in construction_log:
        s += "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ NEW CLUSTER @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
        s+=str(key)+"\n"
        for entry in construction_log[key]:
            if type(entry) == ClusterState:
                s+=entry.stringify_lite()
            if type(entry) == Action:
                s+=entry.stringify()
                dummy = 0
        s+= "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ END CLUSTER @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"

    return s