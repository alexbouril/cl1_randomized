import time
import pickle
DEBUG = False
SLEEP_DEBUG = False
import os
import logging
import cProfile, pstats, io
from pstats import SortKey
import numpy
import heapq
import datetime

def profiler(cmd_str):
    pr = cProfile.Profile()
    pr.enable()
    exec(cmd_str)
    pr.disable()
    s = io.StringIO()
    sortby = SortKey.CUMULATIVE
    ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    ps.print_stats()
    print(s.getvalue())
#
# class TwoN_Add:
#     def __init(self, add_candidate, good_neighbors, prior_cohesiveness, proposed_score, post_cohesiveness):
#         self.add_candidate = add_candidate
#         self.good_neighbors = good_neighbors
#         self.prior_cohesiveness = prior_cohesiveness
#         self.proposed_score = proposed_score
#         self.post_cohesiveness = post_cohesiveness
#
# class checking_2n_add:
#     def __init__(self):
#         self.changes_made =[TwoN_Add]
#         self.final_cohesiveness = None



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
    if SLEEP_DEBUG:
        print(m)
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
    def __init__(self, _in:float, num_edges_to:int, _out:float, num_edges_from:int):
        self._in = _in
        self.num_edges_to = num_edges_to
        self._out = _out
        self.num_edges_from = num_edges_from

    def copy(self):
        return Relationship(self._in, self.num_edges_to, self._out, self.num_edges_from)

    def stringify(self):
        s=""
        s+="weight_to: "+ str(self._in) +"\t"
        s+="weight_from: "+ str(self._out)+"\n"
        s+="num_edges_to: "+ str(self.num_edges_to)+"\t"
        s+="num_edges_from: "+ str(self.num_edges_from)+"\n"
        return s


class ClusterState:
    def __init__(self,
                 current_cluster,
                 add_candidates,
                 remove_candidates,
                 cohesiveness,
                 best_change=None,
                 best_change_score=None,
                 current_cluster_weight_in=None,
                 current_cluster_weight_out=None,
                 neighborhood_2=None,
                 neighborhood_3=None,
                 last_failed_add_round_no = None,
                 last_failed_remove_round_no=None,
                 round_no=None,
                 number_of_shakes=None):
        self.current_cluster = current_cluster.copy()
        self.current_cluster = copy_relationship_dictionary(current_cluster)
        self.add_candidates = add_candidates.copy()
        self.add_candidates = {rel:add_candidates[rel].copy() for rel in add_candidates}
        self.remove_candidates = remove_candidates.copy()
        self.remove_candidates = {rel:remove_candidates[rel].copy() for rel in remove_candidates}
        self.cohesiveness = cohesiveness
        self.neighborhood_2 = neighborhood_2
        self.neighborhood_3 = neighborhood_3

        self.best_change=best_change
        self.best_change_score=best_change_score
        self.current_cluster_weight_in=current_cluster_weight_in
        self.current_cluster_weight_out=current_cluster_weight_out

        self.last_failed_add_round_no = last_failed_add_round_no
        self.last_failed_remove_round_no = last_failed_remove_round_no
        self.round_no = round_no
        self.local_number_of_shakes_remaining = number_of_shakes

    def make_backup(self):
        return ClusterState(self.current_cluster, self.add_candidates, self.remove_candidates, self.cohesiveness)

    def stringify_heavy(self, graph):
        s=">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n"
        s+=f"CURRENT_COHESIVENESS: {self.cohesiveness}\n"
        s += "\n"
        s+="**************************************\n"
        s += "********** CURRENT CLUSTER: **********\n"
        s+="**************************************\n"
        s += "\n"
        for protein in self.current_cluster:
            s+=str(protein)+", "+graph.id_to_name[protein]+"\n"
            s+=self.current_cluster[protein].stringify()
        s += "\n"
        s += "*************************************\n"
        s+="********** ADD CANDIDATES: **********\n"
        s+="*************************************\n"
        s += "\n"

        for protein in self.add_candidates:
            s+=str(protein)+", "+graph.id_to_name[protein]+"\n"
            s+=self.add_candidates[protein].stringify()
        s+="\n"
        s+="****************************************\n"
        s+="********** REMOVE CANDIDATES: **********\n"
        s+="****************************************\n"
        s += "\n"

        for protein in self.remove_candidates:
            s+=str(protein)+", "+graph.id_to_name[protein]+"\n"
            s+=self.remove_candidates[protein].stringify()
        s+=">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n"
        return s


    def stringify_lite(self):
        s=">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n"
        s+=f"CURRENT_COHESIVENESS: {self.cohesiveness}\n"
        s += "********** CURRENT CLUSTER: **********\n"
        s+= str([protein for protein in self.current_cluster])+"\n"
        s+="********** ADD CANDIDATES : **********\n"
        s+= str([protein for protein in self.add_candidates]) +"\n"
        s+="********** REMOVE CANDIDATES: **********\n"
        s+= str([protein for protein in self.remove_candidates]) +"\n"
        s+=">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n"
        return s





class Action:
    def __init__(self, action_name, involved_node = None):
        self.action_name = action_name
        self.involved_node  = involved_node

    def stringify(self):
        s=">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n"
        if self.action_name:
            s+=f"ACTION: {str(self.action_name)}\n"
        if self.involved_node:
            s+=f"INVOLVED NODE: {str(self.involved_node)}\n"
        s+=">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n"

        return s



def stringify_construction_log(construction_log, verbose = False, graph = None):
    s = ""
    for key in construction_log:
        s += "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ NEW CLUSTER @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
        s+=str(key)+"\n"
        for entry in construction_log[key]:
            if type(entry) == ClusterState:
                if verbose:
                    s += entry.stringify_heavy(graph)
                else:
                    s += entry.stringify_lite()
            if type(entry) == Action:
                s+=entry.stringify()
        s+= "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ END CLUSTER @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"

    return s

def stringify_single_cluster_construction_log(construction_log, cluster_key, verbose = False, graph = None):
    s = ""
    s += "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ NEW CLUSTER @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
    s += str(cluster_key) + "\n"
    for entry in construction_log[cluster_key]:
        if type(entry) == ClusterState:
            if verbose:
                s+= entry.stringify_heavy(graph)
            else:
                s += entry.stringify_lite()
        if type(entry) == Action:
            s += entry.stringify()
    s += "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ END CLUSTER @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
    return s


def copy_relationship_dictionary(rd):
    return {rel:rd[rel].copy() for rel in rd}


def get_quality(gold_standard_filename, output_filename):
    import subprocess
    res = subprocess.check_output(["python2",
                                   "cl1_reproducibility/reproducibility/scripts/match_standalone.py",
                                   gold_standard_filename,
                                   output_filename])
    for line in res.splitlines():
        print(line)
        # a = str(line)
        # a = a.replace("b", "").replace("=", "").replace("\'", "").split()


def nice_comment(comment):
    bound = ""
    for i in range(len(comment)+6):
        bound+="#"
    bound+="\n"
    s = bound + "#  " + comment + "  #\n" + bound
    return(s)