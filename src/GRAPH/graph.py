from src.COMMON.cmn import *


def dfs(cl1, current_vertex, ignore_vertex, current_cluster_membership_hashset, visited):
    visited.add(current_vertex)
    for neighbor in cl1.graph.hash_graph[current_vertex]:
        if neighbor not in current_cluster_membership_hashset:
            continue
        elif neighbor in visited:
            continue
        elif neighbor == ignore_vertex:
            continue
        else:
            dfs(cl1, neighbor, ignore_vertex, current_cluster_membership_hashset, visited)

class Graph:
    #  expects a file with every line as:
    # source_name target_name weight
    def __init__(self, original_filename):
        self.original_filename = original_filename
        self.list_of_names = self.get_list_of_names()
        self.num_proteins = len(self.list_of_names)
        self.name_to_id, self.id_to_name  = self.original_protein_names_to_integer_id()
        self.new_filename = self.original_filename.split(".")[0]+"_translated."+self.original_filename.split(".")[1]
        self.translate_graph_original_names_to_id()
        self.hash_graph, self.num_edges = self.create_hash_graph_weighted()

    def create_hash_graph_weighted(self):
        h={}
        f = open(self.new_filename, "r")
        num_edges = 0
        for line in f.readlines():
            num_edges+=1
            source, target, weight = line.split()
            source = int(source)
            target = int(target)
            weight = float(weight)
            if source in h:
                h[source][target] = weight
            else:
                h[source] = {target: weight}
            if target in h:
                h[target][source] = weight
            else:
                h[target] = {source: weight}
        return h, num_edges

    #  expects a file with every line as:
    #       source_id_int target_id_int weight
    def get_list_of_names(self):
        f = open(self.original_filename, "r")
        li = []
        for line in f.readlines():
            source, target, weight = line.split()
            if source not in li:
                li.append(source)
            if target not in li:
                li.append(target)
        return li

    # input:
    #   list_protein_names: a list of protein names
    # output:
    #   name_to_id: a dict that has original names mapped to new protein id numbers
    #   id_to_name: a dict that has protein id numbers mapped to a value of an original name
    def original_protein_names_to_integer_id(self):
        name_to_id = {}
        id_to_name = {}
        for id, name in enumerate(self.list_of_names):
            name_to_id[name] = id
            id_to_name[id] = name
        return name_to_id, id_to_name

    def translate_graph_original_names_to_id(self):
        new_file = open(self.new_filename,"w+")
        old_file = open(self.original_filename, "r")
        for line in old_file.readlines():
            source, target, weight = line.split()
            source = self.name_to_id[source]
            target = self.name_to_id[target]
            new_file.write(str(source)+"\t"+str(target)+"\t"+weight+"\n")
        old_file.close()
        new_file.close()