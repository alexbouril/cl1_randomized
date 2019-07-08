import sys
import heapq
import numpy as np

GLOBAL_DEBUG = False

def debugging(*argv, debug):
    if debug or GLOBAL_DEBUG:
        for s in argv:
            print(s, end="")


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
                h[source].append([target, weight])
            else:
                h[source] = [[target, weight]]
            if target in h:
                h[target].append([source, weight])
            else:
                h[target] = [[source, weight]]
        return h, num_edges

    #  expects a file with every line as:
    # source_id_int target_id_int weight
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

def translate_complex_original_names_to_id(filename, new_filename, name_to_id):
    new_file = open(new_filename, "w+")
    old_file = open(filename, "r")
    s = " "
    for line in old_file.readlines():
        arr = line.split()
        new_arr = []
        for a in arr:
            new_arr.append(str(name_to_id[a]))
        new_file.write(s.join(new_arr)+"\n")
    old_file.close()
    new_file.close()


# input: file with protein ids, each line corresponds to a complex and is a list of space separated protein id numbers
# output: a list representing a solution, containing complexes objects
def get_complex_list_from_translated_file(filename):
    file = open(filename, "r")
    list_of_complexes = []
    for line in file.readlines():
        arr = [int(x) for x in line.split()]
        list_of_complexes.append(Complex(arr))
    return list_of_complexes


"""input:
        l1: a list of integers in sorted order
        l2: a list of integers in sorted order
    output: a float representing the jaccard similarity of the two lists
    """
def jaccard_similarity(l1:list, l2:list)->float:
    set1 = set(l1)
    set2 = set(l2)
    numerator = len(set1.intersection(set2))
    denominator = len(set1.union(set2))
    return numerator/denominator


# input:
#   list_protein_names: a list of protein names
# output:
#   name_to_id: a dict that has original names mapped to new protein id numbers
#   id_to_name: a dict that has protein id numbers mapped to a value of an original name
def proteinNamesToID(list_protein_names:list):
    name_to_id = {}
    id_to_name = {}
    for id, name in enumerate(list_protein_names):
        name_to_id[name] = id
        id_to_name[id] = name
    return name_to_id, id_to_name


# input:
#   protein_counts: a list, each index represents a protein, the value represents the number of its appearances in the solution
# output: a boolean which is True if each protein appears >=1 time, False otherwise
# TODO: extend feasibility to respresent varying fuzzy levels of feasibility
def feasibility(proteinCounts:list)->bool:
    numberMissing = 0
    for count in proteinCounts:
        if count <= 0:
            numberMissing+=1
    return not numberMissing


# input:
#   list_complexes: a list of Complex objects which represents a solution
#   num_proteins: the integer number of proteins in the PPIN
# output: array where the index in a protein ID, and value is the number of complexes the protein appears in
def protein_counts(list_complexes:list, num_proteins:int)->list:
    protein_counts = [0 for x in range(num_proteins)]
    for complex in list_complexes:
        for protein in complex.get_proteins():
            protein_counts[protein] += 1
    return protein_counts


# input:
#   previousProteinCounts: the list of protein counts of the previous solution
#   change_hash: a dict with a list of complexes to be removed at 'remove', & a list of complexes to be added at 'add'
# output: updated list of protein counts
def update_protein_counts(previousProteinCounts: list, change_hash:dict)->list:
    new_protein_counts = previousProteinCounts.copy()
    for complex_to_remove in change_hash['remove']:
        for protein in complex_to_remove:
            new_protein_counts[protein]-=1
    for complex_to_add in change_hash['add']:
        for protein in complex_to_add:
            new_protein_counts[protein] += 1
    return previousProteinCounts


class Complex:
    def __init__(self, li):
        self.proteins = sorted(li)

    def get_proteins(self):
        return self.proteins


class Solution:
    """input:
            score: a score corresponding to the quality of the current solution
            list_complexes: a list representing a solution, containing complexes objects
            protein_count: an array where the index in a protein ID, and value is the number of complexes the protein appears in
    """
    def __init__(self, score, soln_list, protein_count):
        self.score = score
        self.list_complexes = soln_list
        self.protein_count = protein_count

    def pprint_complexes(self):
        for c in self.list_complexes:
            print(c.proteins)


# TODO: make sure that heap only keeps track of 'num_to_heap' elements
    # say that we want to keep track of k best solutions on the Path
    # we can use two heaps (inspired by min heap / max heap / O(1)median) to keep track of top K solutions during traversal
    # then put the solutions into 1 regular heap at the end
# note that comparisons are being made by 'score', which is the user provided measure of quality
# note that there are changes in the class to make the min heap a max heap
class EliteSetHeap:
    def __init__(self, num_to_heap, li):
        self.num_to_heap = num_to_heap
        self.heap = heapq.heapify(li)

    def push(self, score_soln_list):
        negated_score = -1 * score_soln_list[0]
        soln = score_soln_list[1]
        to_push = [negated_score, soln]
        heapq.heappush(self.heap, to_push)

    def pop(self):
        val = heapq.heappop(self.heap)
        val[0] *= -1
        return val

    def top_k(self, k=1):
        top_list = heapq.nsmallest(k, self.heap)
        for element in top_list:
            element[0] *= -1
        return top_list


# TODO: implement quality functions
class Quality:
    quality_types = ['dummy1', 'dummy2', 'cohesiveness_weighted']

    def __init__(self, type):
        self.type = type
        if self.type not in self.quality_types:
            raise Exception('Invalid quality measure type %s entered.' % type)
        if self.type=="dummy1":
            self.quality_function = self.dummy1
        elif self.type=="dummy2":
            self.quality_function = self.dummy2
        elif self.type=="cohesiveness_weighted":
            self.quality_function = self.cohesiveness_weighted

# ''''input:
#         proposed_complex_list: a list of Complexes objects
#         ppin_graph: a _____ representation of the PPIN
#     output: a real number'''
    def dummy1(self, proposed_complex_list: list, ppin_graph) -> float:
        #raise Exception("Quality.dummy1 not implemented yet")
        return -1

    def dummy2(self, proposed_complex_list, ppin_graph):
        raise Exception("Quality.dummy2 not implemented yet")

    # TODO: Speed up
    # TODO: figure out why density results differ from the results given in ClusterONE with Cytoscape
    # TODO: consider moving cohesiveness score into Complex class, so that it is maintained and not recalculated from scratch
    # TODO: figure out what better weighted average scheme can be used
    def cohesiveness_weighted(self, proposed_complex_list, ppin_graph:Graph, debug_bool=False):
        numerator = 0
        denominator = 0
        for complex in proposed_complex_list:
            list_of_internal_proteins = complex.proteins
            list_of_external_proteins = []
            in_weight = 0
            out_weight = 0
            num_internal_edges = 0
            num_external_edges = 0
            for source in list_of_internal_proteins:
                for target, weight in ppin_graph.hash_graph[source]:
                    if target in list_of_internal_proteins:
                        in_weight+=weight
                        num_internal_edges+=1
                    else:
                        out_weight+=weight
                        num_external_edges+=1
                        if target not in list_of_external_proteins:
                            list_of_external_proteins.append(target)
            num_possible_edges = (len(complex.proteins) * (len(complex.proteins)+1) )/2
            num_internal_edges/=2
            in_weight/=2
            density = num_internal_edges/num_possible_edges
            complex_cohesiveness = in_weight / (in_weight + out_weight)
            numerator+=len(list_of_internal_proteins) * complex_cohesiveness
            denominator += len(list_of_internal_proteins)

            debugging(len(list_of_internal_proteins), "\t", density, "\t",complex_cohesiveness, "\n","========================",\
                      in_weight,"\t", out_weight, "\n","\n", debug=debug_bool)

        score = numerator / denominator
        print(score)
        return score




# TODO: implement other relinking strategies
class Strategy:
    strategy_types = ['forward', 'backward']

    def __init__(self, strategy_type):
        self.type = strategy_type
        if self.type not in self.strategy_types :
            raise Exception('Invalid strategy type %s entered.'%strategy_type)
        if self.type=="forward":
            self.strategy_function = self.forward
        elif self.type=="backward":
            self.strategy_function = self.backward

    def backward(self, soln_current, soln_guide, quality_object, heap_object):
        raise Exception("Strategy.backward not implemented yet")


    def match(self, soln_current:Solution, soln_guide:Solution):
        cur = soln_current.list_complexes
        gui = soln_guide.list_complexes
        # TODO: check to see which proteins should even be considered

        # TODO: match # of complexes

        len_cur = len(cur)
        len_gui = len(gui)
        if len_cur < len_gui:
            for i in range(len_gui-len_cur):
                cur.append([])
        else:
            for i in range(len_cur - len_gui):
                gui.append([])
        # TODO: maximum weighted matching

        # TODO: Move one protein at a time until we have relinked.
        #   Choose intermediate solutions greedily
        #       Perform Local Search from each intermediate solution
        #       Keep track of best solutions found
        return None

    def forward_2(self, soln_current: Solution, soln_guide: Solution, quality_object: Quality, heap_object: EliteSetHeap, ppin_graph)->EliteSetHeap:

        return None





    """input: 
            soln_current: a Solution object
            soln_guide: a Solution object
            quality_object: a Quality object
            heap_object: a Heap object
            ppin_graph: 
            num_proteins: an integer representing the number of proteins found in the graph
       output: a heap object
    """
    def forward(self, soln_current: Solution, soln_guide: Solution, quality_object: Quality, heap_object: EliteSetHeap, ppin_graph)->EliteSetHeap:
        C = set()       # a set of complexes common to soln_current and soln_guide
        G = set()       # a set of complexes found in soln_guide but not soln_current
        I = set()       # a set of complexes found in soln_current but no soln_guide
        for i in soln_current.list_complexes:
            for g in soln_guide.list_complexes:
                if jaccard_similarity(i, g) == 1:
                    C.add(i)
        for g in soln_guide.list_complexes:
            if not g in C:
                G.add(g)
        for i in soln_current.list_complexes:
            if not i in C:
                I.add(i)
        while I or G:
            candidates = []
            # -----------------------------------------------------
            # Find the best move that moves a complex from G into C
            # -----------------------------------------------------
            best_g_score = -1 * sys.maxsize
            best_g_move = None
            for g in G:
                # TODO: figure out if proposed should be a list or a set
                proposed = list(I.union(C.union(set([g]))))     # proposed is a list of complexes, NOT a Solution object
                proposed_protein_counts = soln_current.protein_count.copy()   # proposed_protein_counts is a protein count array corresponding to the proposed list of complexes
                proposed_quality = quality_object.quality_function(proposed, ppin_graph)
                for protein_number in g:
                    proposed_protein_counts[protein_number]+=1
                proposed_feasibility = feasibility(proposed_protein_counts)
                if proposed_quality * proposed_feasibility >= best_g_score:
                    best_g_score = proposed_quality
                    best_g_move = g
            # Add the Solution corresponding to the best 'G-move' to a list of candidates
            # TODO: figure out what should be copied
            option1_score = best_g_score.copy()
            option1_list_complexes = list(I.union(C.union(set([best_g_move]))))
            option1_change_hash = {'remove': [], 'add': [best_g_move]}
            option1_protein_counts =  update_protein_counts(soln_current.protein_count.copy(), option1_change_hash)
            option1 = Solution(option1_score, option1_list_complexes, option1_protein_counts)
            candidates.append((option1, option1_change_hash))

            # -----------------------------------------------------
            # Find the best move that removes a complex from I
            # -----------------------------------------------------
            best_i_score = -1 * sys.maxsize
            best_i_move = None
            for i in I:
                proposed = I.union(C)
                proposed.remove(i)
                proposed = list(proposed)   # proposed is a list of complexes
                proposed_protein_counts = soln_current.protein_count.copy()
                proposed_quality = quality_object.quality_function(proposed, ppin_graph)
                for protein_number in i:
                    proposed_protein_counts[protein_number]-=1
                proposed_feasibility = feasibility(proposed_protein_counts)
                if proposed_quality * proposed_feasibility >= best_i_score:
                    best_i_score = proposed_score
                    best_i_move = i
            # Add the Solution corresponding to the best 'I-move' to a list of candidates
            option2_score = best_i_score.copy()
            option2_list_complexes =  I.union(C)
            option2_list_complexes.remove(best_i_move)
            option2_list_complexes = list(option2_list_complexes)
            option2_change_hash =  {'remove': [best_i_move], 'add': []}
            option2_protein_counts =  update_protein_counts(soln_current.protein_count.copy(), option2_change_hash)
            option2 = Solution(option2_score, option2_list_complexes, option2_protein_counts)
            candidates.append((option2, option2_change_hash))

            # -----------------------------------------------------
            # Find the best replacement of a complex in I with one from G
            # -----------------------------------------------------
            best_ig_score = -1 * sys.maxsize
            best_ig_move = None
            for i in I:
                for g in G:
                    proposed = I.union(C).union(set([g]))
                    proposed.remove(i)
                    proposed = list(proposed.copy())
                    proposed_score = quality_object.quality_function(proposed, ppin_graph)
                    proposed_protein_counts = soln_current.protein_count.copy()
                    for protein_number in g:
                        proposed_protein_counts[protein_number] += 1
                    for protein_number in i:
                        proposed_protein_counts[protein_number] -= 1
                    proposed_feasibility = feasibility(proposed_protein_counts)
                    if proposed_score * proposed_feasibility >= best_ig_score:
                        best_ig_score =proposed_score
                        best_ig_move = [i, g]
            # Add the Solution corresponding to the best 'IG-move' to a list of candidates
            option3_score = best_ig_score.copy()
            option3_list_complexes =  I.union(C).union(set([best_ig_move[1]]))
            option3_list_complexes.remove(best_ig_move[0])
            option3_list_complexes = list(option2_list_complexes)
            option3_change_hash = {"remove": [best_ig_move[0]], 'add': [best_ig_move[1]]}
            option3_protein_counts = update_protein_counts(soln_current.protein_count.copy(), option3_change_hash)
            option3 = Solution(option3_score, option3_list_complexes, option3_protein_counts)
            candidates.append((option3, option3_change_hash))

            # Select the best move
            candidates = sorted(candidates, reverse=True, key=lambda x: x[0].score)
            best_solution = candidates[0][0]
            change_hash = candidates[0][1]

            # Update I, C, G properly
            for protein_complex in change_hash['remove']:
                I.remove(protein_complex)
            for protein_complex in change_hash['add']:
                C.add(protein_complex)
                G.remove(protein_complex)

            # Keep track of the best solutions seen so far
            # TODO: make sure that the copy method is working for the Solution object
            soln_copy = best_solution.copy()
            heap_object.push((soln_copy.score, soln_copy))

        return heap_object


class Relinking:
    """ input:
            initial_solution: a list representing a solution, containing complexes objects
            guide_solution: a list representing a solution, containing complexes objects
            num_to_heap:  an integer representing the maximum number of elements that the elite set should contain
            quality_measure_name: a string denoting one of the available options for quality measures (see 'Quality' class)
            num_proteins: an integer representing the number of proteins in the PPI network
            relinking_strategy: a string denoting one of the available options for the relinking strategy (see 'Strategy' class)

    """
    def __init__(self, base_path, ppin_graph_file, initial_solution_file, guide_solution_file, num_to_heap, quality_measure_name, relinking_strategy):
        ppin_graph_file = base_path + ppin_graph_file
        initial_solution_file = base_path + initial_solution_file
        guide_solution_file = base_path + guide_solution_file

        self.ppin_graph = Graph(ppin_graph_file)

        # create quality Quality object that has information about the quality measure used to measure fit of a solution
        self.quality = Quality(quality_measure_name)

        # get the initial protein counts and the initial score, then create the initial Solution object
        tranlated_initial_solution_filename = initial_solution_file.split(".")[0] + "_translated." + initial_solution_file.split(".")[1]
        translate_complex_original_names_to_id(initial_solution_file, tranlated_initial_solution_filename, self.ppin_graph.name_to_id)
        initial_solution = get_complex_list_from_translated_file(tranlated_initial_solution_filename)
        initial_protein_counts = protein_counts(initial_solution, self.ppin_graph.num_proteins)
        initial_score = self.quality.quality_function(initial_solution, self.ppin_graph)
        self.initial_solution = Solution(initial_score, initial_solution, initial_protein_counts)

        # get the guide protein_counts and the guide score, then create the guide Solution object
        tranlated_guide_solution_filename = guide_solution_file.split(".")[0] + "_translated." + guide_solution_file.split(".")[1]
        translate_complex_original_names_to_id(guide_solution_file, tranlated_guide_solution_filename, self.ppin_graph.name_to_id)
        guide_solution = get_complex_list_from_translated_file(tranlated_guide_solution_filename)
        guide_protein_counts = protein_counts(guide_solution, self.ppin_graph.num_proteins)
        guide_score = self.quality.quality_function(guide_solution, self.ppin_graph, debug_bool=True)
        self.guide_solution = Solution(guide_score, guide_solution, guide_protein_counts)

        # create the current Solution object
        self.current_solution = Solution(initial_score, initial_solution, initial_protein_counts)

        # create the EliteSetHeap object
        self.elite_set_heap = EliteSetHeap(num_to_heap, [self.initial_solution])

        # create the Strategy object
        self.relinking_strategy = Strategy(relinking_strategy)




    def relink(self):
        self.relinking_strategy.strategy_function(self.current_solution, self.guide_solution, self.quality, self.elite_set_heap, self.ppin_graph)

    def return_heap(self):
        return self.elite_set_heap.heap


    # TODO: test how weighted average quality corresponds to ground truth matching scores
    # TODO: Talk to Dr. Bui about how well PR works in light of the overlapping nature
    # TODO: Talk to Dr. Bui about the Nascimento paper
    #   should feasibility be maintained on path?
    #   should merges be allowed on path?
    #   should solution only contain complexes present in initial and guide?
    #       if we only take complexes present in initial and guide, it doesn't seem reasonable that you can expand any complex, because each is already locally optimal
    #   say that we have initial and guide solutions.  if we want to find a new (hopefully better) solution, is it reasonable to
    #       union the set of complexes from initial and guide
    #       


























