def cohesiveness(cl1, list_of_proteins, only_allow_good_queries = True) -> float:
    """Returns the cohesiveness of a potential complex.

    Searches the neighbors of each protein in the list, keeping track of weights of
    edges with both endpoints being internal to the complex, as well as weights of
    edges with one endpoint external to proposed complex.
    Then calculate cohesiveness of potential complex.

    :param list_of_proteins: list of protein id's (int) representing a potential complex
    :return: (float) the cohesiveness of the potential complex
    """
    weight_in = 0
    weight_out = 0
    for source in list_of_proteins:
        if not only_allow_good_queries and source not in cl1.graph.hash_graph:
            continue
        for target in cl1.graph.hash_graph[source]:
            edge_weight = cl1.graph.hash_graph[source][target]
            if target in list_of_proteins:
                weight_in += edge_weight/2.0
            else:
                weight_out += edge_weight
    # TODO: check that all is good with the calculation
    return weight_in / \
           ((weight_in + weight_out)+cl1.penalty_value_per_node*len(list_of_proteins))



def num_edges_inside(cl1, list_of_proteins):
    num_edges = 0
    for source in list_of_proteins:
        for target in cl1.graph.hash_graph[source]:
            if target in list_of_proteins:
                num_edges+=1
    num_edges/=2
    return num_edges


def density(cl1, list_of_proteins, only_allow_good_queries = True):
    # TODO: check that density is calculated with degree, not weights
    """Returns the density of a potential complex.

    Searches the neighbors of each protein in the list, keeping track of weights of
    edges with both endpoints being internal to the complex.
    Then calculate density of potential complex.

    :param list_of_proteins: list of protein id's (int) representing a potential complex
    :return: (float) the density of the potential complex
    """
    in_weight = 0
    for source in list_of_proteins:
        if not only_allow_good_queries and source not in cl1.graph.hash_graph:
            continue
        for target in cl1.graph.hash_graph[source]:
            if target in list_of_proteins:
                in_weight += cl1.graph.hash_graph[source][target]
    n = len(list_of_proteins)
    # TODO: check that authors didn't mean n* (n+1)/2
    # TODO: check that authors want to double count edges
    denominator = (n * (n - 1)) / 2
    return in_weight / denominator


def modularity(self, list_of_proteins):
    """NOT YET IMPLEMENTED

    :param list_of_proteins:
    :return:
    """
    return 1


def get_quality(gold_standard_filename,complexes_filename):
    import subprocess
    res = subprocess.check_output(["python2",
                                   "cl1_reproducibility/reproducibility/scripts/match_standalone.py",
                                   gold_standard_filename,
                                   complexes_filename])
    for line in res.splitlines():
        print(line)
        # a = str(line)
        # a = a.replace("b", "").replace("=", "").replace("\'", "").split()
        # self.quality_report[a[0]] = float(a[1])

