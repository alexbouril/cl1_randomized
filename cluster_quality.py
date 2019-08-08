from common import *


def cohesiveness(self, list_of_proteins) -> float:
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
        for target in self.graph.hash_graph[source]:
            if target in list_of_proteins:
                weight_in += self.graph.hash_graph[source][target] / 2.0
            else:
                weight_out += self.graph.hash_graph[source][target]
    # TODO: check that all is good with the calculation
    return weight_in / (weight_in + weight_out)


def density(self, list_of_proteins):
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
        for target in self.graph.hash_graph[source]:
            if target in list_of_proteins:
                in_weight += self.graph.hash_graph[source][target]
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