def jaccard_similarity(l1:list, l2:list)->float:
    set1 = set(l1)
    set2 = set(l2)
    numerator = len(set1.intersection(set2))
    denominator = len(set1.union(set2))
    return numerator/denominator


def overlap_similarity(A, B):
    """Implements the overlap score described in the paper

    :param A: a set
    :param B: a set
    :return: the overlap score
    """
    numerator = len(A.intersection(B)) ** 2
    denominator = len(A) * len(B)
    return numerator / denominator