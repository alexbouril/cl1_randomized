
def sizeThreshold(source, size_thresh):
    retval = [cluster for cluster in source if
                                              len(cluster) > size_thresh]
    return retval


def densityThreshold(hash_graph, li, dens_thresh):
    def checkDensity(cluster_set):
        in_weight = 0
        for source in cluster_set:
            for target in hash_graph[source]:
                if target in cluster_set:
                    in_weight += hash_graph[source][target]
        n = len(cluster_set)
        # TODO: check that authors didn't mean n* (n+1)/2
        denominator = (n * (n - 1)) / 2
        return in_weight / denominator

    retval = []
    for cluster_set in li:
        density = checkDensity(cluster_set)
        if density > dens_thresh:
            retval.append(cluster_set)
    return retval