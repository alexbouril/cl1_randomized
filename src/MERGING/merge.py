import collections

def merge(merge_threshold, source):  # takes a list of clusters
    hash_graph = collections.defaultdict(set)

    def similarity(A, B):
        """Implements the overlap score described in the paper
        """
        # fst, snd = (A, B) if len(A) < len(B) else (B, A)
        # numerator = len(fst.intersection(snd)) **2
        numerator = len(A.intersection(B)) ** 2
        denominator = len(A) * len(B)
        return numerator / denominator

    def dfs(index, local_visited):
        local_visited.add(index)
        for neigbor in hash_graph[index]:
            if neigbor not in local_visited:
                dfs(neigbor, local_visited)

    for i in range(len(source)):
        for j in range(i + 1, len(source)):
            if similarity(source[i], source[j]) > merge_threshold:
                hash_graph[i].add(j)
                hash_graph[j].add(i)

    global_visited = set()
    merge_indices = list()
    for index in range(len(source)):
        if index not in global_visited:
            local_visited = set()
            dfs(index, local_visited)
            # TODO: is the copy necessary?
            # merge_indices.append(local_visited.copy())
            merge_indices.append(local_visited)
            for reached in local_visited:
                global_visited.add(reached)

    new_clusters = list()
    for group in merge_indices:
        new_cluster = set()
        for cluster_index in group:
            new_cluster = new_cluster.union(source[cluster_index])
        new_clusters.append(new_cluster)

    return new_clusters