import collections

def merge(merge_threshold, source):  # takes a list of clusters
    hash_graph = collections.defaultdict(set)
    # TODO clean this up
    source = [set(x) for x in source]
    def dfs(index, local_visited):
        local_visited.add(index)
        for neigbor in hash_graph[index]:
            if neigbor not in local_visited:
                dfs(neigbor, local_visited)

    for i in range(len(source)):
        for j in range(i + 1, len(source)):
            if (len(source[i].intersection(source[j])) **2) / (len(source[i])*len(source[j])) > merge_threshold:
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