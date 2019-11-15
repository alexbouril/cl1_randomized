from src.QUALITY.quality import cohesiveness
from src.CLUSTER_STATE.cluster_state import generate_cluster_state_given_current_cluster_list, ClusterState
from src.CONSTRUCTION.randomized_cluster_grower import rcg
from src.CONSTRUCTION.original_cluster_grower import ocg
# def zipper_merger():  # takes a list of clusters
#     threshold = self.merge_threshold
#     hash_graph = dict()
#
#     def find_one_neighborhood_of_complex(complex):
#         complex_set = set([v for v in complex])
#         neighbor_set = set()
#         for v in complex_set:
#             for u in self.graph.hash_graph[v]:
#                 if u not in complex_set:
#                     neighbor_set.add(u)
#         return complex_set, neighbor_set
#
#     def similarity(A, B):
#         A_complex_set, A_neighbor_set = find_one_neighborhood_of_complex(A)
#         B_complex_set, B_neighbor_set = find_one_neighborhood_of_complex(B)
#
#         """Implements the overlap score described in the paper
#
#         :param A: a set
#         :param B: a set
#         :return: the overlap score
#         """
#         numerator = len(A.intersection(B)) ** 2
#         denominator = len(A) * len(B)
#         return numerator / denominator
#
#     def dfs(index, local_visited):
#         local_visited.add(index)
#         for neigbor in hash_graph[index]:
#             if neigbor not in local_visited:
#                 dfs(neigbor, local_visited)
#
#     for i in range(len(self.initial_cluster_list)):
#         if i not in hash_graph:
#             hash_graph[i] = set()
#         for j in range(i + 1, len(self.initial_cluster_list)):
#             if similarity(self.initial_cluster_list[i], self.initial_cluster_list[j]) > threshold:
#                 if i in hash_graph:
#                     hash_graph[i].add(j)
#                 else:
#                     hash_graph[i] = {j}
#                 if j in hash_graph:
#                     hash_graph[j].add(i)
#                 else:
#                     hash_graph[j] = {i}
#
#     global_visited = set()
#     merge_indices = list()
#     for index in range(len(self.initial_cluster_list)):
#         if index not in global_visited:
#             local_visited = set()
#             dfs(index, local_visited)
#             # TODO: is the copy necessary?
#             merge_indices.append(local_visited.copy())
#             for reached in local_visited:
#                 global_visited.add(reached)
#
#     new_clusters = list()
#     for group in merge_indices:
#         new_cluster = set()
#         for cluster_index in group:
#             new_cluster = new_cluster.union(self.initial_cluster_list[cluster_index])
#         new_clusters.append(new_cluster)
#
#     self.merged_cluster_list = new_clusters


def zipper_merge(cl1, source):  # takes a list of clusters
    hash_graph = dict()

    def find_one_neighborhood_of_complex(complex):
        complex_set = set([v for v in complex])
        neighbor_set = set()
        for v in complex_set:
            for u in cl1.graph.hash_graph[v]:
                if u not in complex_set:
                    neighbor_set.add(u)
        return complex_set, neighbor_set

    def find_fuse(i, j):
        """Implements the overlap score described in the paper
        """
        A_complex_set, A_neighbor_set = neighborhoods_of_sources[i]
        B_complex_set, B_neighbor_set = neighborhoods_of_sources[j]
        fuse_points = A_neighbor_set.intersection(B_neighbor_set)
        if fuse_points and len(fuse_points)>=7:
            complex_union = A_complex_set.union(B_complex_set)
            complex_union_list = list(complex_union)
            new_starting_complex_list = complex_union_list+list(fuse_points)
            new_cluster_state = generate_cluster_state_given_current_cluster_list(cl1,new_starting_complex_list)
            pre = new_cluster_state.best_change_score
            new_cluster = ocg(cl1, new_cluster_state, [])
            new_cluster_list = [x for x in new_cluster]
            post = new_cluster_state.best_change_score
            #TODO use weighted average of cohesiveness scores of original clusters A, B to compare post against
            threshold = (cohesiveness_of_sources[i] * len(source[i]) + cohesiveness_of_sources[j] * len(source[j]))\
                        /(len(source[i])+len(source[j]))
            if post>threshold*3:
                print("pre:", pre)
                print("cohesiveness(A):", cohesiveness_of_sources[i])
                print("cohesiveness(B):", cohesiveness_of_sources[j])
                print("threshold:", threshold)
                print("|A|:", len(A_complex_set), "|B|:", len(B_complex_set), "len(fuse_points)", len(fuse_points),
                      "|new_cluster_list|:", len(new_cluster_list))
                print('############## post:', post)
                print('^^^^^^^^^^^^^^^IMPROVEMENT^^^^^^^^^^^^^^^^^', post-threshold, post -pre)
                retval = new_cluster_list
            else:
                retval = None
        else:
            retval = None
        return retval

    new_clusters = list()
    source = [x for x in source if len(x)<=20]
    neighborhoods_of_sources = list(map(find_one_neighborhood_of_complex,source))
    cohesiveness_of_sources = [cohesiveness(cl1, x) for x in source]
    for i in range(200): #len(source)):
        # TODO calculate partial starting cluster state
        print("##################################################################### ", i)
        for j in range(i+1, 100):#len(source)):
            c = find_fuse(i,j)
            if c:
                new_clusters.append(c)
                print("len(new_clusters):",len(new_clusters))

    return new_clusters