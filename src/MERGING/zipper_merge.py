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

    def find_v_weight_into_outOf_complex(v, complex):
        weight_in = 0
        weight_out = 0
        edge_weights = [cl1.graph.hash_graph[v][neighbor] for neighbor in cl1.graph.hash_graph[v]]
        for weight, neighbor in zip(edge_weights, cl1.graph.hash_graph[v]):
            if neighbor in complex:
                weight_in+=weight
            else:
                weight_out+=weight
        return weight_in, weight_out



    def find_fuse(A, B):
        """Implements the overlap score described in the paper
        """
        A_complex_set, A_neighbor_set = find_one_neighborhood_of_complex(A)
        B_complex_set, B_neighbor_set = find_one_neighborhood_of_complex(B)

        neighbor_intersection = A_neighbor_set.intersection(B_neighbor_set)

        complex_union = A_complex_set.union(B_complex_set)
        complex_union_list = list(complex_union)

        fuse_points = set()

        current_cohesiveness = cohesiveness(cl1, complex_union_list)
        for neighbor in neighbor_intersection:
            # c = cohesiveness(cl1, complex_union_list+[neighbor])
            # if c>current_cohesiveness:
            fuse_points.add(neighbor)

        if fuse_points and len(fuse_points)>=3:
            print(cohesiveness(cl1, list(A_complex_set)),
                  len(A_complex_set),
                  cohesiveness(cl1, list(B_complex_set)),
                  len(B_complex_set),
                  current_cohesiveness,
                  cohesiveness(cl1, list(A_complex_set)+list(B_complex_set)+list(fuse_points)))
            print("len(fuse_points)", len(fuse_points))
            print('hello')
            new_starting_complex_list = complex_union_list+list(fuse_points)
            print('beautiful')
            new_cluster_state = generate_cluster_state_given_current_cluster_list(cl1,new_starting_complex_list)
            print('world')
            # print(new_cluster_state.cohesiveness)
            pre = cohesiveness(cl1, [x for x in new_cluster_state.current_cluster])
            print("pre:",pre)

            new_cluster = [x for x in ocg(cl1, new_cluster_state, [])]
            post = cohesiveness(cl1, new_cluster)
            print('############## post:', post, len(new_cluster))

            if post>pre:
                print('^^^^^^^^^^^^^^^IMPROVEMENT^^^^^^^^^^^^^^^^^', post-pre)
            return new_cluster

    new_clusters = list()

    for i,a in enumerate(source):
        for j in range(i+1, len(source)):
            print(i,j)
            b = source[j]
            c = find_fuse(a,b)
            if c:
                new_clusters.append(c)

    return new_clusters