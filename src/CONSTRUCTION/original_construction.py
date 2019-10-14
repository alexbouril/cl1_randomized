from src.CLUSTER_STATE.cluster_state import *
from src.CONSTRUCTION.original_cluster_grower import *
from src.CONSTRUCTION.initialize import *


def original_construction(cl1):
    considered_vertices = set()
    index = 0
    if cl1.sort_seeds_by == 'degree':
        sorted_seeds = cl1.vertices_by_degree
    else:
        sorted_seeds = cl1.vertices_by_weight

    while index < len(sorted_seeds):
        current_seed = sorted_seeds[index][0]
        if current_seed in considered_vertices:
            index += 1
            continue
        else:
            current_cluster_construction_log = []
            cl1.initial_clustering_seeds.append(current_seed)
            current_cluster_construction_log.append(Action("seed", current_seed))
            ############################################################
            # INITIALIZE THE CURRENT CLUSTER STATE STARTING WITH THE SELECTED SEED
            ############################################################
            current_cluster, remove_candidates, add_candidates, current_score, current_cluster_weight_in, current_cluster_weight_out = \
                initialize_complex(cl1, current_seed)
            current_cluster_construction_log.append(
                ClusterState(current_cluster, add_candidates, remove_candidates, current_score))

            # TODO: clean this block up.  These parameters are ignored in the original_construction method
            last_failed_add_round_no = -777
            last_failed_remove_round_no = -666
            round_no = 0
            local_number_of_shakes_remaining = cl1.number_of_shakes

            cs = ClusterState(current_cluster,
                              add_candidates,
                              remove_candidates,
                              current_score,
                              current_cluster_weight_in=current_cluster_weight_in,
                              current_cluster_weight_out=current_cluster_weight_out,
                              last_failed_add_round_no=last_failed_add_round_no,
                              last_failed_remove_round_no=last_failed_remove_round_no,
                              round_no=round_no,
                              number_of_shakes=local_number_of_shakes_remaining)

            ############################################################
            # Grow the cluster
            ############################################################
            ocg(cl1, cs, current_cluster_construction_log)
            index+=1
            ############################################################
            # If we aren't seeding from each node, mark used nodes
            ############################################################
            if not cl1.seed_from_all:
                add_to_considered = set([v for v in cs.current_cluster])
                considered_vertices = considered_vertices.union(add_to_considered)

