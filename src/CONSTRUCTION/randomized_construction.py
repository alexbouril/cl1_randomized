from src.COMMON.cmn import *
from src.CONSTRUCTION.initialize import *
from src.CLUSTER_STATE.cluster_state import *
from src.CONSTRUCTION.randomized_cluster_grower import rcg


def randomized_construction(cl1):
    """
    :param cl1:
    :return:
    """
    considered_vertices = set()
    index = 0
    if cl1.sort_seeds_by == 'degree':
        sorted_seeds = cl1.vertices_by_degree
    else:
        sorted_seeds = cl1.vertices_by_weight
    ############################################################
    # set the state of the random number generator
    ############################################################
    numpy.random.set_state(cl1.rng_initial_state)
    # cl1.logger.info(f'(RUN:{cl1.run_title}) rng_seed is {cl1.rng_seed}')
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
            # current_cluster_construction_log.append(
            #     ClusterState(current_cluster, add_candidates, remove_candidates, current_score))

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
            """
            Commented out code allows to see memory allocation
            """
            # tracemalloc.start()
            new_cluster = rcg(cl1, cs, current_cluster_construction_log)
            # snapshot = tracemalloc.take_snapshot()
            # top_stats = snapshot.statistics('lineno')
            # print("[ Top 10 ]")
            # for stat in top_stats[:10]:
            #     print(stat)
            ############################################################
            # add the cluster
            ############################################################
            cl1.initial_clustering.append(new_cluster)
            ############################################################
            # If we aren't seeding from each node, mark used nodes
            ############################################################
            if not cl1.seed_from_all:
                add_to_considered = set([v for v in cs.current_cluster])
                considered_vertices = considered_vertices.union(add_to_considered)
                considered_vertices.add(current_seed)
                print("len(considered_vertices): ",len(considered_vertices))
            index+=1

