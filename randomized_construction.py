from common import *
from construction_operations import *

def randomized_construction(self):
    """
    :param self:
    :return:
    """
    considered_vertices = set()
    index = 0
    if self.sort_seeds_by == 'degree':
        sorted_seeds = self.vertices_by_degree
    else:
        sorted_seeds = self.vertices_by_weight

    ############################################################
    # set the state of the random number generator
    ############################################################
    numpy.random.set_state(self.rng_initial_state)
    self.logger.info(f'(RUN:{self.run_title}) rng_seed is {self.rng_seed}')

    while index < len(sorted_seeds):
        current_seed = sorted_seeds[index][0]
        current_seed_degree = sorted_seeds[index][1]
        debug("current_seed: %s" %str(current_seed), "current_seed_degree: %s" %str(current_seed_degree))
        if current_seed in considered_vertices:
            debug("SKIP %s" % str(current_seed))
            index += 1
            continue
        else:
            current_cluster_construction_log = []
            self.initial_clustering_seeds.append(current_seed)
            current_cluster_construction_log.append(Action("seed", current_seed))
            debug("Starting cluster #%s" % str(len(self.initial_clustering)))
            # time.sleep(3)
            # TODO: ignore vertices that have been removed before during construction of current cluster
            ignore_vertices = set()

            ############################################################
            # INITIALIZE THE CURRENT CLUSTER STARTING WITH THE SELECTED SEED
            ############################################################
            current_cluster, remove_candidates, add_candidates, current_score, current_cluster_weight_in, current_cluster_weight_out = \
                initialize_complex(self, current_seed)

            current_cluster_construction_log.append(ClusterState(current_cluster, add_candidates, remove_candidates, current_score))

            last_failed_add_round_no = -777
            last_failed_remove_round_no = -666
            round_no = 0

            local_number_of_shakes_remaining = self.number_of_shakes

            while (add_candidates or remove_candidates) and abs(last_failed_remove_round_no - last_failed_add_round_no) != 1:
                debug("Current cluster #%s" % str(len(self.initial_clustering)))
                decider = numpy.random.rand()

                ############################################################
                # CONSIDER ADDING A VERTEX ON THE BOUNDARY
                ############################################################
                if (decider <= .5 or last_failed_remove_round_no == round_no) and last_failed_add_round_no != round_no:
                    round_no += 1
                    best_change, best_change_score = \
                        find_best_add(self, add_candidates, current_cluster, current_score, current_cluster_weight_in, current_cluster_weight_out)
                    if best_change:
                        current_score = best_change_score
                        current_cluster_weight_in, current_cluster_weight_out = \
                            add(self, add_candidates, current_cluster, remove_candidates, best_change, best_change_score, current_cluster_weight_in, current_cluster_weight_out)
                        current_cluster_construction_log.append(Action("add", best_change))
                        current_cluster_construction_log.append(
                            ClusterState(current_cluster, add_candidates, remove_candidates, current_score))

                    else:
                        debug("\n", "No improvement by ADDING", "\n")
                        current_cluster_construction_log.append(Action("failed to add"))
                        last_failed_add_round_no = round_no

                ############################################################
                # CONSIDER REMOVING A VERTEX ON THE BOUNDARY
                ############################################################
                if (decider > .5 or last_failed_add_round_no == round_no) and last_failed_remove_round_no != round_no:
                    round_no += 1
                    best_change, best_change_score = \
                        find_best_remove(self, remove_candidates, current_cluster, current_cluster_weight_in, current_cluster_weight_out, current_score)
                    if best_change:
                        current_score = best_change_score
                        current_cluster_weight_in, current_cluster_weight_out = \
                            remove(self, remove_candidates, add_candidates, current_cluster, best_change, best_change_score, current_cluster_weight_in, current_cluster_weight_out)
                        current_cluster_construction_log.append(Action("remove", best_change))
                        current_cluster_construction_log.append(
                            ClusterState(current_cluster, add_candidates, remove_candidates, current_score))
                    else:
                        debug("\n", "No improvement by REMOVING", "\n")
                        current_cluster_construction_log.append(Action("failed to remove"))
                        last_failed_remove_round_no = round_no

                ############################################################
                # IF STUCK IN LOCAL OPTIMUM, CONSIDER TAKING A 'BAD' ADD STEP
                ############################################################
                # if local_number_of_shakes_remaining and add_candidates and abs(last_failed_remove_round_no - last_failed_add_round_no) == 1:
                #     local_number_of_shakes_remaining -= 1
                #
                #     current_score, \
                #     current_cluster_weight_in, \
                #     current_cluster_weight_out, \
                #     round_no, \
                #     last_failed_add_round_no = \
                #         add_shake(self, add_candidates, remove_candidates, current_cluster, current_cluster_weight_in, current_cluster_weight_out, round_no, last_failed_add_round_no)
                debug("$$$$$$$$$", last_failed_add_round_no, last_failed_remove_round_no, decider)


            # add current_cluster to the list of clusters
            self.initial_clustering.append(current_cluster)
            self.construction_log[tuple([protein for protein in current_cluster])] = current_cluster_construction_log
            index += 1
            if not self.seed_from_all:
                considered_vertices.add(current_seed)
                for v in current_cluster:
                    considered_vertices.add(v)
            print("CLUSTER #%s: %s" % (str(len(self.initial_clustering)), str([vertex for vertex in current_cluster])))
