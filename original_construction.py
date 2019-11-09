from src.CONSTRUCTION.initialize import *


def oc(self):
    considered_vertices = set()
    index = 0
    if self.sort_seeds_by == 'degree':
        sorted_seeds = self.vertices_by_degree
    else:
        sorted_seeds = self.vertices_by_weight

    while index < len(sorted_seeds):
        current_seed = sorted_seeds[index][0]
        current_seed_degree = sorted_seeds[index][1]
        debug("current_seed: %s" %str(current_seed), "current_seed_degree: %s" %str(current_seed_degree))
        if current_seed in considered_vertices:
            debug("SKIP %s" % str(current_seed))
            index += 1
            continue
        else:
            end_construction = False
            current_cluster_construction_log = []
            self.initial_clustering_seeds.append(current_seed)
            current_cluster_construction_log.append(Action("seed", current_seed))
            debug("Starting cluster #%s" % str(len(self.initial_clustering)))
            # time.sleep(3)
            # TODO: ignore vertices that have been removed before during CONSTRUCTION of current cluster
            ignore_vertices = set()

            ############################################################
            # INITIALIZE THE CURRENT CLUSTER STARTING WITH THE SELECTED SEED
            ############################################################
            current_cluster, remove_candidates, add_candidates, current_score, current_cluster_weight_in, current_cluster_weight_out = \
                initialize_complex(self, current_seed)
            # current_cluster_construction_log.append(ClusterState(current_cluster, add_candidates, remove_candidates, current_score))

            while (add_candidates or remove_candidates) and not end_construction:
                debug("Current cluster #%s" % str(len(self.initial_clustering)))

                best_add_change, best_add_change_score = \
                    find_best_add(self, add_candidates, current_cluster, current_score,
                                                current_cluster_weight_in, current_cluster_weight_out)
                best_remove_change, best_remove_change_score = \
                    find_best_remove(self, remove_candidates, current_cluster, current_cluster_weight_in,
                                     current_cluster_weight_out, current_score)
                print(best_add_change, best_add_change_score,best_remove_change,best_remove_change_score, len(current_cluster), len(add_candidates))
                ############################################
                # ################
                # END CURRENT CLUSTER GROWTH IF NECESSARY
                ############################################################
                if (best_add_change is None) and (best_remove_change is None):
                    end_construction = True
                else:
                    ############################################################
                    # ADD A VERTEX ON THE BOUNDARY
                    ############################################################
                    if best_add_change_score > best_remove_change_score:
                            current_score = best_add_change_score
                            current_cluster_weight_in, current_cluster_weight_out = \
                                add(self,
                                    add_candidates,
                                    current_cluster,
                                    remove_candidates,
                                    best_add_change,
                                    best_add_change_score,
                                    current_cluster_weight_in,
                                    current_cluster_weight_out)
                            current_cluster_construction_log.append(Action("adding", best_add_change))
                            # current_cluster_construction_log.append(
                            #     ClusterState(current_cluster, add_candidates, remove_candidates, current_score))

                    ############################################################
                    # REMOVE A VERTEX ON THE BOUNDARY
                    ############################################################
                    else:
                        print("==================================================== REMOVING =====================================================")
                        current_score = best_remove_change_score
                        current_cluster_weight_in, current_cluster_weight_out = \
                            remove(self, remove_candidates, add_candidates, current_cluster, best_remove_change, best_remove_change_score, current_cluster_weight_in, current_cluster_weight_out)
                        current_cluster_construction_log.append(Action("removing", best_remove_change))
                        # current_cluster_construction_log.append(
                        #     ClusterState(current_cluster, add_candidates, remove_candidates, current_score))

            # add current_cluster to the list of clusters
            self.initial_clustering.append(current_cluster)
            self.construction_log[tuple([protein for protein in current_cluster])] = current_cluster_construction_log
            index += 1
            if not self.seed_from_all:
                considered_vertices.add(current_seed)
                for v in current_cluster:
                    considered_vertices.add(v)
            print("CLUSTER #%s: %s" % (str(len(self.initial_clustering)), str([vertex for vertex in current_cluster])))
