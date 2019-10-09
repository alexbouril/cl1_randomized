# import src.cl1_randomized.cl1_randomized
# from src.common.common import *
#
# def randomized_cluster_grower(locals):
#     current_cluster_construction_log = []
#     self.initial_clustering_seeds.append(current_seed)
#     current_cluster_construction_log.append(Action("seed", current_seed))
#     debug("Starting cluster #%s" % str(len(self.initial_clustering)))
#     # time.sleep(3)
#     # TODO: ignore vertices that have been removed before during construction of current cluster
#     ignore_vertices = set()
#
#     ############################################################
#     # INITIALIZE THE CURRENT CLUSTER STARTING WITH THE SELECTED SEED
#     ############################################################
#     current_cluster, remove_candidates, add_candidates, current_score, current_cluster_weight_in, current_cluster_weight_out = \
#         initialize_complex(self, current_seed)
#
#     current_cluster_construction_log.append(
#         ClusterState(current_cluster, add_candidates, remove_candidates, current_score))
#
#     last_failed_add_round_no = -777
#     last_failed_remove_round_no = -666
#     round_no = 0
#
#     ############################################################
#     # Set the number of shakes available to the current cluster construction
#     ############################################################
#     local_number_of_shakes_remaining = self.number_of_shakes
#
#     ############################################################
#     # initialize the backup in case add_shake gets us stuck in a lower local optima
#     ############################################################
#     backup_current_cluster = current_cluster.copy()
#     backup_add_candidates = add_candidates.copy()
#     backup_remove_candidates = remove_candidates.copy()
#     backup_current_score = current_score
#     remove_counter = 0
#     while (add_candidates or remove_candidates) and \
#             abs(last_failed_remove_round_no - last_failed_add_round_no) != 1 \
#             and remove_counter < 100:
#         debug("Current cluster #%s" % str(len(self.initial_clustering)))
#         decider = numpy.random.rand()
#
#         ############################################################
#         # CONSIDER ADDING A VERTEX ON THE BOUNDARY
#         ############################################################
#         if (decider <= .5 or last_failed_remove_round_no == round_no) and last_failed_add_round_no != round_no:
#             # #@#
#             # best_change, best_change_score = \
#             #     find_best_add(self, add_candidates, current_cluster, current_score,
#             #                   current_cluster_weight_in, current_cluster_weight_out)
#             # #@#
#             # if best_change is None:
#             if len(current_cluster) > 5 and round_no % 4 == 1:
#                 round_no += 1
#                 ######################################################################
#                 # JUST FOR COMPARISON
#                 #####################################################################
#                 # @#
#                 best_change, best_change_score = \
#                     find_best_add(self, add_candidates, current_cluster, current_score,
#                                   current_cluster_weight_in, current_cluster_weight_out)
#                 if best_change is None:
#                     print("NO REGULAR ADD AVAILABLE")
#                 if best_change:
#                     print("REGULAR BEST CHANGE AVAILABLE WITH SCORE OF: %s" % str(best_change_score))
#                 ######################################################################
#                 # NOW FOR THE REAL CALL
#                 #####################################################################
#                 best_change, best_change_score = \
#                     careful_find_best_2neighborhood_add(self, add_candidates, current_cluster, current_score,
#                                                         current_cluster_weight_in, current_cluster_weight_out)
#                 # find_best_add(self, add_candidates, current_cluster, current_score, current_cluster_weight_in, current_cluster_weight_out)
#
#                 if best_change:
#
#                     print("*******", current_score)
#                     ##################################
#                     # sanity check
#                     ##################################
#                     # new_cluster = [v for v in current_cluster]+[best_change]
#                     # current_score_check = cohesiveness(self, new_cluster)
#                     # if abs(best_change_score-current_score_check)>.001:
#                     #     print(abs(best_change_score-current_score_check), best_change_score, current_score_check)
#
#                     print("adding: ", best_change)
#                     current_score = best_change_score
#                     current_cluster_weight_in, current_cluster_weight_out = \
#                         add(self, add_candidates, current_cluster, remove_candidates, best_change,
#                             best_change_score, current_cluster_weight_in, current_cluster_weight_out)
#                     current_cluster_construction_log.append(Action("adding", best_change))
#                     current_cluster_construction_log.append(
#                         ClusterState(current_cluster, add_candidates, remove_candidates, current_score))
#                     print(current_score)
#                     #################################################
#                     # update the backup
#                     #################################################
#                     if current_score > backup_current_score:
#                         backup_current_cluster = current_cluster.copy()
#                         backup_add_candidates = add_candidates.copy()
#                         backup_remove_candidates = remove_candidates.copy()
#                         backup_current_score = current_score
#
#                     for counter in range(4):
#                         round_no += 1
#                         best_change, best_change_score = \
#                             find_best_suboptimal_add(self,
#                                                      add_candidates,
#                                                      current_cluster,
#                                                      current_cluster_weight_in,
#                                                      current_cluster_weight_out)
#                         if best_change:
#                             current_score = best_change_score
#                             current_cluster_weight_in, current_cluster_weight_out = \
#                                 add(self, add_candidates, current_cluster, remove_candidates, best_change,
#                                     best_change_score, current_cluster_weight_in, current_cluster_weight_out)
#                             current_cluster_construction_log.append(Action("adding", best_change))
#                             current_cluster_construction_log.append(
#                                 ClusterState(current_cluster, add_candidates, remove_candidates, current_score))
#                             #################################################
#                             # update the backup
#                             #################################################
#                             if current_score > backup_current_score:
#                                 backup_current_cluster = current_cluster.copy()
#                                 backup_add_candidates = add_candidates.copy()
#                                 backup_remove_candidates = remove_candidates.copy()
#                                 backup_current_score = current_score
#                         else:
#                             debug("\n", "No improvement by ADDING", "\n")
#                             current_cluster_construction_log.append(
#                                 Action("failed to add, current cohesiveness: %s" % str(current_score)))
#                             last_failed_add_round_no = round_no
#                     for counter in range(2):
#                         round_no += 1
#                         best_change, best_change_score = \
#                             find_best_add(self, add_candidates, current_cluster, current_score,
#                                           current_cluster_weight_in, current_cluster_weight_out)
#                         if best_change:
#                             current_score = best_change_score
#                             current_cluster_weight_in, current_cluster_weight_out = \
#                                 add(self, add_candidates, current_cluster, remove_candidates, best_change,
#                                     best_change_score, current_cluster_weight_in, current_cluster_weight_out)
#                             current_cluster_construction_log.append(Action("adding", best_change))
#                             current_cluster_construction_log.append(
#                                 ClusterState(current_cluster, add_candidates, remove_candidates, current_score))
#                             #################################################
#                             # update the backup
#                             #################################################
#                             if current_score > backup_current_score:
#                                 backup_current_cluster = current_cluster.copy()
#                                 backup_add_candidates = add_candidates.copy()
#                                 backup_remove_candidates = remove_candidates.copy()
#                                 backup_current_score = current_score
#
#                         else:
#                             debug("\n", "No improvement by ADDING", "\n")
#                             current_cluster_construction_log.append(
#                                 Action("failed to add, current cohesiveness: %s" % str(current_score)))
#                             last_failed_add_round_no = round_no
#
#                         print(".......................", current_score)
#
#                 else:
#                     debug("\n", "No improvement by ADDING", "\n")
#                     current_cluster_construction_log.append(
#                         Action("failed to add, current cohesiveness: %s" % str(current_score)))
#                     last_failed_add_round_no = round_no
#             else:
#                 round_no += 1
#                 # @#
#                 best_change, best_change_score = \
#                     find_best_add(self, add_candidates, current_cluster, current_score,
#                                   current_cluster_weight_in, current_cluster_weight_out)
#
#                 if best_change:
#                     print("-------------------------------------------------------------", current_score)
#                     current_score = best_change_score
#                     print("-<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<", current_score)
#
#                     current_score_check = cohesiveness(self, [v for v in current_cluster] + [best_change])
#                     ##################################
#                     # sanity check
#                     ##################################
#                     if abs(best_change_score - current_score_check) > .001:
#                         print("----------------", current_score, current_score_check)
#                     current_cluster_weight_in, current_cluster_weight_out = \
#                         add(self, add_candidates, current_cluster, remove_candidates, best_change,
#                             best_change_score, current_cluster_weight_in, current_cluster_weight_out)
#                     current_cluster_construction_log.append(Action("adding", best_change))
#                     current_cluster_construction_log.append(
#                         ClusterState(current_cluster, add_candidates, remove_candidates, current_score))
#                     #################################################
#                     # update the backup
#                     #################################################
#                     if current_score > backup_current_score:
#                         backup_current_cluster = current_cluster.copy()
#                         backup_add_candidates = add_candidates.copy()
#                         backup_remove_candidates = remove_candidates.copy()
#                         backup_current_score = current_score
#
#                 else:
#                     debug("\n", "No improvement by ADDING", "\n")
#                     current_cluster_construction_log.append(
#                         Action("failed to add, current cohesiveness: %s" % str(current_score)))
#                     last_failed_add_round_no = round_no
#
#         ############################################################
#         # CONSIDER REMOVING A VERTEX ON THE BOUNDARY
#         ############################################################
#         if (decider > .5 or last_failed_add_round_no == round_no) and last_failed_remove_round_no != round_no:
#             round_no += 1
#             best_change, best_change_score = \
#                 find_best_remove(self, remove_candidates, current_cluster, current_cluster_weight_in,
#                                  current_cluster_weight_out, current_score)
#             if best_change:
#                 remove_counter += 1
#                 print(
#                     "==================================================== REMOVING %s=====================================================" % str(
#                         best_change))
#                 current_score = best_change_score
#                 current_cluster_weight_in, current_cluster_weight_out = \
#                     remove(self, remove_candidates, add_candidates, current_cluster, best_change, best_change_score,
#                            current_cluster_weight_in, current_cluster_weight_out)
#                 current_cluster_construction_log.append(Action("removing", best_change))
#                 current_cluster_construction_log.append(
#                     ClusterState(current_cluster, add_candidates, remove_candidates, current_score))
#                 #################################################
#                 # update the backup
#                 #################################################
#                 if current_score > backup_current_score:
#                     backup_current_cluster = current_cluster.copy()
#                     backup_add_candidates = add_candidates.copy()
#                     backup_remove_candidates = remove_candidates.copy()
#                     backup_current_score = current_score
#             else:
#                 debug("\n", "No improvement by REMOVING", "\n")
#                 current_cluster_construction_log.append(
#                     Action("failed to remove, current cohesiveness: %s" % str(current_score)))
#                 last_failed_remove_round_no = round_no
#
#         ############################################################
#         # IF STUCK IN LOCAL OPTIMUM, CONSIDER TAKING A 'BAD' ADD STEP
#         ############################################################
#         if local_number_of_shakes_remaining and add_candidates and abs(
#                 last_failed_remove_round_no - last_failed_add_round_no) == 1:
#             local_number_of_shakes_remaining -= 1
#             current_cluster_construction_log.append(
#                 Action("trying add_shake, current cohesiveness: %s" % str(current_score)))
#             current_score, \
#             current_cluster_weight_in, \
#             current_cluster_weight_out, \
#             round_no, \
#             last_failed_add_round_no = \
#                 add_shake(self,
#                           add_candidates,
#                           remove_candidates,
#                           current_cluster,
#                           current_cluster_weight_in,
#                           current_cluster_weight_out,
#                           round_no,
#                           last_failed_add_round_no)
#             current_cluster_construction_log.append(
#                 ClusterState(current_cluster,
#                              add_candidates,
#                              remove_candidates,
#                              current_score))
#             # TODO: improve the granularity of the backup for add_shake
#             #   what if in the middle of add shake we see an improvement that is lost by the end of add_shake
#             #################################################
#             # update the backup
#             #################################################
#             if current_score > backup_current_score:
#                 backup_current_cluster = current_cluster.copy()
#                 backup_add_candidates = add_candidates.copy()
#                 backup_remove_candidates = remove_candidates.copy()
#                 backup_current_score = current_score
#
#         debug("$$$$$$$$$", last_failed_add_round_no, last_failed_remove_round_no, decider)
#
#     #########################################
#     # in the case that the current_cluster is not the best one that we saw, revert to the best one that we saw
#     ##########################################
#     if backup_current_score > current_score:
#         current_cluster_construction_log.append(Action("reverting to previous state"))
#         current_cluster = backup_current_cluster.copy()
#         add_candidates = backup_add_candidates.copy()
#         remove_candidates = backup_remove_candidates.copy()
#         current_score = backup_current_score
#         current_cluster_construction_log.append(
#             ClusterState(current_cluster, add_candidates, remove_candidates, current_score))
#
#     # add current_cluster to the list of clusters
#     self.initial_clustering.append(current_cluster)
#     self.construction_log[tuple([protein for protein in current_cluster])] = current_cluster_construction_log
#     index += 1
#     if not self.seed_from_all:
#         considered_vertices.add(current_seed)
#         for v in current_cluster:
#             considered_vertices.add(v)
#     print("CLUSTER #%s: %s" % (str(len(self.initial_clustering)), str([vertex for vertex in current_cluster])))