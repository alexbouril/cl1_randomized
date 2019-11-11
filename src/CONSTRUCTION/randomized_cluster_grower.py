from src.COMMON.cmn import *
from src.QUALITY.quality import *
from src.CONSTRUCTION.find import find_best_add, find_best_suboptimal_add, careful_find_best_2neighborhood_add,find_best_mixed_measure_add2,find_best_remove, find_best_mixed_measure_add, find_best_density_add
from src.CONSTRUCTION.modify import add,remove

def rcg(cl1, cs, current_cluster_construction_log):
    def update_backup_cluster(cs, best_seen_cs):
        if cs.cohesiveness > best_seen_cs.cohesiveness:
            # TODO make this step faster by updating references, don't need to copy back objects
            best_seen_cs.cohesiveness = cs.cohesiveness
            del best_seen_cs.current_cluster
            del best_seen_cs.add_candidates
            del best_seen_cs.remove_candidates
            best_seen_cs.current_cluster = cs.current_cluster.copy()
            best_seen_cs.add_candidates = cs.add_candidates.copy()
            best_seen_cs.remove_candidates = cs.remove_candidates.copy()

    def try_step(cs, best_seen_cs, search_function, change_function, current_cluster_construction_log, change_type):
        cs.round_no+=1
        search_function(cl1, cs)
        if cs.best_change:
            cs.cohesiveness = cs.best_change_score
            change_function(cl1, cs)
            update_backup_cluster(cs, best_seen_cs)
            current_cluster_construction_log.append(Action(change_type, cs.best_change))
            # current_cluster_construction_log.append(cs.make_backup())
            return True
        else:
            current_cluster_construction_log.append(Action("failed to %s, current cohesiveness: %s" % (change_type, str(cs.cohesiveness))))
            if change_type == "adding":
                cs.last_failed_add_round_no = cs.round_no
            else:
                cs.last_failed_remove_round_no = cs.round_no
            return False
        ##############################################
        # SANITY CHECK
        ##############################################
        # a = cs.cohesiveness
        # b = cohesiveness(cl1, [v for v in cs.current_cluster])
        # offby=abs(a-b )
        # if offby >.005:
        #     print("math problem", a, b, offby)


    ############################################################
    # initialize the backup in a shake heuristic gets us stuck in a lower local optima
    ############################################################
    best_seen_cs = cs.make_backup()
    remove_counter = 0
    number_of_calls_to_find_best_2_neighborhood_add = 0
    tries_2_neighborhood = 0
    max_tries = 2
    while (cs.add_candidates or cs.remove_candidates) and \
            abs(cs.last_failed_remove_round_no - cs.last_failed_add_round_no) != 1 \
            and remove_counter < 100:
        decider = np.random.rand()
        ############################################################
        # CONSIDER ADDING A VERTEX ON THE BOUNDARY
        ############################################################
        if  (decider <= .5 or cs.last_failed_remove_round_no == cs.round_no) and cs.last_failed_add_round_no != cs.round_no:
            if cs.round_no %7==1:
                for counter in range(4):
                    try_step(cs, best_seen_cs,
                             find_best_suboptimal_add,
                             add,
                             current_cluster_construction_log,
                             "adding")
            elif cs.round_no%5==0 and len(cs.current_cluster)>4:
                try_step(cs, best_seen_cs,
                         find_best_mixed_measure_add2,
                         add,
                         current_cluster_construction_log,
                         "adding")
                for counter in range(3):
                    try_step(cs, best_seen_cs,
                             find_best_suboptimal_add,
                             add,
                             current_cluster_construction_log,
                             "adding")

            # elif len(cs.current_cluster) > 5 and cs.round_no % 20 == 1:
            #     # ######################################################################
            #     # # JUST FOR COMPARISON
            #     # #####################################################################
            #     # find_best_add(cl1, cs)
            #     # if cs.best_change is None:
            #     #     print("NO REGULAR ADD AVAILABLE")
            #     # if cs.best_change:
            #     #     print("REGULAR BEST CHANGE AVAILABLE WITH SCORE OF: %s" % str(cs.best_change_score))
            #     ######################################################################
            #     # NOW FOR THE REAL CALL
            #     #####################################################################
            #     number_of_calls_to_find_best_2_neighborhood_add+=1
            #     try_step(cs, best_seen_cs,
            #              careful_find_best_2neighborhood_add,
            #              add,
            #              current_cluster_construction_log,
            #              "adding")
            #     if cs.best_change:
            #         ######################################################################
            #         # TRY TO ADD SUBOPTIMALLY A FEW TIMES
            #         #####################################################################
            #         for counter in range(4):
            #             try_step(cs, best_seen_cs,
            #                     find_best_suboptimal_add,
            #                     add,
            #                     current_cluster_construction_log,
            #                     "adding")
            #         ######################################################################
            #         # TRY TO ADD REGULARLY A FEW TIMES
            #         #####################################################################
            #         for counter in range(2):
            #             try_step(cs, best_seen_cs,
            #                     find_best_add,
            #                     add,
            #                     current_cluster_construction_log,
            #                     "adding")
            else:
                try_step(cs, best_seen_cs,
                         find_best_add,
                         add,
                         current_cluster_construction_log,
                         "adding")
        ############################################################
        # CONSIDER REMOVING A VERTEX ON THE BOUNDARY
        ############################################################
        if (decider > .5 or cs.last_failed_add_round_no == cs.round_no) and cs.last_failed_remove_round_no != cs.round_no:
            try_step(cs, best_seen_cs,
                     find_best_remove,
                     remove,
                     current_cluster_construction_log,
                     "removing")
        # if tries_2_neighborhood < max_tries and  abs(cs.last_failed_remove_round_no - cs.last_failed_add_round_no) == 1:
        #     tries_2_neighborhood+=1
        #     # ######################################################################
        #     # # JUST FOR COMPARISON
        #     # #####################################################################
        #     # find_best_add(cl1, cs)
        #     # if cs.best_change is None:
        #     #     print("NO REGULAR ADD AVAILABLE")
        #     # if cs.best_change:
        #     #     print("REGULAR BEST CHANGE AVAILABLE WITH SCORE OF: %s" % str(cs.best_change_score))
        #     ######################################################################
        #     # NOW FOR THE REAL CALL
        #     #####################################################################
        #     number_of_calls_to_find_best_2_neighborhood_add+=1
        #     try_step(cs, best_seen_cs,
        #              careful_find_best_2neighborhood_add,
        #              add,
        #              current_cluster_construction_log,
        #              "adding")
        #     if cs.best_change:
        #         ######################################################################
        #         # TRY TO ADD SUBOPTIMALLY A FEW TIMES
        #         #####################################################################
        #         for counter in range(4):
        #             try_step(cs, best_seen_cs,
        #                     find_best_suboptimal_add,
        #                     add,
        #                     current_cluster_construction_log,
        #                     "adding")
        #         ######################################################################
        #         # TRY TO ADD REGULARLY A FEW TIMES
        #         #####################################################################
        #         for counter in range(2):
        #             try_step(cs, best_seen_cs,
        #                     find_best_add,
        #                     add,
        #                     current_cluster_construction_log,
        #                     "adding")
        # # TODO: make Bad Step work with the new object schema
    #########################################
    # in the case that the current_cluster is not the best one that we saw, revert to the best one that we saw
    ##########################################
    if best_seen_cs.cohesiveness > cs.cohesiveness:
        current_cluster_construction_log.append(Action("reverting to previous state"))
        current_cluster_construction_log.append(best_seen_cs)
    ##############################################
    # add current_cluster to the list of clusters
    #############################################
    cl1.construction_log[tuple([protein for protein in best_seen_cs.current_cluster])] = current_cluster_construction_log
    print("best_seen", best_seen_cs.cohesiveness,[v for v in best_seen_cs.current_cluster])
    print("CLUSTER #%s: %s" % (str(len(cl1.initial_clustering)), str([vertex for vertex in best_seen_cs.current_cluster])))
    print("COHESIVENESS: ", best_seen_cs.cohesiveness)
    print("LENGTH: ", len(cs.current_cluster))
    print("number_of_calls_to_find_best_2_neighborhood_add: ",number_of_calls_to_find_best_2_neighborhood_add)
    return best_seen_cs.current_cluster