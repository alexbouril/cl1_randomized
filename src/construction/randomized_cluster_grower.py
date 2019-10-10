from  src.cl1_randomized.cl1_randomized import *
from src.common.common import *
from src.construction.construction_operations import *

def randomized_cluster_grower(cl1: CL1_Randomized,cs:ClusterState, current_cluster_construction_log):
    def update_backup_cluster(cs, best_seen_cs):
        if cs.cohesiveness > best_seen_cs.cohesiveness:
            best_seen_cs = cs.make_backup()

    def try_step(cs, best_seen_cs, search_function, change_function, current_cluster_construction_log, change_type):
        print(search_function)
        print(change_function)
        cs.round_no+=1
        search_function(cl1, cs)
        if cs.best_change:
            print("yeah")
            change_function(cl1, cs)
            update_backup_cluster(cs, best_seen_cs)
            current_cluster_construction_log.append(Action(change_type, cs.best_change))
            current_cluster_construction_log.append(cs.make_backup())
        else:
            print("nope")
            current_cluster_construction_log.append(Action("failed to %s, current cohesiveness: %s" % (change_type, str(cs.cohesiveness))))
            if change_type == "adding":
                cs.last_failed_add_round_no = cs.round_no
            else:
                cs.last_failed_remove_round_no = cs.round_no

    ############################################################
    # initialize the backup in a shake heuristic gets us stuck in a lower local optima
    ############################################################
    best_seen_cs = cs.make_backup()
    remove_counter = 0
    while (cs.add_candidates or cs.remove_candidates) and \
            abs(cs.last_failed_remove_round_no - cs.last_failed_add_round_no) != 1 \
            and remove_counter < 100:
        decider = numpy.random.rand()
        ############################################################
        # CONSIDER ADDING A VERTEX ON THE BOUNDARY
        ############################################################
        if (decider <= .5 or cs.last_failed_remove_round_no == cs.round_no) and cs.last_failed_add_round_no != cs.round_no:
            if len(cs.current_cluster) > 5 and cs.round_no % 4 == 1:
                cs.round_no += 1
                # ######################################################################
                # # JUST FOR COMPARISON
                # #####################################################################
                # find_best_add(cl1, cs)
                # if cs.best_change is None:
                #     print("NO REGULAR ADD AVAILABLE")
                # if cs.best_change:
                #     print("REGULAR BEST CHANGE AVAILABLE WITH SCORE OF: %s" % str(cs.best_change_score))
                ######################################################################
                # NOW FOR THE REAL CALL
                #####################################################################
                try_step(cs,
                          best_seen_cs,
                          careful_find_best_2neighborhood_add(),
                          add,
                          current_cluster_construction_log,
                          "adding")
                if cs.best_change:
                    for counter in range(4):
                        try_step(cs, best_seen_cs,
                                find_best_suboptimal_add,
                                add,
                                current_cluster_construction_log,
                                "adding")
                    for counter in range(2):
                        try_step(cs, best_seen_cs,
                                find_best_add,
                                add,
                                current_cluster_construction_log,
                                "adding")
            else:
                try_step(cs,
                          best_seen_cs,
                          find_best_add,
                          add,
                          current_cluster_construction_log,
                          "adding")
        ############################################################
        # CONSIDER REMOVING A VERTEX ON THE BOUNDARY
        ############################################################
        if (decider > .5 or cs.last_failed_add_round_no == cs.round_no) and cs.last_failed_remove_round_no != cs.round_no:
            try_step(cs,
                      best_seen_cs,
                      find_best_add,
                      remove,
                      current_cluster_construction_log,
                      "removing")

        # TODO: make Bad Step work with the new object schema

    #########################################
    # in the case that the current_cluster is not the best one that we saw, revert to the best one that we saw
    ##########################################
    if best_seen_cs.cohesiveness > cs.cohesiveness:
        current_cluster_construction_log.append(Action("reverting to previous state"))
        current_cluster_construction_log.append(best_seen_cs)

    # add current_cluster to the list of clusters
    cl1.initial_clustering.append(best_seen_cs.current_cluster)
    cl1.construction_log[tuple([protein for protein in best_seen_cs.current_cluster])] = current_cluster_construction_log
    print("CLUSTER #%s: %s" % (str(len(cl1.initial_clustering)), str([vertex for vertex in best_seen_cs.current_cluster])))