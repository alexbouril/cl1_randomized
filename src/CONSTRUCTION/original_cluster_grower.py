from src.COMMON.cmn import *
from src.QUALITY.quality import cohesiveness
from src.CONSTRUCTION.find import find_best_add, find_best_remove
from src.CONSTRUCTION.modify import add,remove

def ocg(cl1, cs, current_cluster_construction_log):

    def take_step(cs, change_function, current_cluster_construction_log, change_type):
        cs.cohesiveness = cs.best_change_score
        change_function(cl1, cs)
        current_cluster_construction_log.append(Action(change_type, cs.best_change))
        current_cluster_construction_log.append(cs.make_backup())
        ##############################################
        # SANITY CHECK
        ##############################################
        # a = cs.cohesiveness
        # b = cohesiveness(cl1, [v for v in cs.current_cluster])
        # offby=abs(a-b )
        # if offby >.005:
        #     print("math problem", a, b, offby)

    def figure_out_which_move_is_best(c1, cs):
        ############
        # reset
        ############
        cs.best_change=None
        cs.best_change_score=None
        ############
        # find the best add step
        ############
        best_add, best_add_score = find_best_add(cl1, cs)
        ############
        # reset
        ############
        cs.best_change=None
        cs.best_change_score=None
        ############
        # find the best remove step
        ############
        best_remove, best_remove_score = find_best_remove(cl1, cs)

        if not best_add and not best_remove:
            return None
        elif best_add and not best_remove:
            cs.best_change = best_add
            cs.best_change_score = best_add_score
            return 1
        elif best_remove and not best_add:
            cs.best_change = best_remove
            cs.best_change_score = best_remove_score
            return 2
        elif best_add_score > best_remove_score:
            cs.best_change = best_add
            cs.best_change_score = best_add_score
            return 1
        else:
            cs.best_change = best_remove
            cs.best_change_score = best_remove_score
            return 2


    while True:
        choice = figure_out_which_move_is_best(cl1, cs)
        if choice is None:
            break
        ############################################################
        # CONSIDER ADDING A VERTEX ON THE BOUNDARY
        ############################################################
        if choice == 1:
            take_step(cs, add, current_cluster_construction_log,"adding")
        ############################################################
        # CONSIDER REMOVING A VERTEX ON THE BOUNDARY
        ############################################################
        else:
            take_step(cs, remove, current_cluster_construction_log, "removing")
    ##############################################
    # add current_cluster to the list of clusters
    #############################################
    cl1.construction_log[tuple([protein for protein in cs.current_cluster])] = current_cluster_construction_log
    print("CLUSTER #%s: %s" % (str(len(cl1.initial_clustering)), str([vertex for vertex in cs.current_cluster])))
    print("COHESIVENESS: ", cs.cohesiveness)
    print("LENGTH: ", len(cs.current_cluster))

    return cs.current_cluster.copy()

