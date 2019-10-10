
def copy_relationship_dictionary(rd):
    return {rel:rd[rel].copy() for rel in rd}


class Basic_Cluster_State:
    def __init__(self, current_cluster, add_candidates, remove_candidates, cohesiveness):
        self.current_cluster = copy_relationship_dictionary(current_cluster)
        self.add_candidates = copy_relationship_dictionary(add_candidates)
        self.remove_candidates = copy_relationship_dictionary(remove_candidates)
        self.cohesiveness = cohesiveness


    def stringify_lite(self):
        s=">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n"
        s+=f"CURRENT_COHESIVENESS: {self.cohesiveness}\n"
        s += "********** CURRENT CLUSTER: **********\n"
        s+= str([protein for protein in self.current_cluster])+"\n"
        s+="********** ADD CANDIDATES : **********\n"
        s+= str([protein for protein in self.add_candidates]) +"\n"
        s+="********** REMOVE CANDIDATES: **********\n"
        s+= str([protein for protein in self.remove_candidates]) +"\n"
        s+=">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n"
        return s


class ClusterState(Basic_Cluster_State):
    __output_level = "verbose"

    def __init__(self,
                 current_cluster,
                 add_candidates,
                 remove_candidates,
                 cohesiveness,
                 best_change=None,
                 best_change_score=None,
                 current_cluster_weight_in=None,
                 current_cluster_weight_out=None,
                 neighborhood_2=None,
                 neighborhood_3=None,
                 last_failed_add_round_no = None,
                 last_failed_remove_round_no=None,
                 round_no=None,
                 number_of_shakes=None):
        self.neighborhood_2 = neighborhood_2
        self.neighborhood_3 = neighborhood_3

        self.best_change=best_change
        self.best_change_score=best_change_score
        self.current_cluster_weight_in=current_cluster_weight_in
        self.current_cluster_weight_out=current_cluster_weight_out

        self.last_failed_add_round_no = last_failed_add_round_no
        self.last_failed_remove_round_no = last_failed_remove_round_no
        self.round_no = round_no
        self.local_number_of_shakes_remaining = number_of_shakes

        super(ClusterState,self).__init__(current_cluster=current_cluster,
                                          add_candidates=add_candidates,
                                          remove_candidates=remove_candidates,
                                          cohesiveness=cohesiveness)

    def make_backup(self):
        return ClusterState(self.current_cluster, self.add_candidates, self.remove_candidates, self.cohesiveness)

    def stringify_heavy(self, graph):
        s=">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n"
        s+=f"CURRENT_COHESIVENESS: {self.cohesiveness}\n"
        s += "\n"
        s+="**************************************\n"
        s += "********** CURRENT CLUSTER: **********\n"
        s+="**************************************\n"
        s += "\n"
        for protein in self.current_cluster:
            s+=str(protein)+", "+graph.id_to_name[protein]+"\n"
            s+=self.current_cluster[protein].stringify()
        s += "\n"
        s += "*************************************\n"
        s+="********** ADD CANDIDATES: **********\n"
        s+="*************************************\n"
        s += "\n"
        for protein in self.add_candidates:
            s+=str(protein)+", "+graph.id_to_name[protein]+"\n"
            s+=self.add_candidates[protein].stringify()
        s+="\n"
        s+="****************************************\n"
        s+="********** REMOVE CANDIDATES: **********\n"
        s+="****************************************\n"
        s += "\n"
        for protein in self.remove_candidates:
            s+=str(protein)+", "+graph.id_to_name[protein]+"\n"
            s+=self.remove_candidates[protein].stringify()
        s+=">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n"
        return s

    def __str__(self):
        if self.__output_level =="verbose":
            return self.stringify_heavy(self)
        else:
            return self.stringify_lite(self)
