
def n_w(fname_prefix):
    old = open(fname_prefix+".txt", "r")
    new = open(fname_prefix+"+naively_weighted.txt", 'w')
    for line in old.readlines():
        print(repr(line))
        to_write = line.replace("\n","\t1.0\n")
        new.write(to_write)




n_w('biogrid_yeast_physical_unweighted')