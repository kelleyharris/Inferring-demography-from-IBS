from demographic_function_builder import *

## The optimization_bounds vector needs to have the same length as the third argument of tract_abundance

def optimization_bounds():
    return [[0.05, 20.0], [0.05, 20.0], [0.05, 40.0], [0.001, 0.5], [0.001, 0.5], [0.001, 0.5]]



## The time arguments of the tract length distribution need to be time differences/ interval lengths rather than ordered absolute times
## e.g. if you want tm < ts, specify the two time arguments tm, ts_diff, where ts = tm+ts_diff

def tract_abundance(L, (uncoalesced_pop1, uncoalesced_pop2, uncoalesced_between), (N1, N2, N, tm, ts_diff, f)):
    prob_list = []
    prob_list = initialize_pop(L, prob_list, uncoalesced_pop1, N1)
    prob_list = initialize_pop(L, prob_list, uncoalesced_pop2, N2)
    uncoalesced_pop1 = time_lapse(uncoalesced_pop1, N1, tm)
    uncoalesced_pop2 = time_lapse(uncoalesced_pop2, N2, tm)
    prob_list, uncoalesced_pop1, uncoalesced_pop2, uncoalesced_between = two_way_admixture(L, prob_list, uncoalesced_pop1, uncoalesced_pop2, uncoalesced_between, N1, N2, tm, f, 0)
    uncoalesced_pop1 = time_lapse(uncoalesced_pop1, N1, ts_diff)
    uncoalesced_pop2 = time_lapse(uncoalesced_pop2, N2, ts_diff)
    ts=tm+ts_diff
    prob_list, uncoalesced_merged = pop_merge(L, prob_list, uncoalesced_pop1, uncoalesced_pop2, uncoalesced_between, N1, N2, N, ts)
    return math.fsum(prob_list)


def describe_history((N1, N2, N, tm, ts_diff, f), gen_time, standard_popsize):
    output = 'Time of most recent gene flow: '+str(tm*gen_time*standard_popsize*0.001)+' kya\n'
    output += 'Divergence time: '+str((tm+ts_diff)*gen_time*standard_popsize*0.001)+' kya\n'
    output += 'Size of population 1: '+str(N1*standard_popsize)+'\n'
    output += 'Size of population 2: '+str(N2*standard_popsize)+'\n'
    output += 'Ancestral population size: '+str(N*standard_popsize)+'\n'
    return output
