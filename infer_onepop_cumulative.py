from parse_tract_file import get_spectrum_counts
import math
import scipy.optimize
import random
from math import exp

import sys
pop1_filename=sys.argv[1]
min_tract_length=int(sys.argv[2])
num_size_changes=len(sys.argv)-3

from demographic_function_builder import *

def tract_lengths(L,t_diff_vec, N_vec): # returns a history with len(t_diff_vec) population size changes
    uncoal = 1
    prob_list = []
    prob_list = initialize_pop(L,prob_list, uncoal, N_vec[0])
    time=0
    for i in range(num_size_changes):
        time+=t_diff_vec[i]
        uncoal = time_lapse(uncoal, N_vec[i], t_diff_vec[i])
        prob_list, uncoal = popsize_change(L, prob_list, uncoal, N_vec[i], N_vec[i+1],time)
    return math.fsum(prob_list)

pop1_lengths, pop1_non_cumul, pop1_cumul, pop1_total_length = get_spectrum_counts(open(pop1_filename))

L_series=[1]
for a in range(1,int(math.log(10**7)/math.log(2))):
    if int(1*2**a)>L_series[0]:
        L_series.insert(0,int(1*2**a))

min_tract_index=int(math.log(min_tract_length)/math.log(2))
print min_tract_index, len(L_series)

L_series_index=0
pop1_lengths_binned=[0]


for i in range(1, len(pop1_lengths)+1):
    if pop1_lengths[-i]>=L_series[L_series_index]:
        pop1_lengths_binned[-1]+=pop1_non_cumul[-i]
    else:
        while pop1_lengths[-i]<L_series[L_series_index]:
            L_series_index+=1
            pop1_lengths_binned.append(0)
        pop1_lengths_binned[-1]=pop1_cumul[-i]

def neg_log_poisson(k,lamb):
    if k==0:
        return lamb
    else:
        return lamb-k*math.log(lamb)+0.5*math.log(2*math.pi*k)+k*math.log(k)-k

def likelihood(x0):
    like=0
    last_expected=0
    for i in range(len(L_series)-min_tract_index):
        L=L_series[i]

        expected = tract_lengths(L,x0[:num_size_changes],x0[num_size_changes:])*pop1_total_length
        like+=neg_log_poisson(pop1_lengths_binned[i],expected)
#        last_expected = expected
    return like

bds = []
for i in range(num_size_changes):
    bds.append([0.001,float(sys.argv[i+3])])
for i in range(num_size_changes+1):
    bds.append([0.05,10])

x0=[]

for i in range(len(bds)):
    x0.append(random.uniform(bds[i][0],bds[i][1]))

print x0


print scipy.optimize.fmin_l_bfgs_b(likelihood, x0, approx_grad=True, bounds=bds, factr=10000)

    
