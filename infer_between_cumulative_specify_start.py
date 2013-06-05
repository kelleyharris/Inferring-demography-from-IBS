from parse_tract_file import get_spectrum_counts
import math
import scipy.optimize
import random
from math import exp

import sys
pop1_filename=sys.argv[1]
min_tract_length=int(sys.argv[2])
#min_tract_index=int(sys.argv[2])
num_size_changes=(len(sys.argv)-2)/3

### Change this parameter to a value higher than 1 if your data show evidence of variable mutation rate (see manual). var_mutrate_factor=39 is recommended for human data.
var_mutrate_factor=1
#min_tract_length=5
###

from calc_ibs_backcoal_varmu import prob_L_from_mut_precise_varmu
def prob_L_from_mut_precise(L,ts,N):
    return prob_L_from_mut_precise_varmu(L,ts,ts*2,1,N)


from demographic_function_builder import *

def tract_lengths(L,t_diff_vec, N_vec): # returns a history with len(t_diff_vec) population size changes
    uncoal = 1
    prob_list = [prob_L_from_mut_precise(L,t_diff_vec[0],N_vec[0])]
    time = t_diff_vec[0]
    for i in range(1,num_size_changes):
        time += t_diff_vec[i]
        uncoal = time_lapse(uncoal, N_vec[i-1], t_diff_vec[i])
        prob_list, uncoal = popsize_change(L, prob_list, uncoal, N_vec[i-1], N_vec[i],time)
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
for i in range(num_size_changes):
    bds.append([0.05,10])

x0=[]

for i in range(len(bds)):
    x0.append(float(sys.argv[i+3+num_size_changes])*random.uniform(0.9,1.1))

"""
x0=[0.01730113, 0.04280518, 1.0, 0.9445281, 1.83794258, 0.20138663, 0.89519478, 1.50150796, 1.09407398]

for i in range(len(bds)):
    x0[i]*=random.uniform(0.95,1.05)
"""

print x0


print scipy.optimize.fmin_l_bfgs_b(likelihood, x0, approx_grad=True, bounds=bds, factr=10000)

    
