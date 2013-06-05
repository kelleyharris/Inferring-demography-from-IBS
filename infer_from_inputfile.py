from parse_tract_file import get_spectrum_counts
import math
import scipy.optimize
import random
from math import exp

import sys
pop1_filename=sys.argv[1]
pop2_filename=sys.argv[2]
pop1_pop2_filename=sys.argv[3]
min_tract_length=int(sys.argv[4])

from demographic_inputfile import *

pop1_lengths, pop1_non_cumul, pop1_cumul, pop1_total_length = get_spectrum_counts(open(pop1_filename))
pop2_lengths, pop2_non_cumul, pop2_cumul, pop2_total_length = get_spectrum_counts(open(pop2_filename))
pop1_pop2_lengths, pop1_pop2_non_cumul, pop1_pop2_cumul, pop1_pop2_total_length = get_spectrum_counts(open(pop1_pop2_filename))

print 'total lengths: pop1, pop2, between = ',pop1_total_length, pop2_total_length, pop1_pop2_total_length

L_series=[1]
for a in range(1,18):
    if int(1*2**a)>L_series[0]:
        L_series.insert(0,int(1*2**a))

min_tract_index=int(math.log(min_tract_length)/math.log(2))

L_series_index=0
pop1_lengths_binned=[0]
pop1_pop2_lengths_binned=[0]
pop2_lengths_binned=[0]


for i in range(1, len(pop1_lengths)+1):
    if pop1_lengths[-i]>=L_series[L_series_index]:
        pop1_lengths_binned[-1]+=pop1_non_cumul[-i]
    else:
        while pop1_lengths[-i]<L_series[L_series_index]:
            L_series_index+=1
            pop1_lengths_binned.append(0)
        pop1_lengths_binned[-1]+=pop1_non_cumul[-i]

L_series_index=0

for i in range(1, len(pop1_pop2_lengths)+1):
    if pop1_pop2_lengths[-i]>=L_series[L_series_index]:
        pop1_pop2_lengths_binned[-1]+=pop1_pop2_non_cumul[-i]
    else:
        while pop1_pop2_lengths[-i]<L_series[L_series_index]:
            L_series_index+=1
            pop1_pop2_lengths_binned.append(0)
        pop1_pop2_lengths_binned[-1]+=pop1_pop2_non_cumul[-i]

L_series_index=0

for i in range(1, len(pop2_lengths)+1):
    if pop2_lengths[-i]>=L_series[L_series_index]:
        pop2_lengths_binned[-1]+=pop2_non_cumul[-i]
    else:
        while pop2_lengths[-i]<L_series[L_series_index]:
            L_series_index+=1
            pop2_lengths_binned.append(0)
        pop2_lengths_binned[-1]+=pop2_non_cumul[-i]

def neg_log_poisson(k,lamb):
    if k==0:
        return lamb
    else:
        return lamb-k*math.log(lamb)+0.5*math.log(2*math.pi*k)+k*math.log(k)-k

def likelihood(x0):
    like=0
    last_expected_pop1=0
    last_expected_pop1_pop2=0
    last_expected_pop2=0
    for i in range(len(L_series)-min_tract_index):
        L=L_series[i]

        expected_pop1 = tract_abundance(L,(1,0,0),x0)*pop1_total_length
        like+=neg_log_poisson(pop1_lengths_binned[i],expected_pop1-last_expected_pop1)
        last_expected_pop1=expected_pop1
        
        expected_pop2 = tract_abundance(L,(0,1,0),x0)*pop2_total_length
        like+=neg_log_poisson(pop2_lengths_binned[i],expected_pop2-last_expected_pop2)
        last_expected_pop2=expected_pop2
        
        expected_pop1_pop2 = tract_abundance(L,(0,0,1),x0)*pop1_pop2_total_length
        like+=neg_log_poisson(pop1_pop2_lengths_binned[i],expected_pop1_pop2-last_expected_pop1_pop2)
        last_expected_pop1_pop2=expected_pop1_pop2
    return like

x0=[]

bds=optimization_bounds()

for i in range(len(bds)):
    x0.append(random.uniform(bds[i][0], bds[i][1]))

print x0


print scipy.optimize.fmin_l_bfgs_b(likelihood, x0, approx_grad=True, bounds=bds, factr=10000)

    
