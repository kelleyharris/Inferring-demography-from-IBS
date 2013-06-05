##### 
standard_popsize = 10000 ## Your mutation and recombination rates should be scaled with respect to this population size 
#####

## This parameter specifies that we will look at IBS tracts longer than min_tract_length base pairs

from calc_ibs_backcoal_varmu import prob_L_from_mut_precise_varmu
from parse_tract_file import get_spectrum_counts
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from math import exp
from math import log

import pylab
from matplotlib import rc


from demographic_inputfile import *

import sys
non_african_tractfile=sys.argv[1]
african_tractfile=sys.argv[2]
between_tractfile=sys.argv[3]
result_file=sys.argv[4]
output_file=sys.argv[5]
gen_time=int(sys.argv[6])
min_tract_length = int(sys.argv[7])

from demographic_inputfile import *
fig,axs = plt.subplots(nrows=2, ncols=2, sharex=True,sharey=True)
ax0,ax1,ax2,ax3=axs[0,0],axs[0,1],axs[1,0],axs[1,1]


YRI_lengths, YRI_non_cumul, YRI_cumul, YRI_total_length = get_spectrum_counts(open(non_african_tractfile))

L_series=[1]
for a in range(1,int(log(5*20**6)/log(1.25))):
    if int(1*1.25**a)>L_series[0]:
        L_series.insert(0,int(1*1.25**a))
#print L_series
total_length=YRI_total_length

#print 'total CEU length: ',YRI_total_length

last_prob_0, last_prob_1, last_prob_2, last_prob_3, last_prob_4, last_prob_5 = 0,0,0,0,0,0
y0, y1, y2, y3,y4,y5 = [],[],[], [],[],[]

results=open(result_file)
lines=results.readlines()
results.close()

while len(lines[-1])<2:
    lines.pop()
s1=lines[-1].lstrip('[').split(']')
s2=s1[0].split(', ')

params=[]
for entry in s2:
    if len(entry)>1:
        params.append(float(entry))

print describe_history(params, gen_time, standard_popsize)

for i in range(len(L_series)-1):
    L=L_series[i]
    prob5=tract_abundance(L,(1,0,0),params)
    y5.append(-(prob5-last_prob_5)/(L_series[i]-L_series[i-1]))
    last_prob_5=prob5


L_series_index=0
YRI_lengths_binned=[0]
for i in range(1, len(YRI_lengths)+1):
    if YRI_lengths[-i]>=L_series[L_series_index]:
        YRI_lengths_binned[-1]+=YRI_non_cumul[-i]
    else:
        while YRI_lengths[-i]<L_series[L_series_index]:
            L_series_index+=1
            YRI_lengths_binned.append(0)
        YRI_lengths_binned[-1]+=YRI_non_cumul[-i]


def bin_dataset(data_file):
    data_lengths, data_non_cumul, data_cumul, total_length= get_spectrum_counts(data_file)
    L_series_index=0
    data_lengths_binned=[0]
    for i in range(1, len(data_lengths)+1):
        if data_lengths[-i]>=L_series[L_series_index]:
            data_lengths_binned[-1]+=data_non_cumul[-i]
        else:
            while data_lengths[-i]<L_series[L_series_index]:
                L_series_index+=1
                data_lengths_binned.append(0)
            data_lengths_binned[-1]+=data_non_cumul[-i]
    for i in range(1,len(L_series)):
        data_lengths_binned[i]*=1.0/(total_length*(-L_series[i]+L_series[i-1]))
    data_lengths_binned[0]*=1.0/total_length
    return data_lengths_binned, total_length

YRI_lengths_binned,total_length=bin_dataset(open(non_african_tractfile))
ax1.loglog(L_series, YRI_lengths_binned, label=non_african_tractfile)
CEU_lengths_binned=YRI_lengths_binned
ax1.loglog(L_series[:-1], y5, label='Inferred history')
ax1.legend(loc='lower left')

y5=[]
last_prob_5=0
for i in range(len(L_series)-1):
    L=L_series[i]
    prob5=tract_abundance(L, (0,1,0),params)
    y5.append(-(prob5-last_prob_5)/(L_series[i]-L_series[i-1]))
    last_prob_5=prob5

YRI_lengths_binned,total_length=bin_dataset(open(african_tractfile))
real_YRI_lengths_binned=YRI_lengths_binned
ax2.loglog(L_series, YRI_lengths_binned, label=african_tractfile)
ax2.loglog(L_series[:-1], y5, label='Inferred history')
ax2.legend(loc='lower left')

y5=[]
last_prob_5=0
for i in range(len(L_series)-1):
    L=L_series[i]
    prob5 = tract_abundance(L, (0,0,1), params)
    y5.append(-(prob5-last_prob_5)/(L_series[i]-L_series[i-1]))
    last_prob_5=prob5

YRI_lengths_binned,total_length=bin_dataset(open(between_tractfile))
CEU_YRI_lengths_binned=YRI_lengths_binned
ax3.loglog(L_series, YRI_lengths_binned, label=between_tractfile)
ax3.loglog(L_series[:-1], y5, label='Inferred history')
ax3.legend(loc='lower left')

ax0.loglog(L_series,CEU_YRI_lengths_binned, label=between_tractfile)
ax0.loglog(L_series,CEU_lengths_binned, label=non_african_tractfile)
ax0.loglog(L_series,real_YRI_lengths_binned, label=african_tractfile)
ax0.legend(loc='lower left')

plt.xlabel('IBS tract length (L)')

plt.ylim(ymin=10**(-16))
plt.xlim(xmin=min_tract_length)
plt.xlim(xmax=10**6)
plt.xscale('log')
plt.yscale('log')

F=pylab.gcf()
F.set_size_inches(6.83,7)
rc('font',**{'family':'sans-serif','sans-serif':['Arial'],'size':5})
plt.savefig(output_file,dpi=300, format='pdf')
    
              
    
