from parse_tract_file import get_spectrum_counts
from math import log
import matplotlib.pyplot as plt
import sys

import pylab
from matplotlib import rc


L_series=[1]
for a in range(1,int(log(5*20**6)/log(1.25))):
    if int(1*1.25**a)>L_series[0]:
        L_series.insert(0,int(1*1.25**a))
print L_series


def bin_dataset(data_file):
    data_lengths, data_non_cumul, data_cumul, total_length= get_spectrum_counts(data_file)
    print data_lengths[-1]
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
    print 'total length:', total_length
    print 'heterozygosity: ',1.0*data_cumul[0]/total_length
    return data_lengths_binned, total_length

plot_title=sys.argv[-1]
for i in range(1,len(sys.argv)-1):
    lengths_binned, total_length=bin_dataset(open(sys.argv[i]))
    plt.loglog(L_series,lengths_binned,label=sys.argv[i])

plt.legend(loc='lower left')
plt.xlabel('IBS tract length (L)')
plt.ylabel('Frequency of L-base IBS tracts')
plt.xlim(xmax=10**7)
plt.xlim(xmin=0)

F=pylab.gcf()
F.set_size_inches(4.5,4.5)
rc('font',**{'family':'sans-serif','sans-serif':['Arial'],'size':10})
plt.savefig(plot_title,dpi=300)


              
    
