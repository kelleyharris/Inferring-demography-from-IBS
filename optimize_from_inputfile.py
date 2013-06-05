import sys
import subprocess

african_filename=sys.argv[1]
non_african_filename=sys.argv[2]
between_filename=sys.argv[3]
outfile_name=sys.argv[4]
min_tract_length=sys.argv[5]

output=''

from demographic_inputfile import *

def extract_params_recent(output):
    cut1=output.split('])')
    cut2=cut1[0].split('([')
    if len(cut2)>1:
        s=cut2[1].split(',')

        like_split=cut1[1].split(',')
        like=float(like_split[1])

        params=[]
        for i in range(len(s)):
            params.append(float(s[i]))

        return params, like
    else:
        return [],float('inf')

like_dict1=dict({})
for j in range(20):
    proc1=subprocess.Popen(['nice','python','infer_from_inputfile.py',african_filename,non_african_filename,between_filename, min_tract_length],stdout=subprocess.PIPE)
    out1=proc1.communicate()[0]
    bneck_params, like =extract_params_recent(out1)
    like_dict1[like]=bneck_params
    output+='['
    for param in extract_params_recent(out1):
        output+=str(param)+', '
    output+='], '+str(like)+'\n'

old_params=like_dict1[min(like_dict1.keys())]
output+='['
param_ind=0
bds=optimization_bounds()

for param in old_params:
    output+=str(param)+', '
    if param==bds[param_ind][0]:
        print 'Warning: parameter '+str(param_ind)+' hit lower bound'
    elif param==bds[param_ind][1]:
        print 'Warning: parameter '+str(param_ind)+'  hit upper bound. Consider increasing upper bound and optimizing again.'
    param_ind+=1
output+='], '+str(like)+'\n'
outfile=open(outfile_name,'w')
outfile.write(output)
outfile.close()


