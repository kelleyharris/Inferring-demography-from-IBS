import sys
import subprocess
import random

from math import exp

tract_filename = sys.argv[1]
min_tract_length = sys.argv[2]
num_reps=int(sys.argv[3])
generation_time = int(sys.argv[4])


def optimize_random_start(upper_bounds):
    command_line=['nice','python','infer_onepop_cumulative.py',tract_filename, min_tract_length]
    for i in range(len(upper_bounds)):
        command_line.append(str(upper_bounds[i]))
    proc=subprocess.Popen(command_line,stdout=subprocess.PIPE)
    output=proc.communicate()[0]
    s=output.split('])')
    s1=s[0].split('([')
    s2=s1[1].split(',')
    for i in range(len(s2)):
        s2[i]=float(s2[i])
    return s2

def randomize_and_optimize(upper_bounds, x0):
    for i in range(len(x0)):
        x0[i]*=random.uniform(0.95,1.05)
    command_line=['nice','python','infer_onepop_cumulative_specify_start.py',tract_filename, min_tract_length]
    for i in range(len(upper_bounds)):
        command_line.append(str(upper_bounds[i]))
    for i in range(len(x0)):
        command_line.append(str(x0[i]))
    proc=subprocess.Popen(command_line,stdout=subprocess.PIPE)
    output=proc.communicate()[0]
    s=output.split('])')
    print s
    s1=s[0].split('([')
    return s1
    

def optimize_prescribed_start(upper_bounds, x0):
    s1=[]
    while len(s1)<2: 
        s1=randomize_and_optimize(upper_bounds,x0)
    s2=s1[1].split(',')
    print s2
    for i in range(len(s2)):
        s2[i]=float(s2[i])
    return s2

def choose_next_params(upper_bounds, x0):
    hitbounds=False
    new_upper_bounds=[]
    for i in range(len(upper_bounds)):
        if abs(upper_bounds[i]-x0[i])<0.001:
            hitbounds=True
            new_upper_bounds.append(1.5*upper_bounds[i])
        else:
            new_upper_bounds.append(upper_bounds[i])
    if hitbounds:
        return new_upper_bounds, x0
    else:
        num_size_changes=len(upper_bounds)
        prob_coal=[1-exp(-x0[0]/x0[num_size_changes])]
        margin=1-prob_coal[0]
        index_max_prob=0
        for i in range(1,num_size_changes):
            prob_no_coal=exp(-x0[i]/x0[num_size_changes+i])
            prob_coal.append(margin*(1-prob_no_coal))
            margin*=prob_no_coal
            if prob_coal[-1]>prob_coal[index_max_prob]:
                index_max_prob=i
        prob_coal.append(margin)
        if prob_coal[-1]>prob_coal[index_max_prob]:
            index_max_prob=num_size_changes
        if index_max_prob<num_size_changes:
            x0.insert(num_size_changes+index_max_prob, x0[num_size_changes+index_max_prob])
            x0[index_max_prob]*=0.5
            x0.insert(index_max_prob,x0[index_max_prob])
        else:
            x0.append(x0[-1])
            x0.insert(num_size_changes,x0[num_size_changes-1])
        new_upper_bounds=[]
        for i in range(num_size_changes+1):
            new_upper_bounds.append(2*x0[i])
        return new_upper_bounds, x0

def parse_history(history):
    output='Population sizes (present to past) :'
    for i in range(len(history)/2, len(history)-1):
        output+=str(10000*history[i])+', '
    output+=str(10000*history[-1])+'\n'
    output+='Population size change times (kya) :'
    time=0
    for i in range(len(history)/2-1):
        time+=history[i]
        output+=str(generation_time*20*time)+', '
    if len(history)/2-1>0:
        time+=history[len(history)/2-1]
        output+=str(generation_time*20*time)
    return output
     

history=optimize_random_start([])
history.append(history[0])
history.insert(0,1.0)
history=optimize_prescribed_start(history,[1.0])
upper_bounds=[1.0]
print ''
print 'Best history so far: '
print parse_history(history)
print ''
print 'Parameters for plotting: '
print history
print ''
print ''

for rep_num in range(num_reps):
    upper_bounds, history = choose_next_params(upper_bounds, history)
#    print upper_bounds, history
    history = optimize_prescribed_start(upper_bounds, history)
    print ''
    print 'Best history so far: '
    print parse_history(history)
    print ''
    print 'Parameters for plotting: '
    print history
    print ''
    print ''
