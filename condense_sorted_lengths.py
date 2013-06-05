import sys

filename=sys.argv[1]

infile=open(filename+'_sorted.txt')
lines =infile.readlines()
infile.close()
i=0
while not lines[i][0] in '0123456789':
    i+=1

curr_length=1
count=0
outfile=open(filename+'.ibs', 'w')
while i<len(lines):
    length=int(float(lines[i].strip('\n')))
    if length==curr_length:
        count+=1
    else:
        outfile.write(str(count)+'\t'+str(curr_length)+'\n')
        count=1
        curr_length=length
    i+=1
outfile.close()
