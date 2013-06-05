import sys

infile=open(sys.argv[1])
seq_length=int(sys.argv[2])
lines=infile.readlines()
infile.close()

for line in lines:
    if line.startswith('positions:'):
        line.strip('\n')
        seg_sites=line.split(' ')
        for i in range(2,len(seg_sites)-1):
            length=seq_length*(float(seg_sites[i])-float(seg_sites[i-1]))
            length_out+=length
            if length>=1:
                lengths.append(int(length))
lengths.sort()
outfile=open(sys.argv[3],'w')

i=1 
count=0
curr_length=1
for l in lengths:
    if l==curr_length:
        count+=1
    else:
        outfile.write(str(count)+'\t'+str(curr_length)+'\n')
        curr_length=l
        count=1
outfile.close()
