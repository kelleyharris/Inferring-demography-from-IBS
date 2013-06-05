import sys
import subprocess
import os

input1, input2, output = sys.argv[1],sys.argv[2], sys.argv[3]

subprocess.Popen(['mkdir',output])

chromlist=os.listdir(input1)

for i in range(1,23):

for filename in chromlist:
    outfile=open(output+'/'+filename,'w')
    infile1=open(input1+'/'+filename)
    infile2=open(input2+'/'+filename)
    lines1=infile1.readlines()
    lines2=infile2.readlines()
    infile1.close()
    infile2.close()
    index1=0
    index2=0
    while index1<len(lines1) or index2<len(lines2):
        while  index1<len(lines1) and len(lines1[index1])<2:
            index1+=1
        while index2<len(lines2) and len(lines2[index2])<2:
            index2+=1
        if index1==len(lines1):
            outfile.write(lines2[index2])
            index2+=1
        elif index2==len(lines2):
            outfile.write(lines1[index1])
            index1+=1
        else:
            line1=lines1[index1]
            line2=lines2[index2]
            s1=line1.split('\t')
            start1, end1=int(s1[0]), int(s1[1].strip('\n'))
            s2=line2.split('\t')
            start2, end2= int(s2[0]), int(s2[1].strip('\n'))
            if end1<start2:
                outfile.write(lines1[index1])
                index1+=1
            elif end2<start1:
                outfile.write(lines2[index2])
                index2+=1
            else:
                start=min(start1, start2)
                end=max(end1, end2)
                while (start1<=end and index1<len(lines1)) or (start2<=end and index2<len(lines2)):
                    if start1<=end and index1<len(lines1):
                        end=max(end1, end)
                        index1+=1
                        while index1<len(lines1) and len(lines1[index1])<2:
                            index1+=1
                    else:
                        end=max(end2, end)
                        index2+=1
                        while index2<len(lines2) and len(lines2[index2])<2:
                            index2+=1
                    if index1<len(lines1):
                        line1=lines1[index1]
                        s1=line1.split('\t')
                        start1, end1=int(s1[0]), int(s1[1].strip('\n'))
                    if index2<len(lines2):
                        line2=lines2[index2]
                        s2=line2.split('\t')
                        start2, end2= int(s2[0]), int(s2[1].strip('\n'))
#                    print index1, len(lines1), index2, len(lines2)
                print i, start, end
                outfile.write(str(start)+'\t'+str(end)+'\n')
    outfile.close()
