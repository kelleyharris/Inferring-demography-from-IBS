### Input: two VCF files with data mapped to the same reference genome with different sets of variable sites 
### Output: two vcf files with the same number of lines corresponding to the same set of SNPs

import sys

pop1, pop2 = sys.argv[1], sys.argv[2]

if len(sys.argv)>3:
    label1, label2 = sys.argv[3], sys.argv[4]
else:
    label1, label2 = pop1, pop2

race1_infile=open(pop1+'.popdata')
race2_infile=open(pop2+'.popdata')

race1_lines=race1_infile.readlines()
race2_lines=race2_infile.readlines()

race1_infile.close()
race2_infile.close()

race1_split=race1_lines[0].split('\t')
race2_split=race2_lines[0].split('\t')

race1_numseqs=len(race1_split[3])
race2_numseqs=len(race2_split[3])

index1=0
index2=0

race1_output=''
race2_output=''

while index1<len(race1_lines) or index2<len(race2_lines):
    split1=race1_lines[index1].split('\t')
    split2=race2_lines[index2].split('\t')
    if split1[1]==split2[1] and split1[0]==split2[0]:
        race1_output+=race1_lines[index1]
        race2_output+=race2_lines[index2]
        index1+=1
        index2+=1
    elif int(split1[0])<int(split2[0]) or (split1[0]==split2[0] and int(split1[1])<int(split2[1])):
        race2_output+=split1[0]+'\t'+split1[1]+'\t'+split1[2]+'\t'
        for i in range(race2_numseqs):
            race2_output+=split1[2]
        race2_output+='\t\n'
        race1_output+=race1_lines[index1]
        index1+=1
    else:
        race1_output+=split2[0]+'\t'+split2[1]+'\t'+split2[2]+'\t'
        for i in range(race1_numseqs):
            race1_output+=split2[2]
        race1_output+='\t\n'
        race2_output+=race2_lines[index2]
        index2+=1
    while index1==len(race1_lines) and index2<len(race2_lines):
        split2=race2_lines[index2].split('\t')
        race1_output+=split2[0]+'\t'+split2[1]+'\t'+split2[2]+'\t'
        for i in range(race1_numseqs):
            race1_output+=split2[2]
        race1_output+='\t\n'
        race2_output+=race2_lines[index2]
        index2+=1

    while index2==len(race2_lines) and index1<len(race1_lines):
        split1=race1_lines[index1].split('\t')
        race2_output+=split1[0]+'\t'+split1[1]+'\t'+split1[2]+'\t'
        for i in range(race2_numseqs):
            race2_output+=split1[2]
        race2_output+='\t\n'
        race1_output+=race1_lines[index1]
        index1+=1

race1_outfile=open(label1+'_for_'+label2+'.popdata','w')
race2_outfile=open(label2+'_for_'+label1+'.popdata','w')

race1_outfile.write(race1_output)
race2_outfile.write(race2_output)
            
race1_outfile.close()
race2_outfile.close()
