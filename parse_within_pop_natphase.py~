import sys
import subprocess

## Required arguments: race1, race2, maskFolder, percentMissing (if not using a maskFolder, enter maskFolder='None')
## Optional arguments: hapstart1, hapend1, hapstart2, hapend2 (need to enter all these arguments or none of them)

## If the input files are not in the same directory as this script, enter the path to their directory here.
filepath=''

race1=sys.argv[1]

maskFolder, percentMissing = sys.argv[2], float(sys.argv[3])

if maskFolder=='None':
    masked=False
    gap_start, gap_end, next_gap_start, next_gap_end = float('inf'), float('inf'), float('inf'), float('inf')
else:
    masked=True


infile_race1=open(filepath+race1+'.popdata')
race1_lines=infile_race1.readlines()
race2_lines=race1_lines

split1=race1_lines[0].split('\t')
split2=race2_lines[0].split('\t')

num_haps1=len(split1[3])
num_haps2=num_haps1

if len(sys.argv)>4:
    hapstart1, hapend1 = int(sys.argv[4]), int(sys.argv[5])
    hapstart2, hapend2 = hapstart1, hapend1
else:
    hapstart1, hapend1, hapstart2, hapend2 = 0, num_haps1, 0, num_haps2

lengths=[]

length_locations=[]

def update_gap(gaps,gap_index):
    if gap_index>=len(gaps):
        return float('inf'), float('inf')
    else:
        while len(gaps[gap_index])<2:
            gap_index+=1
        gsplit=gaps[gap_index].split('\t')
        gap_start=int(gsplit[0])
        gap_end=int(gsplit[1])
        return gap_start,gap_end

for i_over_2 in range(hapstart1/2, hapend1/2):
    i=2*i_over_2
    for j in range(i+1, i+2):
        last_chrom=0
        line_index=0
        last_pos=0
        count_missing=0
        while line_index < len(race1_lines):
            line1=race1_lines[line_index]
            line2=race2_lines[line_index]
            lsplit1=line1.split('\t')
            lsplit2=line2.split('\t')
            chrom=int(lsplit1[0])
            alleles1=lsplit1[3]
            alleles2=lsplit2[3]
            if alleles1[i]=='N' or alleles2[j]=='N':
                count_missing+=1
            position=int(lsplit1[1])
            if not chrom==last_chrom:
                if masked:
                    gapfile=open(maskFolder+'/chrom'+str(chrom)+'.txt')
                    gaps=gapfile.readlines()
                    gap_index=0
                    gap_start, gap_end=update_gap(gaps,gap_index)
                    skipped_gap=True
                last_pos=0
                count_missing=0
                last_chrom=chrom
#            if masked:
#                next_gap_start, next_gap_end=update_gap(gaps,gap_index+1)
            if position<gap_start and position>last_pos:
                if position-last_pos>0 and not (alleles1[i]==alleles2[j]) and not alleles1[i]=='N' and not alleles2[j]=='N':
                    if float(count_missing)/(position-last_pos+1)<percentMissing and not skipped_gap:
                        lengths.append(position-last_pos)
                        length_locations.append(str(chrom)+'\t'+str(last_pos)+'\t'+str(position-last_pos)+'\n')
                    last_pos=position
                    skipped_gap=False
                    count_missing=0
            else:
                while position>=gap_start:
                    last_pos=gap_end
                    count_missing=0
                    gap_index+=1
                    gap_start,gap_end = update_gap(gaps,gap_index)
                    next_gap_start, next_gap_end = update_gap(gaps,gap_index+1)
                    skipped_gap=True
            line_index+=1

output=''
for l in lengths:
    output+=str(l)+'\n'
outfile=open(filepath+race1+'_lengths_unsorted.txt', 'w')
outfile.write(output)
outfile.close()

output_getinfo=''
for l in length_locations:
    output_getinfo+=l
outfile_getinfo=open(filepath+race1+'_position_info.txt', 'w')
outfile_getinfo.write(output_getinfo)
outfile_getinfo.close()

pair_name=race1
sort_lengths=subprocess.Popen(['nice','sort','-n',filepath+pair_name+'_lengths_unsorted.txt'],stdout=subprocess.PIPE)
sort_out=sort_lengths.communicate()[0]
sortfile=open(filepath+pair_name+'_lengths_sorted.txt','w')
sortfile.write(sort_out)
sortfile.close()
condense=subprocess.Popen(['nice','python','condense_sorted_lengths.py',filepath+pair_name+'_lengths'])

sort_pos_info=subprocess.Popen(['nice','sort','-k3n',filepath+race1+'_position_info.txt'],stdout=subprocess.PIPE)
sort_pos_out=sort_pos_info.communicate()[0]
sort_pos_lines=sort_pos_out.split('\n')
sort_pos_file=open(filepath+pair_name+'_position_info_sorted.txt','w')
sort_pos_file.write('Chromosome\t position\t Tract length\n')
for ind in range(1,len(sort_pos_lines)+1):
    sort_pos_file.write(sort_pos_lines[-ind]+'\n')
sort_pos_file.close()
