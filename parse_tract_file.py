def get_spectrum_freqs(infile):
    lines=infile.readlines()
    lengths=[]
    non_cumul=[]
    cumul=[]
    total_length=0
    for line in lines:
        s=line.split('\t')
        lengths.append(int(s[1]))
        non_cumul.append(int(s[0]))
        total_length+=lengths[-1]*non_cumul[-1]
        cumul.append(non_cumul[-1])
    for i in range(1,len(cumul)):
        cumul[-i-1]+=cumul[-i]
    for i in range(len(cumul)):
        cumul[i]*=1.0/total_length
        non_cumul[i]*=1.0/total_length
    return lengths, non_cumul, cumul

def get_spectrum_counts(infile):
    lines=infile.readlines()
    lengths=[]
    non_cumul=[]
    cumul=[]
    total_length=0
    for line in lines:
        s=line.split('\t')
        lengths.append(int(float(s[1])))
        non_cumul.append(int(s[0]))
        cumul.append(non_cumul[-1])
        total_length+=int(float(s[0]))*int(float(s[1]))
    for i in range(1,len(cumul)):
        cumul[-i-1]+=cumul[-i]
    return lengths, non_cumul, cumul, total_length
