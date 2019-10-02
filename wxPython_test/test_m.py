#import re
#re.findall("GGTT", "ACGAGGTTTTCCGGAGGTTTTGGCAGGTTAA")

s = r'ACATCATCTATCTATACAATAAAAACTATCCCCTAACTACTACACTACTATCATCACATCATATCACTTTATATCCTAC'

def find_mot_func(st, m, k):

    kmers = []
    for i in range(0,len(m)):
        if i + k <= len(m):
            m1 = m[i:i+k]
            kmers.append(m1) 

    match_profiles = [ st.index(a) for a in kmers if a in st ]

    t0 = match_profiles[0]
    first_mismatch_t = st[t0-1 : t0 + k]
    mis_index = {}    

    n0 = m.find(first_mismatch_t[1:])
    mut_seq = first_mismatch_t[1:]

    for i in range(n0-1, 0, -1):
        mut_seq  = m[i] + mut_seq
        if mut_seq not in st:
            mut_seq = st[st.find(mut_seq[1:]) -1] + mut_seq[1:]
            mis_index[i] = st[st.find(mut_seq[1:]) -1]
        else:
            continue

    motif_result = ""
    for i in range(len(m)):
        if i in mis_index:
            motif_result += mis_index[i]
        else:
            motif_result += m[i]

        if motif_result not in st:
            break

    return motif_result, st.index(motif_result)

## https://bioinformatics.stackexchange.com/questions/594/what-motif-finding-software-is-available-for-multiple-sequences-10kb
## https://github.com/ToHanwei/motif_search/blob/master/search_motif.py
## https://github.com/recervictory/program-repo/blob/master/Python/GreedyMotifSearchWithPseudocounts.py

#x = find_mot_func(s, "CTATCCCCTAACTACTACACTACTA", 9)
motif, motif_site = find_mot_func(s, "CTATCGCTTAACTACTACACTACTA", 9)

print("{0}\t{1}".format(motif_site, motif))

