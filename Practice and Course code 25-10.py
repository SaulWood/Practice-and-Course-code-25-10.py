print("Hello World")
#Cant remember how to use enter ffs
#input('Press ENTER to exit', )


def total_income(salary, principle,capital_growth):
    capital_appreciation = principle * ((capital_growth+100)/100)
    total = capital_appreciation + salary
    return total

print(total_income(25000,70000,10))

k=3
dna_seq = 'TAGCTATTATCGACTATCTATCTATACGTACGTACGGCTAGCTAGTAGAGCTTCGATGATGATGCTGCTAGCTAGCTAGCAGCGTAGCTAGCGCTGCTAGTAGTGTCGACGAGCTACGATCGATTA'
TAG_list = []
for i in range(0,len(dna_seq)-k+1,3):
   window = dna_seq[i:i+3]
   if window == 'TAG':
       TAG_list.append(i)

print(TAG_list)
print('new_ex')
#1. Strings and indexing

dnA_seq = "ATGACCATGATTACGGATTCACTG"
#Take first 5 bases.
print(dnA_seq[:5])
#Take last 3 bases
print(dnA_seq[-3:])
#take 1st to 7th base.
print(dnA_seq[0:7])
#1st argument sets start of indexc 2nd argument sets end of index 3rd argument sets step value.
print(dnA_seq[0:10:2])
#Minus applied to 1st argument-start value counted from reverse direction.
print(dnA_seq[1])
print(dnA_seq[-1])
print(dnA_seq[-5])
print(dnA_seq[-2])
print(dnA_seq[-4])
print(dnA_seq[-3:])
A_Count = 0
for i in range(0,len(dnA_seq)-1,1):
    window = dnA_seq[i]
    if window == 'A':
        A_Count += 1

print(A_Count)
#Write code to split it into codons (triplets):
seq = "ATGACCATGATTACGGATTCACTG"
codon_list = []
k=3
for i in range(0,len(seq)-k+1,3):
    window = seq[i:i+k]
    codon_list.append(window)

print(codon_list)

#Count how many "ATG" codons appear in-frame only.
ATG_codon_count = 0
for i in range(0,len(seq)-k+1,3):
    window = seq[i:i+k]
    if window == 'ATG':
       ATG_codon_count += 1

print(ATG_codon_count)

#Count each nucleotide, store in a dictionary.
#Range function dosent include last index hence range(0,len(seq)) will go from i=0 to i=23 not up to i=24
seq = "ATGACCATGATTACGGATTCACTG"
print(len(seq))
nucleotide_dict = {'A':0, 'C':0, 'G':0, 'T':0}
for i in range(0,len(seq),1):
    window = seq[i]
    if window == 'A':
        nucleotide_dict['A'] += 1
    elif window == 'C':
        nucleotide_dict['C'] += 1
    elif window == 'G':
        nucleotide_dict['G'] += 1
    elif window == 'T':
        nucleotide_dict['T'] += 1

print(nucleotide_dict)

print(nucleotide_dict['A'])

#Working with pandas

import Bio
print(Bio.__version__)
#Reversing dna seq.
dna = "ATGACTAGCATCGATCG"
rev_dna_seq = dna[::-1]
print(rev_dna_seq)

#Transcribe DNA → RNA
dna_seq ='ATCGTAGCTATATCGCATATCGCGTAGTATCGTATATAGC'
dict_replace = {'T':'U'}
rna_list = []
rna_transcription = ''
for i in range(0,len(dna_seq),1):
   nuc = dna_seq[i]
   if nuc == 'T':
       nuc  = dict_replace['T']

   rna_list.append(nuc)

rna_transcription = ''.join(rna_list)
print(rna_list)
print(rna_transcription)


#>Python for Biologists Course: Coursera Course




