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
#Translating nucleotides to aa:

from Bio import SeqIO
codon_to_aa_table = {
    "ATA":"I", "ATC":"I", "ATT":"I", "ATG":"M",
    "ACA":"T", "ACC":"T", "ACG":"T", "ACT":"T",
    "AAC":"N", "AAT":"N", "AAA":"K", "AAG":"K",
    "AGC":"S", "AGT":"S", "AGA":"R", "AGG":"R",
    "CTA":"L", "CTC":"L", "CTG":"L", "CTT":"L",
    "CCA":"P", "CCC":"P", "CCG":"P", "CCT":"P",
    "CAC":"H", "CAT":"H", "CAA":"Q", "CAG":"Q",
    "CGA":"R", "CGC":"R", "CGG":"R", "CGT":"R",
    "GTA":"V", "GTC":"V", "GTG":"V", "GTT":"V",
    "GCA":"A", "GCC":"A", "GCG":"A", "GCT":"A",
    "GAC":"D", "GAT":"D", "GAA":"E", "GAG":"E",
    "GGA":"G", "GGC":"G", "GGG":"G", "GGT":"G",
    "TCA":"S", "TCC":"S", "TCG":"S", "TCT":"S",
    "TTC":"F", "TTT":"F", "TTA":"L", "TTG":"L",
    "TAC":"Y", "TAT":"Y", "TAA":"_", "TAG":"_",
    "TGC":"C", "TGT":"C", "TGA":"_", "TGG":"W",}

import os

print("Running from:", os.getcwd())
print("Files in this directory:", os.listdir())

#Parse takes two arguments, the file name (1) and the format the file needs to be dealt with in.
#Seq_Record is the variable name we assign to the E Coli sequence file we are working with.
#Parse() converts raw data to a 'defined object' meaning the data is stored in a defined format with parts named.
#SeqIO is a module in the biopython library-it is for writing and reading sequence files.
#SeqI0.parse() is a way to create a type of data object called SeqRecord which: #Reads file line by line.
#Stores header and seq separately under separate variable names, record.id=seq name from FASTA header and record.seq=sequence itself.
#SeqIO.parse() does not create SeqRecord objects immediately.
#It sets up a generator that knows how to create them only if you ask, hence we must use next() despite there being only one seq in our file.
#'next' takes the first sequence generated, i.e. the first Seq_Record data object. Since many seq files have multiple sequences,
#genes, etc the SeqIO.parse() code will create multiple data objects. We want the first one so we use next(). In our case its not necessary
#since there is only one seq in our file.

Seq_Record = next(SeqIO.parse('Data/Ecoli_K12_50kb.fasta', 'fasta'))
#This line of code takes the sequence part of the SeqRecord object and stores it in the E_Coli_seq variable.
E_Coli_seq = Seq_Record.seq
#E_Coli_seq is now a seq data object(no metadata), it has only the sequence information no ID hence it is not a SeqRecord data object.
print(Seq_Record.id)
print(len(E_Coli_seq))
#Translating the EColi Seq file:
print(E_Coli_seq.translate())
#Too many stop codons present, including a few aas in which can't be correct for a random stretch of genomic DNA, hence reading frame is wrong.
#I can find real ORFs by finding stretches from a start codon to a stop codon. These are gene regions, other regions could well be non-coding.
k=3
ORF_1_string = ""
#This code assumes the start codon is at i=0 which is not necessarily true, highly unlikely to be true.
'''for i in range(0, len(E_Coli_seq)-k+1,k):
    codon = E_Coli_seq[i:i+k]
    if codon == 'ATG' or codon == 'TTG' or codon == 'GTG':
        string += codon
    elif codon == 'TAA' or codon == 'TAG' or codon == 'TGA':
        break
print(ORF_string)'''






#>Python for Biologists Course: Coursera Course




