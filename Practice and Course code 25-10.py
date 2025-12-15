
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
#Code to find GC content for this E_coli Sequence.
from Bio.Seq import Seq
print(Bio.__version__)
from Bio import Seq
import Bio
print('new')
E_Coli_seq = E_Coli_seq.upper()

A = E_Coli_seq.count('A')
C = E_Coli_seq.count('C')
G = E_Coli_seq.count('G')
T = E_Coli_seq.count('T')
#Sequence is 1913bp long therefore we use the equation GC% = (G + C)/(A + T + G + C) * 100
GC_fraction = (G + C) / (len(E_Coli_seq))
GC_percentage= GC_fraction * 100

print(A)
print(C)
print(G)
print(T)

print('i')
print(E_Coli_seq)



from Bio.Seq import Seq
print(Bio.__version__)


print(E_Coli_seq.translate())




#FINDING ORFS IN A WGS OF ECOLI:


#Parse takes two arguments, the file name (1) and the format the file needs to be dealt with in.
#Seq_Record is the variable name we assign to the E Coli sequence file we are working with.
#Parse() converts raw data to a 'defined object' meaning the data is stored in a defined format with parts named.
#SeqIO is a module in the biopython library-it is for writing and reading sequence files.
#SeqI0.parse() is a way to create a type of data object called SeqRecord which: Reads file line by line.
#Stores header and seq separately under separate variable names, record.id=seq name from FASTA header and record.seq=sequence itself.
#SeqIO.parse() does not create SeqRecord objects immediately.
#It sets up a generator that knows how to create them only if you ask, hence we must use next() despite there being only one seq in our file.
#'next' takes the first sequence generated, i.e. the first Seq_Record data object. Since many seq files have multiple sequences,
#genes, etc the SeqIO.parse() code will create multiple data objects. We want the first one so we use next(). In our case its not necessary
#since there is only one seq in our file.

EColi_WGS = next(SeqIO.parse('Data/EColi_WGS.fasta', 'fasta'))
#This line of code takes the sequence part of the SeqRecord object and stores it in the E_Coli_seq variable.
E_Coli_seq = EColi_WGS.seq
#E_Coli_seq is now a seq data object(no metadata), it has only the sequence information no ID hence it is not a SeqRecord data object.
print(EColi_WGS.id)
print(len(E_Coli_seq))
#Translating the EColi Seq file:
#Too many stop codons present, including a few aas in which can't be correct for a random stretch of genomic DNA, hence reading frame is wrong.
#I can find real ORFs by finding stretches from a start codon to a stop codon. These are gene regions, other regions could well be non-coding.
print(E_Coli_seq.translate())


#A RF is the way that the DNA could be read for an entire genome-there are 3 of these possible for the forward strand and 3 possible for the reverse complement.
#ORF on the other hand refers to the sub-regions of the RF that are between a start and a stop codon that could thus be translated.

#Hence to find all (pretty much if we discount overlapping genes) of the genes this E Coli can produce we can translate all 6 RFs and then find the ORfs within each Rf which simply requires finding the
# regions between the start and stop codons for each RF...We can use BLAST and visual inspection to decipher if a given AA seq is likely to produce a protein that is actually used in the E Coli
#cell.

from Bio import SeqIO

def codons_in_frame(seq: str, frame: int):
    """
    Return a list of codons (triplets) for a given frame (0, 1, or 2).
    Drops any trailing bases that don't make a full codon.
    """
    if frame not in (0, 1, 2):
        raise ValueError("frame must be 0, 1, or 2")

    seq = seq.upper()
    start = frame
    end = len(seq) - ((len(seq) - start) % 3)  # trim to full codons
#The range function has arguments: range(start,stop,step), hence why we set end to what we do...
    return [seq[i:i+3] for i in range(start, end, 3)]

# ---- load your E. coli genome FASTA (single record) ----
#This is just for the forward strand
E_Coli_forward_seq = str(E_Coli_seq)

# Get codons for RF 0, +1, +2 (forward strand)
rf0  = codons_in_frame(E_Coli_forward_seq, 0)
rf1  = codons_in_frame(E_Coli_forward_seq, 1)
rf2  = codons_in_frame(E_Coli_forward_seq, 2)

print("RF0 first 20 codons:", rf0[:20])
print("RF+1 first 20 codons:", rf1[:20])
print("RF+2 first 20 codons:", rf2[:20])
print("Lengths (codons):", len(rf0), len(rf1), len(rf2))

'''
print(rf0.translate())
print(rf1.translate())
print(rf2.translate())'''
#This won't work since the rf0,1,2 are at the moment lists of codons rather than strings of nucleotides. Thus the translate() function within BioPython won't work-since it is not compatible with lists
#it needs a string to work on...
print('test')
from Bio.Seq import Seq

rf0_for_trans = "".join(rf0)
rf1_for_trans = "".join(rf1)
rf2_for_trans = "".join(rf2)
#The translate() only works on seq objects, it doesnt work on python strings hence why I need Seq()

rf0_F_final=Seq(rf0_for_trans).translate()
print(rf0_F_final)

rf1_F_final=Seq(rf1_for_trans).translate()
print(rf1_F_final)

rf2_F_final=Seq(rf2_for_trans).translate()
print(rf2_F_final)

#Now to find the ORFs:

#Start codons in Bacteria:[ATG],[GTG],[TTG].
#Stop codons in Bacteria:[TAA],[TAG],[TGA].

START_CODONS = {"ATG","GTG","TTG"}
STOP_CODONS = {"TAA", "TAG","TGA" }
ORF_string = ''
for i in range(0, len(rf0_for_trans)-2,3):
    codon = rf0_for_trans[i:i+3]

    if window == 'ATG' or window == 'GTG' or window == 'TTG':
        ORF_string = ''.join(rf0_for_trans[i:i+3])

    elif:
    window == 'TAA', 'TAG', 'TGA'
    break




