
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
from Bio import SeqIO
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

#Hence to find all (pretty much if we discount overlapping genes) of the genes this E Coli can produce we can translate all 6 RFs and find the ORfs within each Rf which simply requires finding the
# regions between the start and stop codons for each RF...We can use BLAST and visual inspection to decipher if a given AA seq is likely to produce a protein that is actually used in the E Coli
#cell.

from Bio import SeqIO


#Function to create a list of codons for each forward frame 0,1 and 2 within the E Coi Genome.

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
#The translate() only works on seq objects, it doesn't work on python strings hence why I need Seq()
#f=Forward strand

rf0_F_final=Seq(rf0_for_trans).translate()
print(rf0_F_final)

rf1_F_final=Seq(rf1_for_trans).translate()
print(rf1_F_final)

rf2_F_final=Seq(rf2_for_trans).translate()
print(rf2_F_final)

#Now to find the ORFs:

#Start codons in Bacteria:[ATG],[GTG],[TTG].
#Stop codons in Bacteria:[TAA],[TAG],[TGA].
'''
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
'''
#This code above is wrong for various reasons...
#As a reminder the goal is to find all the genes that are expressed in the E Coli genome...
#The properties of

#This code is designed to find all of the rf0 ORFs ie all the open reading frames within the frame shift 0 of the E Coli Genome.

rf0_for_trans = "".join(rf0)
seq = rf0_for_trans
START = {"ATG","GTG","TTG"}
STOP = {"TAA", "TAG","TGA" }
ORF_rf0_list = []
in_rf0_ORF = False
current_rf0_ORF = ""
for i in range(0, len(seq)-2,3):
    codon = seq[i:i+3]

    if codon in START:
        in_rf0_ORF = True
        current_rf0_ORF = codon
    else:
        current_rf0_ORF += codon
    if codon in STOP:
        ORF_rf0_list.append(current_rf0_ORF)
        current_orf = ""
        in_orf = False

print("Number of ORFs found:", len(ORF_rf0_list))
print("First ORF:", ORF_rf0_list[0])


print("first 20", ORF_rf0_list[0:20])
#I need to translate this list of ORFs as nucleotide strings in RF0 into aa sequences-this requires the nt strings becoming seq objects so translate() function can work...

#%%

aa_ORF0 = [str(Seq(orf_nt).translate()) for orf_nt in ORF_rf0_list]
orf0_long_aa=[]
for i in aa_ORF0:
    if len(i) > 300:
        orf0_long_aa.append(i)
print(orf0_long_aa)
#%%
