import Bio
from Bio import SeqIO


#Find longer and shorter sequences 

id_len_dict = {}
id_seq_dict = {}

def sort_dict(non_sorted):
    sorted_dict = {k: v for k, v in sorted(non_sorted.items(), key=lambda x: x[1])}
    return sorted_dict

def maxmin_dict(non_sorted):
    key_max, val_max = max(non_sorted.items(), key=lambda x: x[1])
    key_min, val_min = min(non_sorted.items(), key=lambda x: x[1])
    return key_max, val_max, key_min, val_min

def sequence_extract_fasta(fasta_files):
    id_len_dict = {}
    id_seq_dict = {}
    with open("dna2.fasta") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            id_len_dict[record.id] = len(record.seq)
            id_seq_dict[record.id] = record.seq
    return sort_dict(id_len_dict), id_seq_dict

#%%
def find_orfs(sequence, frame):
    start_codon = "ATG"
    stop_codons = ["TAA", "TAG", "TGA"]
    orfs = []
    single_orf=""

    # Find all start codon positions
    start_positions = [i for i in range(frame,len(sequence),3) if sequence[i:i+3] == start_codon]

    # For each start position, look for stop codons
    for start in start_positions:
        for i in range(start, len(sequence), 3):
            codon = sequence[i:i+3]
            single_orf+=codon
            # print(len(single_orf))
            if codon in stop_codons:
                # print(codon)
                orfs.append(single_orf)
                single_orf=""
                break
            
            
    return orfs,start_positions
# orfs,start_positions = find_orfs(id_seq_dict['gi|142022655|gb|EQ086233.1|16'],2)

#%%

def find_largest_orf_general(id_seq_dict,frame):
    largest_orf_seq = ''
    largest_orf_id=''
    largest_orf_length = 0
    for key, value in id_seq_dict.items():
        orfs, start_positions = find_orfs(id_seq_dict[key],frame)
        for orf_seq in orfs:
            if len(orf_seq)>largest_orf_length:
                largest_orf_id= key
                largest_orf_seq = orf_seq
                largest_orf_length = len(orf_seq)
    return largest_orf_id, largest_orf_seq, largest_orf_length

def find_largest_orf_id(seq_id, id_seq_dict,frame):
    largest_orf_seq = ''
    largest_orf_start=0
    largest_orf_length = 0        
    orfs, start_positions = find_orfs(id_seq_dict[seq_id],frame)
    for orf_seq,start in zip(orfs,start_positions):
        if len(orf_seq)>largest_orf_length:
            largest_orf_seq = orf_seq
            largest_orf_length = len(orf_seq)
            largest_orf_start = start

    return largest_orf_seq, largest_orf_length, largest_orf_start


#Find longer and shorter sequences 
id_len_dict, id_seq_dict = sequence_extract_fasta("dna.example.fasta")
key_max, val_max, key_min, val_min = maxmin_dict(id_len_dict)

# #identify all ORFs
largest_orf_id, largest_orf_seq, largest_orf_length = find_largest_orf_general(id_seq_dict,0)
print(largest_orf_seq )
print(largest_orf_id, largest_orf_length)
print('----')
largest_orf_seq, largest_orf_length, largest_orf_start = find_largest_orf_id('gi|142022655|gb|EQ086233.1|16', id_seq_dict,2)
print(largest_orf_seq)
print(largest_orf_length, largest_orf_start)
# largest_orf_seq, largest_orf_length, largest_orf_start = find_largest_orf_id(seq_id, id_seq_dict)

#find all repeats
# find_seq_in_all_seqs():

# %%
from collections import defaultdict
def find_all_substrings(sequence, substring):
    positions = []
    index = sequence.find(substring)
    while index != -1:
        positions.append(index)
        index = sequence.find(substring, index + 1)

    return positions


sequence = 'ATGCGGGCCATGCTCCTGCATCGCCGCCTTTCGTTCCACCCGGGCCGGCATCGAGTGATGCCGGCGTTGACGTTTTCGTGGAGTGAGTCAGATGAATCACGCAGCGAATCCCGCCGATC'
length=7
repeat_occurrences = defaultdict(int)
for seq_id, sequence in id_seq_dict.items():
    for start in range(len(sequence)-length-1):
        repeat = sequence[start:start+length]
        repeat_index = find_all_substrings(sequence,repeat)
        if any([i < start for i in repeat_index]):
            pass
        else:
            repeat_occurrences[repeat] = repeat_occurrences[repeat] + len(repeat_index)
maxmin_dict(repeat_occurrences)


# repeat = 'TGC'
# i=0
# sequence = 'ATGCGGGCCATGCTCCTGCATCGCCGCCTTTCGTTCCACCCGGGCCGGCATCGAGTGATGCCGGCGTTGACGTTTTCGTGGAGTGAGTCAGATGAATCACGCAGCGAATCCCGCCGATC'
# while i<len(sequence):
#     sequence.index(repeat, start)


# %%
