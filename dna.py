import sys

complement = {'A':'T', 
              'C':'G', 
              'G':'C', 
              'T':'A'}

def reverse_complement(string):
    sequence = string.upper()
    return ''.join([complement[base] for base in sequence[::-1]])


genetic_code = {'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L',
                'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
                'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M',
                'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
                'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S',
                'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
                'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
                'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
                'TAT':'Y', 'TAC':'Y', 'TAA':'*ochre', 'TAG':'*amber',
                'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q',
                'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K',
                'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E',
                'TGT':'C', 'TGC':'C', 'TGA':'*opal', 'TGG':'W',
                'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
                'AGT':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R',
                'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G'}

def translate(string):
    nucleotides = string.upper()
    nucleotides = ''.join(nucleotides.split('-'))
    amino_acids = ''.join([genetic_code[nucleotides[(3*i):(3*i+3)]]
                           for i in range(len(nucleotides)/3)])
    return amino_acids.split('*')[0] 

def transcribe(string):
    dna = string.upper()
    return dna.replace('T','U')

def reverse_transcribe(string):
    rna = string.upper()
    return rna.replace('U','T')

def read_fasta(filename):
    fasta_file = open(filename, 'r')
    sequences = {}
    sequence_name = ''
    for line in fasta_file:
        if line[0] == '>':
            sequence_name = line[1:-1]
            sequences[sequence_name] = ''
        elif sequence_name == '':
            sys.exit('File is not in proper FASTA format.')
        else:
            sequences[sequence_name] += line[:-1].upper()
    fasta_file.close()
    return sequences

def write_fasta(sequences, filename):
    fasta_file = open(filename, 'w')
    for sequence in sequences:
        fasta_file.write('>%s\n%s\n\n' % (sequence, sequences[sequence]))
    fasta_file.close()
    return
