import sys
import sqlite3

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
                'TAT':'Y', 'TAC':'Y', 'TAA':'*', 'TAG':'*',
                'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q',
                'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K',
                'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E',
                'TGT':'C', 'TGC':'C', 'TGA':'*', 'TGG':'W',
                'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
                'AGT':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R',
                'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G'}

def translate(string):
    nucleotides = string.upper()
    nucleotides = ''.join(nucleotides.split('-'))
    amino_acids = ''.join([genetic_code[nucleotides[(3*i):(3*i+3)]]
                           for i in range(len(nucleotides)/3)])
    return amino_acids 

def transcribe(string):
    dna = string.upper()
    return dna.replace('T','U')

def reverse_transcribe(string):
    rna = string.upper()
    return rna.replace('U','T')

def read_fasta(filename, short_id=True):
    fasta_file = open(filename, 'r')
    sequences = {}
    sequence_name = ''
    for line in fasta_file:
        if line[0] == '>':
            if short_id:
                sequence_name = line[1:-1].split()[0]
            else:
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

def retrieve_sequence(coordinates, sequences, seq_id=''):
    if len(sequences) == 1:
        sequence = sequences.values()[0]
    else:
        try:
            sequence = sequences[seq_id]
        except:
            sys.exit('Specify a sequence ID from the FASTA file.')
    if coordinates[0] < coordinates[1]:
        return sequence[(coordinates[0]-1):coordinates[1]]
    else:
        return reverse_complement(sequence[(coordinates[1]-1):coordinates[0]])

def within_interval(query, reference):
    [left, right] = [False, False]
    if query[0] >= min(reference) and query[0] <= max(reference):
        left = True
    if query[1] >= min(reference) and query[1] <= max(reference):
        right = True
    return [left, right]

def sqlify_fasta(filename, dbname, tblname='', short_id=True):
    if tblname == '':
        tblname = filename.split('.')[0]
    conn = sqlite3.connect(dbname)
    cursor = conn.cursor()
    cursor.execute('''CREATE TABLE ''' + tblname + '''(
                        sequence_name text, sequence text
                      );''')
    fasta_file = open(filename, 'r')
    sequence_name = ''
    sequence = ''
    for line in fasta_file:
        if line[0] == '>':
            if sequence != '':
                cursor.execute('''INSERT INTO ''' + tblname + '''(sequence_name, sequence) 
                                  VALUES (
                                    \'''' + sequence_name + '''\',
                                    \'''' + sequence + '''\'
                                  );''')
            if short_id:
                sequence_name = line[1:-1].split()[0]
            else:
                sequence_name = line[1:-1]
        elif sequence_name == '':
            sys.exit('File is not in proper FASTA format.')
        else:
            sequence += line[:-1].upper()
    cursor.execute('''INSERT INTO ''' + tblname + '''(sequence_name, sequence) 
                      VALUES (
                        \'''' + sequence_name + '''\',
                        \'''' + sequence + '''\'
                      );''')
    cursor.execute('''CREATE UNIQUE INDEX id ON ''' + tblname + '''(sequence_name);''')
    conn.commit()
    conn.close()
    return dbname


            
