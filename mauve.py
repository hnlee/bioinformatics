import os
import sys
import sqlite3
import re

def mauve_coordinates(coordinates, dbname, tblname):
    conn = sqlite3.connect(dbname)
    cursor = conn.cursor()
    contigs = dict(list((row[2], (row[0], len(row[1]))) for row in
                        cursor.execute('SELECT * FROM ' + tblname)))
    contig_coordinates = reduce(lambda x,y: x + contigs[y][1],
                                range(1, len(contigs)+1))
    fasta_coordinates = []
    for i in coordinates:
        for (j, k) in enumerate(contig_coordinates):
            if i <= k:
                if j == 1:
                    fasta_coordinates += [(contigs[j][0], k)]
                else:
                    fasta_coordinates += [(contigs[j][0], k-contig_coordinates[(j-1)])]
    return fasta_coordinates

def sqlify_xmfa(filename, dbname, tblname=''):
    if tblname == '':
        tblname = filename.split('.')[0]
    conn = sqlite3.connect(dbname)
    cursor = conn.cursor()
    cursor.execute("""CREATE IF NOT EXISTS " + tblname + "(
                      block integer,
                      sequence_name text,
                      sequence text,
                      start integer,
                      end integer,
                      )""")
    xmfa_file = open(filename, 'r')
    block = 1
    sequence_names = {}
    sequence_name = ''
    sequence = ''
    (start, end) = (0, 0)
    for line in xmfa_file:
        if line[0] == '#':
            if 'File' in line:
                (number, path) = line[:-1].split()
                name = re.split('__|\.', re.split('\\\|\/', path)[-1])[0]       
                sequence_names[number[9:-4]] = name
        elif line[0] == '>':
            if sequence != '' and (start, end) != (0, 0):
                cursor.execute('INSERT INTO ' + tblname + ' VALUES (?,?,?,?,?,?)',
                               (block, sequence_name, sequence, start, end, strand))
                conn.commit()
            metadata = line[:-1].split()[1:2]
            sequence_name = sequence_names[metadata[0].split(':')[0]]
            if metadata[1] = '+':
                start = int(metadata[0].split(':')[1].split('-')[0])
                end = int(metadata[0].split(':')[1].split('-')[1])
            else:
                end = int(metadata[0].split(':')[1].split('-')[0])
                start = int(metadata[0].split(':')[1].split('-')[1])
        elif line[0] == '=' and (start, end) != (0, 0):
            cursor.execute('INSERT INTO ' + tblname + ' VALUES (?,?,?,?,?,?)',
                           (block, sequence_name, sequence, start, end, strand))
            conn.commit()
            block += 1
        else:
            sequence += line[:-1].upper()
    for p in range(1, block+1):
        cursor.execute('SELECT sequence_name FROM ' + tblname + ' WHERE block=?',
                       (p,)
        block_sequences = [str(row[0]) for row in cursor.fetchall()]
        for q in sequence_names:
            if q not in block_sequences:
                cursor.execute('INSERT INTO ' + tblname + ' VALUES (?,?,?,?,?,?)',
                               (p, q, '', 0, 0, '+'))
                conn.commit()
    cursor.execute('CREATE INDEX id_' + tblname + ' ON ' + tblname + '(block)')
    conn.commit()
    conn.close()
    return dbname

def call_polymorphisms():


#alignment_file = '/Volumes/Kakuichi/bergelson/xanthomonas/alignments/by_clade/clade3-4/clade3-4_mauve_noref.fsa'
#alignment_file = '/Volumes/Kakuichi/bergelson/xanthomonas/alignments/by_clade/clade1-2/clade1-2_mauve_noref.fsa'
alignment_file = '/Volumes/Kakuichi/bergelson/xanthomonas/alignments/all_mauve_noref.fsa'
#alignment_file = '/Volumes/Kakuichi/bergelson/xanthomonas/test/test_mauve.fsa'
alignment = open(alignment_file, 'r')
#output1 = open('clade3-4_blocks.txt', 'w')
#output1 = open('clade1-2_blocks.txt', 'w')
output1 = open('all_blocks.txt', 'w')
#output2 = open('clade3-4_clades.txt', 'w')
#output2 = open('clade1-2_clades.txt', 'w')
output2 = open('all_clades.txt', 'w')

output1.write('Block\tLength\tClades\n')
output2.write('Clade\tBlocks\n')

print 'Reading alignment'

index = {}
block = {}
seq_id = ''
block_num = 0
clades = {}
for line in alignment:
    if line[:9] == '#Sequence' and 'File' in line:
        strain_num = line.split('File')[0][9:]
        strain_name = line.split('__')[1][:4]
        index[strain_num] = strain_name
    elif line[0] == '=':
        block_num += 1
        block_length = max(list(len(block[seq]) for seq in block))
        output1.write('%i\t%i\t%i\n' % (block_num,
                                        block_length,
                                        len(block)))
        for seq in block:
            if seq in clades:
                clades[seq] += 1
            else:
                clades[seq] = 1
        block = {}
    elif line[0] == '>':
        seq_id = ':'.join(line.split(' ')[1:3])
        if seq_id.split(':')[1] == '0-0':
            seq_id = ''
            continue
        else:
            block[seq_id] = ''
    else:
        if seq_id == '':
            continue
        else:
            block[seq_id] += line.rstrip()

for clade in clades:
    output2.write('%s\t%i\n' % (clade, clades[clade]))

output1.close()
output2.close()
alignment.close()
import os
import sys
import pdb

genome_dir = '/Volumes/Kakuichi/bergelson/xanthomonas/alignments/by_clade/clade3-4/'
#genome_dir = '/Volumes/Kakuichi/bergelson/xanthomonas/alignments/by_clade/clade1-2/'
#genome_dir = '/Volumes/Kakuichi/bergelson/xanthomonas/alignments/contigs/'
#genome_dir = '/Volumes/Kakuichi/bergelson/xanthomonas/test/'
genome_files = list(file for file in os.listdir(genome_dir)
                    if '__BR' in file)
                    #if '_40MRL.fa.fas' in file)

def read_genome(file_path):
    print file_path
    genome = open(file_path, 'r')
    genome_seq = reduce(lambda x, y: x + y, 
                        list(line[:-1].upper() for line in genome
                             if line[0] != '>'))
    genome.close()
    return(genome_seq)

print 'Loading genomes'

genomes = dict(list((file_name.split('__')[1][:-4], 
#genomes = dict(list((file_name.split('_')[0], 
                     read_genome(genome_dir + file_name))
                    for file_name in genome_files))
strain_list = sorted(genomes)

alignment_file = '/Volumes/Kakuichi/bergelson/xanthomonas/alignments/by_clade/clade3-4/clade3-4_mauve_noref.fsa'
#alignment_file = '/Volumes/Kakuichi/bergelson/xanthomonas/alignments/by_clade/clade1-2/clade1-2_mauve_noref.fsa'
#alignment_file = '/Volumes/Kakuichi/bergelson/xanthomonas/alignments/all_mauve_noref.fsa'
#alignment_file = '/Volumes/Kakuichi/bergelson/xanthomonas/test/test_mauve.fsa'
alignment = open(alignment_file, 'r')
output = open('clade3-4_snps.txt', 'w')
#output = open('clade1_snps.txt', 'w')
#output = open('all_snps.txt', 'w')
output.write('\t'.join(map(lambda x: x + '_base\t' + x + '_pos', strain_list)))
output.write('\n')

REV_COMPL = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', '-': '-'}

def find_snps(block, strain_index, strain_genomes):
    snps = []
    for pos, bases in enumerate(zip(*list(block[seq] for seq in block))):
        if len(set(bases) - set('-')) <= 1:
            continue
        snp = {}
        for seq_id in block:
            strain = strain_index[seq_id.split(':')[0]]
            [start, end] = seq_id.split(':')[1].split('-')
            if block[seq_id][pos] == '-':
                snp_pos = 'NA'
                snp_base = '-'
            elif seq_id.split(':')[2] == '-':
                offset = block[seq_id][:pos + 1].count('-')
                snp_pos = int(end) - 1 - pos + offset 
                snp_base = REV_COMPL[strain_genomes[strain][snp_pos]]
            else:
                offset = block[seq_id][: pos + 1].count('-')
                snp_pos = int(start) - 1 + pos - offset
                snp_base = strain_genomes[strain][snp_pos]
            snp[strain] = (snp_base, snp_pos)
        snps += [snp]
    return(snps)


def write_snps(snp_list, output, strain_list):
    for snp in snp_list:
        for strain in strain_list:
            if strain != strain_list[0]:
                output.write('\t')
            if strain in snp:
                output.write('\t'.join(map(lambda x: str(x), snp[strain])))
            else:
                output.write('-\tNA')
        output.write('\n')

print 'Reading alignment'

index = {}
block = {}
seq_id = ''
block_num = 0
for line in alignment:
    if line[:9] == '#Sequence' and 'File' in line:
        strain_num = line.split('File')[0][9:]
        strain_name = line.split('__')[1][:4]
        index[strain_num] = strain_name
    elif line[0] == '=':
        block_num += 1
        print "Block", block_num
        write_snps(find_snps(block, index, genomes), output, strain_list)
        block = {}
    elif line[0] == '>':
        seq_id = ':'.join(line.split(' ')[1:3])
        if seq_id.split(':')[1] == '0-0':
            seq_id = ''
            continue
        else:
            block[seq_id] = ''
    else:
        if seq_id == '':
            continue
        else:
            block[seq_id] += line.rstrip()

output.close()
alignment.close()
