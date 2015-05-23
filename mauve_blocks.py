import os
import sys
import pdb

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
