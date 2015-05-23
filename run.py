import dna
import os


rast_key = dict(list((line.split()[0], line[:-1].split()[-1])
    for line in open('../xanthomonas/annotations/RAST_ID.txt', 'r')
    if 'BR' in line))
#strains = ['DM2_1_12_02', 'ME_Cv_P30', 'RMX815_1a',
#           'Knox623a', 'Knox652c', 'LP217a', 'RMX_24_a_1']
#genomes = [i for i in os.listdir('/Users/hanalee/xanthomonas/alignments/clade3-4/')
#    if '.fas' in i and 'sslist' not in i]    
#for genome in genomes:
#    dna.sqlify_fasta('../xanthomonas/alignments/clade3-4/' + genome,
#        '../xanthomonas/genomes/xanthomonas.sqlite',
#        genome.split('__')[0],
#        short_id=False)
#for strain in rast_key:
#    if strain in ['LMC_P11', 'MEDV_A37', 'MEDV_P39', 'NL_P121', 'NL_P172']:
#        dna.sqlify_rast('../Documents/myRAST/' + rast_key[strain], 
#            '../xanthomonas/genomes/xanthomonas.sqlite', 
#            strain)
#    else:
#        continue
#dna.make_gbk('../psyringae/genomes/psyringae_a5.sqlite',
#     'LP205a',
#     '../psyringae/LP205_a5.gbk',
#     'Pseudomonas syringae')
