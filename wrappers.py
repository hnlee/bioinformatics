import dna
import os
import subprocess

def run_blast(sequences, db, path='', blast_type='blastn', num_hits=1, evalue=10):
    command = [path + blast_type,
               '-query', 'tmp.fasta',
               '-db', db,
               '-culling_limit', str(num_hits),
               #'-best_hit_overhang', '0.1',
               #'-best_hit_score_edge', '0.1',
               # fix this later to work with more than just best hit
               #'-max_target_seqs', str(num_hits),
               '-evalue', str(evalue),
               '-outfmt', '6 qseqid qlen qstart qend sseqid length sstart send evalue']
    dna.write_fasta(sequences, 'tmp.fasta')
    blast = subprocess.Popen(command, stdout=subprocess.PIPE)
    blast_out, blast_err = blast.communicate()
    blast_results = [hit.split("\t") for hit in blast_out.split("\n")]
    os.remove('tmp.fasta')
    return blast_results

def run_prank(sequences, path='', tree='', codon=False, translate=False):
    command = [path + 'prank',
               '-d=tmp.fasta',
               '-o=tmp',
               '-F']
    if tree != '':
        command += ['t=' + tree]
    if codon:
        command += ['-codon']
    if translate:
        command += ['-translate']
    dna.write_fasta(sequences, 'tmp.fasta')
    prank = subprocess.Popen(command, stdout=subprocess.PIPE)
    prank_out, prank_err = prank.communicate()
    prank_results = dna.read_fasta('tmp.best.fas')
    os.remove('tmp.fasta')
    os.remove('tmp.best.fas')
    return prank_results

def convert_prank(sequences, path='', dna_input='', input_format=''):
    command = [path + 'prank',
               '-convert',
               '-d=tmp.fasta',
               '-o=tmp',
               '-keep']
    formats = ['fasta','phylipi','phylips','paml','nexus','raxml']
    if input_format != '':
        if input_format in format:
            command += ['-f=' + input_format]
        else:
            sys.exit('Specify one of the following: %s' % (', '.join(formats)))
    if dna_input != '':
        command += ['dna=' + dna_input]
    else:
        command += ['-translate']
    dna.write_fasta(sequences, 'tmp.fasta')
    prank = subprocess.Popen(command, stdout=subprocess.PIPE)
    prank_out, prank_err = prank.communicate()
    prank_results = dna.read_fasta('tmp.fas')
    os.remove('tmp.fasta')
    os.remove('tmp.fas')
    return prank_results
    


def run_fsa(sequences, path=''):
    command = ['fsa', 'tmp.fasta']
    dna.write_fasta(sequences, 'tmp.fasta')
    fsa = subprocess.Popen(command, stdout=subprocess.PIPE)
    fsa_out, fsa_err = fsa.communicate()
    fsa_file = open('tmp.aln', 'w')
    fsa_file.write(fsa_out)
    fsa_results = dna.read_fasta('tmp.aln')
    os.remove('tmp.fasta')
    os.remove('tmp.aln')
    return fsa_results
