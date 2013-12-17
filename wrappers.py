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

def run_muscle(sequences, path=''):
    return

def run_fsa(sequences, path=''):
    command = [fsa, 'tmp.fasta']
    dna.write_fasta(sequences, 'tmp.fasta')
    fsa = subprocess.Popen(command, stdout=subprocess.PIPE)
    fsa_out, fsa_err = fsa.communicate()
    fsa_results = fsa_out.split("\n")
    os.remove('tmp.fasta')
    return fsa_results
