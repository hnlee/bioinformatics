import dna
import wrappers
import math
import scipy
import scipy.stats

def codon_alignment(sequences):
    translated_seqs = dict(list((seq, dna.translate(sequences[seq]))
                                for seq in sequences))
    protein_alignment = wrappers.run_fsa(translated_seqs)
    nucleotide_alignment = {}
    for sequence in protein_alignment:
        nucleotide_alignment[sequence] = ''
        for pos, aa in enumerate(protein_alignment[sequence]):
            if aa == '-':
                nucleotide_alignment[sequence] += '---'
            else:
                nucleotide_alignment[sequence] += sequences[sequence][(3*pos):(3*pos+3)]
    return (nucleotide_alignment, protein_alignment)

def pairwise(iterable):
    pairs = list(set([tuple(sorted([x, y])) 
                      for x in iterable
                      for y in iterable
                      if x != y]))
    return pairs

def tajima_D(alignment):

    return

synonymous_sites = {'TTT': 0+0+1.0/3, 'TTC': 0+0+1.0/3,
                    'TTA': 1.0/3+0+1.0/3, 'TTG': 1.0/3+0+1.0/3,
                    'CTT': 0+0+1.0, 'CTC': 0+0+1.0,
                    'CTA': 1.0/3+0+1, 'CTG': 1.0/3+0+1.0,
                    'ATT': 0+0+2.0/3, 'ATC': 0+0+2.0/3,
                    'ATA': 0+0+2.0/3, 'ATG': 0+0+0.0,
                    'GTT': 0+0+1.0, 'GTC': 0+0+1.0,
                    'GTA': 0+0+1.0, 'GTG': 0+0+1.0,
                    'TCT': 0+0+1.0, 'TCC': 0+0+1.0,
                    'TCA': 0+0+1.0, 'TCG': 0+0+1.0,
                    'CCT': 0+0+1.0, 'CCC': 0+0+1.0,
                    'CCA': 0+0+1.0, 'CCG': 0+0+1.0,
                    'ACT': 0+0+1.0, 'ACC': 0+0+1.0,
                    'ACA': 0+0+1.0, 'ACG': 0+0+1.0,
                    'GCT': 0+0+1.0, 'GCC': 0+0+1.0,
                    'GCA': 0+0+1.0, 'GCG': 0+0+1.0,
                    'TAT': 0+0+1.0/3, 'TAC': 0+0+1.0/3,
                    'TAA': 0+0+0.0, 'TAG': 0+0+0.0,
                    'CAT': 0+0+1.0/3, 'CAC': 0+0+1.0/3,
                    'CAA': 0+0+1.0/3, 'CAG': 0+0+1.0/3,
                    'AAT': 0+0+1.0/3, 'AAC': 0+0+1.0/3,
                    'AAA': 0+0+1.0/3, 'AAG': 0+0+1.0/3,
                    'GAT': 0+0+1.0/3, 'GAC': 0+0+1.0/3,
                    'GAA': 0+0+1.0/3, 'GAG': 0+0+1.0/3,
                    'TGT': 0+0+1.0/3, 'TGC': 0+0+1.0/3,
                    'TGA': 0+0+0.0, 'TGG': 0+0+0.0,
                    'CGT': 0+0+1.0, 'CGC': 0+0+1.0,
                    'CGA': 1.0/3+0+1, 'CGG': 1.0/3+0+1,
                    'AGT': 0+0+1.0/3, 'AGC': 0+0+1.0/3,
                    'AGA': 1.0/3+0+1.0/3, 'AGG': 1.0/3+0+1.0/3,
                    'GGT': 0+0+1.0, 'GGC': 0+0+1.0,
                    'GGA': 0+0+1.0, 'GGG': 0+0+1.0}

def substitution_rate(codon_x, codon_y):
    pathways = [(x,y,z) for x in range(3) 
                        for y in range(3) 
                        for z in range(3) 
                        if (x != y and y != z and x != z)]
    (sd, nd) = (0.0, 0.0)
    num_paths = len(pathways)
    for path in pathways:
        prev_step = codon_x
        for substitution in path:
            codon = [base for base in prev_step]
            codon[substitution] = codon_y[substitution]
            next_step = ''.join(codon)
            if dna.genetic_code[next_step] == '*':
                num_paths -= 1
                break
            elif next_step == prev_step:
                continue
            elif dna.genetic_code[prev_step] == dna.genetic_code[next_step]:
                sd += 1
            else:  
                nd += 1
            prev_step = next_step
    (sd, nd) = (sd/num_paths, nd/num_paths)
    return (sd, nd)

def nei_gojobori(alignment):
    pairs = pairwise(alignment)
    (dS, dN, dS_dN) = (0.0, 0.0, 0.0)
    for combination in pairs:
        sequence_x = alignment[combination[0]]
        sequence_y = alignment[combination[1]]
        if sequence_x != sequence_y:
            sys.exit("Error in alignment.")
        seq_length = len(sequence_x)
        (S, Sd, Nd, num_codons) = (0.0, 0.0, 0.0, range(seq_length)/3)
        for i in range(seq_length/3):
            codon_x = sequence_x[(3*i):(3*i+3)]
            codon_y = sequence_y[(3*i):(3*i+3)]
            if dna.genetic_code[codon_x] == "*" or dna.genetic_code[codon_y] == "*":
                num_codons -= 1
                next
            if '-' in codon_x or '-' in codon_y:
                num_codons -= 1
                next
            S += sum(synonymous_sites[codon_x], synonymous_sites[codon_y])/2
            (sd, nd) = substitution_rate(codon_x, codon_y)
            Sd += sd
            Nd += nd
        N = 3*num_codons - S
        ds = -(3.0/4)*math.log(1-(4.0/3)*(Sd/S))
        dn = -(3.0/4)*math.log(1-(4.0/3)*(Nd/N))
        dn_ds = dn/ds
        dS += ds
        dN += dn
        dN_dS += dn_ds
    (dS, dN, dN_dS) = map(lambda x: x/len(pairs), (dS, dN, dN_dS))
    return (dS, dN, dN_dS)

def mcdonald_kreitman(alignment, within, between):
    length = max(len(alignment[strain]) for strain in alignment)
    outgroup = [alignment[strain] for strain in between
                if strain in alignment and len(alignment[strain]) == length]
    ingroup = [alignment[strain] for strain in within
               if strain in alignment and len(alignment[strain]) == length]
    (Pn, Ps, Dn, Ds) = (0, 0, 0, 0)
    for codon in range(length/3):
        in_codons = [sequence[(3*codon):(3*codon+3)] for sequence in ingroup]
        out_codons = [sequence[(3*codon):(3*codon+3)] for sequence in outgroup]
        if in_codons == [] or out_codons == []:
            continue
        elif sum(map(lambda x: x not in dna.genetic_code, out_codons)) > 0:
            continue
        elif sum(map(lambda x: x not in dna.genetic_code, in_codons)) > 0:
            continue
        elif len(set(in_codons)) == 1:
            aa = set([dna.genetic_code[bases] for bases in out_codons])
            if in_codons[0] in out_codons:
                continue
            elif dna.genetic_code[in_codons[0]] in aa:
                Ds += 1
            else:
                Dn += 1
        else:
            aa = set([dna.genetic_code[bases] for bases in in_codons])
            if len(aa) > 1:
                Pn += 1
            else:
                Ps += 1
    (NI, fisher_p) = scipy.stats.fisher_exact([[Ds, Ps],[Dn, Pn]])   
    return (Pn, Ps, Dn, Ds, fisher_p, NI)

def F_st():
    return

def HKA():
    return
