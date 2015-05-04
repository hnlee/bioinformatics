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
    cursor.execute("""CREATE TABLE """ + tblname + """(
                      sequence_name text, 
                      sequence text, 
                      sequence_order integer
                      )""")
    fasta_file = open(filename, 'r')
    sequence_name = ''
    sequence_order = 0
    sequence = ''
    for line in fasta_file:
        if line[0] == '>':
            if sequence != '':
                cursor.execute('INSERT INTO ' + tblname + ' VALUES (?,?,?)',
                               (sequence_name, sequence, sequence_order)) 
                conn.commit()
            if short_id:
                sequence_name = line[1:-1].split()[0]
            else:
                sequence_name = line[1:-1]
            sequence_order += 1
        elif sequence_name == '':
            sys.exit('File is not in proper FASTA format.')
        else:
            sequence += line[:-1].upper()
    cursor.execute('INSERT INTO ' + tblname + ' VALUES (?,?,?)',  
                   (sequence_name, sequence, sequence_order))
    conn.commit()
    cursor.execute('CREATE INDEX id_' + tblname + ' ON ' + tblname + '(sequence_name)')
    conn.commit()
    conn.close()
    return dbname

def retrieve_sql(coordinates, dbname, tblname, seq_id):
    conn = sqlite3.connect(dbname)
    cursor = conn.cursor()
    cursor.execute('SELECT * FROM ' + tblname + ' WHERE sequence_name = ?',
                   (seq_id,))
    row = cursor.fetchone()
    if row == None:
        sys.exit('Specify a sequence ID from the FASTA file.')
    sequence = row[1]
    conn.close()
    if coordinates[0] < coordinates[1]:
        return sequence[(coordinates[0]-1):coordinates[1]]
    else:
        return reverse_complement(sequence[(coordinates[1]-1):coordinates[0]])

def sqlify_rast(rastdir, dbname, tblname):
    peg = read_fasta(rastdir + '/peg.fa', short_id=False) 
    rna = read_fasta(rastdir + '/rna.fa', short_id=False)
    functions_tbl = open(rastdir + 'functions.tbl', 'r')
    peg_functions = dict(list((line.split()[0], (line.split()[1], line.split()[-1]))
                              for line in functions_tbl))
    functions_tbl.close()
    rna_tbl = open(rastdir + 'rna.tbl', 'r')
    rna_functions = dict(list((line.split()[0], line.split()[1:])
                               for line in rna_tbl))
    rna_tbl.close()
    annotations = {}
    for i in sorted(peg):
        peg_id = i.split()[0]
        [start, end] = list(int(j) for j in i.split('_')[-2:])
        contig = i.split()[1].split('_')[0]
        aa_seq = peg[i]
        nt_seq = retrieve_sql((start, end), dbname, tblname, contig)    
        if peg_id in peg_functions:
            (fig_id, function) = peg_functions[peg_id]
        else:
            (fig_id, function) = ('', '')
        annotations[peg_id] = (fig_id, contig, start, end, function, nt_seq, aa_seq)
    for i in sorted(rna):
        rna_id = i
        [start, end] = list(int(j) for j in rna_functions[rna_id][1:2])
        contig = rna_functions[rna_id][0]
        function = rna_functions[rna_id][-1]
        aa_seq = ''
        nt_seq = rna[i]
        annotations[rna_id] = (function, contig, start, end, function, nt_seq, aa_seq)
    conn = sqlite3.connect(dbname)
    cursor = conn.cursor()
    cursor.execute("""CREATE TABLE genes_""" + tblname + """(
                      feature_id text, 
                      fig_id text,
                      contig text,
                      start integer,
                      end integer, 
                      function text,
                      nt_seq text,
                      aa_seq text
                      )""")
    cursor.executemany("""INSERT INTO genes_""" + tblname + """(
                          feature_id, 
                          fig_id, 
                          contig, 
                          start, 
                          end, 
                          function, 
                          nt_seq,
                          aa_seq
                          )""", list([i] + annotations[i] for i in sorted(annotations)))
    conn.commit()
    cursor.execute('CREATE INDEX id_genes_' + tblname + ' ON ' + tblname + '(contig)')
    conn.close()
    return annotations

def make_gbk(dbname, tblname, outputname, organism):
    output = open(outputname, 'w')
    conn = sqlite3.connect()
    cursor = conn.cursor()
    contigs = dict(list((row[0], len(row[1])) for row in
                        cursor.execute('SELECT sequence_name, sequence FROM ' + tblname)
    for contig in contigs:
        output.write('FEATURES             Location/Qualifiers\n')
        output.write('     source          1..%i\n' % (contigs[contig])
        output.write('                     /organism=\"%s\"\n' % (organism)
        output.write('                     /mol_type=\"genomic DNA\"\n')
        output.write('                     /strain=\"%s\"\n' % (tblname)
        output.write('                     /chromosome=\"%s\"\n' % (contig)

        for peg in cursor.execute("""SELECT feature_id, start, end, function, aa_seq 
                                     FROM genes_""" + tblname + """ WHERE 
                                     contig=? AND 
                                     feature_id GLOB \'prot_*\'""", (contig,)):
            (feature_id, start, end, function, aa_seq) = peg
            if start <= end: 
                output.write('     CDS             %i..%i\n' % (start, end))
            else:
                output.write('     CDS             complement(%i..%i)\n' % (end, start))
            output.write('                     /protein_id=\"%s\"\n' % (gene))
            output.write('                     /function=\"%s\"\n' % (function))
            if '(EC ' in function:
                ec_number = function.split('EC ')[1].split(')')[0]
                output.write('                     /EC_number=\"%s\"\n' % (ec_number))
            output.write('                     /db_xref=\"NMPDR:%s\"\n' % (fig_id))
            output.write('                     /translation=\"NMPDR:%s\"\n' % (aa_seq))
        for rna in cursor.execute("""SELECT feature_id, start, end, function
                                     FROM genes_""" + tblname + """ WHERE
                                     contig = ? AND
                                     feature_id NOT GLOB \'prot_*\'""", (contig,)):
            (feature_id, start, end, function) = rna
            if 'tRNA' in function:
                rna_type = 'tRNA'
            else:
                rna_type = 'rRNA'
            if start <= end: 
                output.write('     %s             %i..%i\n' % (rna_type, start, end))
            else:
                output.write('     %s             complement(%i..%i)\n' % (rna_type, end, start))
            output.write('                     /function=\"%s\"\n' % (function))
        cursor.execute('SELECT sequence FROM ' + tblname + ' WHERE sequence_name=?', (contig,))
        (contig_seq) = cursor.fetchone()
        output.write('ORIGIN\n%s\n//\n' % (contig_seq))
    return outputname
