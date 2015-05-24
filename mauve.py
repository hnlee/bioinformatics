import os
import sys
import sqlite3
import re

def mauve_coordinates(coordinates, dbname, tblname):
    conn = sqlite3.connect(dbname)
    cursor = conn.cursor()
    contigs = dict(list((row[2], (row[0], len(row[1]))) for row in
                        cursor.execute('SELECT * FROM ' + tblname)))
    contig_lengths = [contigs[n][1] for n in range(1, len(contigs)+1)]
    contig_coordinates = [1] + map(lambda x: sum(contig_lengths[:x+1])+1, range(len(contigs)))
    fasta_coordinates = []
    for i in coordinates:
        contig = 0
        coordinate = 0
        for (j, k) in enumerate(contig_coordinates):
            if i < k:
                contig = contigs[j][0]
                coordinate = i - contig_coordinates[j-1] + 1
                break
        fasta_coordinates += [(contig, coordinate)]
    return fasta_coordinates

def sqlify_xmfa(filename, dbname, tblname=''):
    if tblname == '':
        tblname = filename.split('.')[0]
    conn = sqlite3.connect(dbname)
    cursor = conn.cursor()
    cursor.execute("""CREATE TABLE IF NOT EXISTS """ + tblname + """(
                      block integer,
                      sequence_name text,
                      sequence text,
                      start integer,
                      end integer
                      )""")
    xmfa_file = open(filename, 'r')
    block = 1
    sequence_names = {}
    sequence_name = ''
    sequence = ''
    (start, end) = (0, 0)
    for line in xmfa_file:
        if line[0] == '#':
            if 'File' in line and 'Sequence' in line:
                (number, path) = line[:-1].split()
                name = re.split('__|\.', re.split('\\\|\/', path)[-1])[0]       
                sequence_names[number[9:-4]] = name
        elif line[0] == '>':
            if sequence != '' and (start, end) != (0, 0):
                cursor.execute('INSERT INTO ' + tblname + ' VALUES (?,?,?,?,?)',
                               (block, sequence_name, sequence, start, end))
                conn.commit()
            sequence = ''
            metadata = line[:-1].split()[1:3]
            sequence_name = sequence_names[metadata[0].split(':')[0]]
            if metadata[1] == '+':
                start = int(metadata[0].split(':')[1].split('-')[0])
                end = int(metadata[0].split(':')[1].split('-')[1])
            else:
                end = int(metadata[0].split(':')[1].split('-')[0])
                start = int(metadata[0].split(':')[1].split('-')[1])
        elif line[0] == '=' and (start, end) != (0, 0):
            cursor.execute('INSERT INTO ' + tblname + ' VALUES (?,?,?,?,?)',
                           (block, sequence_name, sequence, start, end))
            conn.commit()
            block += 1
            sequence = ''
        else:
            sequence += line[:-1].upper()
    print "Saved all alignment data"
    for p in range(1, block+1):
        cursor.execute('SELECT sequence_name FROM ' + tblname + ' WHERE block=?',
                       (p,))
        block_sequences = [row[0] for row in cursor.fetchall()]
        for q in sequence_names:
            if sequence_names[q] not in block_sequences:
                cursor.execute('INSERT INTO ' + tblname + ' VALUES (?,?,?,?,?)',
                               (p, sequence_names[q], '', 0, 0))
                conn.commit()
    cursor.execute('CREATE INDEX id_' + tblname + ' ON ' + tblname + '(block)')
    conn.commit()
    print "Indexed alignment table"
    conn.close()
    return dbname

def call_polymorphisms(dbname, tblname):
    conn = sqlite3.connect(dbname)
    cursor = conn.cursor()
    cursor.execute('SELECT DISTINCT sequence_name FROM ' + tblname)
    sequence_names = [row[0] for row in cursor.fetchall()]
    columns = ' text, '.join(sequence_names)
    value_placeholders = ', '.join(list('?' for n in range(len(sequence_names))))
    cursor.execute('CREATE TABLE IF NOT EXISTS snps_' + tblname + """(
                    block text,
                    position integer,
                    """ + columns + """ text
                    )""")
    cursor.execute('CREATE TABLE IF NOT EXISTS indels_' + tblname + """(
                    block text,
                    position integer,
                    """ + columns + """ text
                    )""")
    cursor.execute('CREATE TABLE IF NOT EXISTS blocks_' + tblname + """(
                    block text,
                    """ + columns + """ text
                    )""")
    cursor.execute('SELECT DISTINCT block FROM ' + tblname)
    blocks = [row[0] for row in cursor.fetchall()]
    for p in blocks:
        cursor.execute('SELECT sequence_name, sequence FROM ' + tblname + ' WHERE block=?',
                       (p,))
        alignment = dict(list((row[0], row[1]) for row in cursor.fetchall()))
        print alignment
        alignment_length = max(list(len(alignment[sequence_name]) 
                                    for sequence_name in alignment))
        absent = [sequence_name for sequence_name in alignment
                  if alignment[sequence_name] == '']
        cursor.execute("""INSERT INTO blocks_""" + tblname + """
                          VALUES (?, """ + value_placeholders + """)""",
                        tuple([p] + list(int(sequence_name in absent) 
                                         for sequence_name in sequence_names)))
        conn.commit()
        aligned_positions = zip(*list(alignment[sequence_name] 
                                      for sequence_name in sequence_names
                                      if sequence_name not in absent))
        for q in range(alignment_length):
            alleles = set(aligned_positions[q])
            if len(alleles) == 1:
                continue
            elif len(alleles) == 2 and '-' in alleles:
                cursor.execute("""INSERT INTO indels_""" + tblname + """
                                  VALUES (?, ?, """ + value_placeholders + """)""",
                               tuple([p, q+1] + list(alignment[sequence_name][q]
                                                     for sequence_name in sequence_names)))
                conn.commit()
            else:
                cursor.execute("""INSERT INTO snps_""" + tblname + """
                                  VALUES (?, ?, """ + value_placeholders + """)""",
                               tuple([p, q+1] + list(alignment[sequence_name][q]
                                                     for sequence_name in sequence_names)))
                conn.commit()
    return dbname

