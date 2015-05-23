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
        contig = 0
        coordinate = 0
        for (j, k) in enumerate(contig_coordinates):
            if i > k:
                contig = contigs[j-1][0]
                coordinate = i - contig_coordinates[j-1] + 1
                break
        fasta_coordinates += [(contig, coordinate)]
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
            if metadata[1] == '+':
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
                       (p,))
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

def call_polymorphisms(dbname, tblname, outputname):
    return outputname

