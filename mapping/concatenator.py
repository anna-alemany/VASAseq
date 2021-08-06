#!/usr/bin/env python3
# Reads R1.fastq and R2.fastq files;
# selects reads with proper cell barcode;
# produces a new _cbc.fastq.gz file.
import sys, os
import itertools as it
import argparse as argp
import numpy as np
import gzip
import pandas as pd
from pandas.io.parsers import read_csv
from collections import Counter
import glob

#### function to identify cells from barcodes, allowing some edit distances ####
def find_compatible_barcodes(barcode, HDmax = 0):
    """Given a barcode sequence and a maximum Hammin distance, it returns a list of compatible barcode sequences"""
    nt = ['N'] if HDmax == 0 else ['N','C','T','G','A']
    HDmax = 1 if HDmax == 0 else HDmax

    compatible_barcodes = set([barcode])
    for hd in range(1, HDmax+1):
        comb = [''.join(l) for l in it.product(nt, repeat = hd)]
        for c in comb:
            for p in it.permutations(range(len(barcode)), hd):
                s0 = barcode
                for x, l in zip(p, c):
                    s0 = s0[:x] + l + s0[x+1:]
                compatible_barcodes.add(s0)
    return list(compatible_barcodes)

#### check input variables ####
parser = argp.ArgumentParser(description = 'Concatenates bcread to bioread qname.')
parser.add_argument('--fqf', help = 'Fastq files names, without _Rx.fastq.gz')
parser.add_argument('--bcread', '-bcr', help = 'read where to find the barcode (umi+cell)', choices = ['R1', 'R2'], default = 'R1')
parser.add_argument('--bioread', '-bior', help = 'read where to find biological information', choices = ['R1', 'R2'], default = 'R2')
parser.add_argument('--demux', '-dx', help = 'print different fastq file for each barcode', action = 'store_true')
parser.add_argument('--lencbc', '-lcbc', help = 'cell barcode length (integer)', type = int, default = 8)
parser.add_argument('--lenumi', '-lumi', help = 'umi length (integer)', type = int, default = 6)
parser.add_argument('--umifirst', help = 'logical variable: umi before cel barcode', action = 'store_true')
parser.add_argument('--cbcfile', '-cbcf', help = 'cell specific barcode file. Please, provide full name')
parser.add_argument('--cbchd', help = 'collapse cell barcodes with the given hamming distance', type = int, default = 0)
parser.add_argument('--outdir', help = 'output directory for cbc.fastq.gz and log files', type = str, default = './')
args = parser.parse_args()

fqr = args.fqf
bcread = args.bcread
bioread = args.bioread
lcbc = args.lencbc
lumi = args.lenumi
umifirst = args.umifirst
cbcfile = args.cbcfile
hd = args.cbchd
outdir = args.outdir
demux = args.demux

#### Find input fastq files ####
fq1s = sorted(glob.glob(fqr + '*_R1*.fastq.gz'))
fq2s = sorted(glob.glob(fqr + '*_R2*.fastq.gz'))
print(fq1s, fq2s)

if len(fq1s) != len(fq2s):
    sys.exit("Please, different number of input and output fastq files")

if len(fq1s) == len(fq2s) == 0:
    sys.exit('fastq files not found')

#### Read barcodes and expand set according to input hamming distance ####
if not os.path.isfile(cbcfile):
    sys.exit("Barcode file not found")

bc_df = read_csv(cbcfile, sep = '\t', names = ['bc','cellID'], index_col = 0)
print(bc_df.head())
bc_df['compatible_bcs'] = bc_df.apply(lambda x: find_compatible_barcodes(x.name, hd), axis = 1)
cnt_allbcs = Counter([x for idx in bc_df.index for x in bc_df.loc[idx, 'compatible_bcs']])
allbc_df = pd.DataFrame({x: {'cellID': bc_df.loc[idx,'cellID'], 'original': idx} for idx in bc_df.index for x in bc_df.loc[idx, 'compatible_bcs'] if cnt_allbcs[x]==1}).T

### Create output directory if it does not exist ####
if not os.path.isdir(outdir):
    os.system('mkdir '+outdir)

#### Read fastq files and assign cell barcode and UMI ####
if not demux:
    fout = open(outdir + '/' + fqr + '_cbc.fastq', 'w')
else:
    fout = {idx:  open(outdir + '/' + fqr + '_' + str(bc_df.loc[idx, 'cellID']).zfill(3) + '_cbc.fastq', 'w') for idx in bc_df.index}

ns = 0; nt = 0
for fq1, fq2 in zip(fq1s, fq2s):
    print(fq1, fq2)
    with gzip.open(fq1) as f1, gzip.open(fq2) as f2: 
        for idx, (l1, l2) in enumerate(zip(f1, f2)):
            try:
                l1, l2 = str(l1.rstrip().rsplit()[0], 'utf-8'), str(l2.rstrip().rsplit()[0], 'utf-8')
            except:
                print('check line '+str(idx))
                l1 = ''; l2 = ''
            l = np.mod(idx,4)
            if l == 0:
                n1, n2 = l1, l2
                if not n1 == n2:
                    print (n1, n2)
                    sys.exit('fastq files not syncrhonized (@name)')
            if l == 1:
                s1, s2 = l1, l2
            if l == 2:
                p1, p2 = l1[0], l2[0]
                if not p1 == p2: # == '+':
                    print(l1, l2, p1, p2)
                    sys.exit('fastq files not synchronized (+)')
            if l == 3 and len(l1) > 0:
                q1, q2 = l1, l2
                nt += 1
                if len(q1) != len(s1) or len(q2) != len(s2):
                    print('phred and read length do not match for '+n1)

                if bcread == 'R1':
                    bcseq = s1[:lumi+lcbc]
                    bcphred = q1[:lumi+lcbc]
                    s1 = s1[lumi+lcbc:]
                    q1 = q1[lumi+lcbc:]
                elif bcread == 'R2':
                    bcseq = s2[:lumi+lcbc]
                    bcphred = q2[:lumi+lcbc]
                    s2 = s2[lumi+lcbc:]
                    q2 = q2[lumi+lcbc:]
                if not umifirst:
                    cellbcseq = bcseq[:lcbc]
                    umiseq = bcseq[lcbc:]
                    cellbcphred = bcphred[:lcbc]
                    umiphred = bcphred[lcbc:]
                else:
                    cellbcseq = bcseq[lumi:]
                    umiseq = bcseq[:lumi]
                    cellbcphred = bcphred[lumi:]
                    umiphred = bcphred[:lumi]

                try:
                    cellID, originalBC = allbc_df.loc[cellbcseq]
                    ns += 1
                    cellbcphred = ''.join([chr(ord(c)+32) for c in cellbcphred])
                    umiphred = ''.join([chr(ord(c)+32) for c in umiphred])

                    name = ';'.join([n1] + [':'.join(x) for x in zip(['SS','CB','QT','RX','RQ','SM'], [cellbcseq, originalBC, cellbcphred, umiseq, umiphred, str(cellID).zfill(3)])])
                    s, q = (s1, q1) if bioread == 'R1' else (s2, q2)
                    if not demux: 
                        fout.write( '\n'.join([name, s, '+', q, '']))
                    else: 
                        fout[cellbcseq].write('\n'.join([name, s, '+', q, '']))
                except: 
                    continue

#nt = (idx+1)/4
if not demux:
    fout.close()
else:
    for cbc in fout:
        fout[cbc].close()

#### LOG ####
fout = open(outdir + '/' + fqr + '.log', 'w')
fout.write('=> to generate cbc file <=\n')
fout.write(', '.join(['fastq file:', str(fqr),'\n']))
fout.write(', '.join(['full barcode in:', str(bcread),'\n']))
fout.write(', '.join(['biological read in:', str(bioread), '\n']))
fout.write(', '.join(['cell specific barcode length:', str(lcbc), '\n']))
fout.write(', '.join(['umi length:', str(lumi), '\n']))
fout.write(', '.join(['umi goes first:', str(umifirst),'\n']))
fout.write(', '.join(['total sequenced reads:', str(nt), '\n']))
fout.write(', '.join(['reads with proper barcodes:', str(ns), str(1.0*ns/nt), '\n']))
fout.close()

#### zip fastq file ####
if not demux: 
    os.system('gzip '+ outdir + '/' + fqr + '_cbc.fastq')
else: 
    os.system('gzip '+ outdir + '/' + fqr + '*_cbc.fastq')


