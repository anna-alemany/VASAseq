#!/usr/bin/env python3
import sys, os
import numpy as np
import pandas as pd
import pysam
from collections import Counter

try:
    bamfile = sys.argv[1]
    stranded = sys.argv[2]
    output = sys.argv[3]
except:
    sys.exit("Please, give: \n(1) input bam file;\n(2) stranded (y/n);\n(3) prefix for output file")

bam = pysam.AlignmentFile(bamfile)
bam_out = pysam.AlignmentFile(output + '.Ribo.bam', template = bam, mode = 'wb')
nreads = 0
nunmapped = 0
nmapped = Counter()
fout = open(output + '.nonRibo.fastq', 'w')
for i, r in enumerate(bam.fetch(until_eof = True)):
    if i == 0:
        reads = [r]
        nreads += 1
    else:
        if r.qname == reads[0].qname:
            reads.append(r)
        else: 
            if all([rs.is_unmapped for rs in reads]): # all umapped => fastq file
                r0 = reads[0]
                r0.qname = '@' + r0.qname
                fout.write('\n'.join([r0.qname, r0.seq, '+', r0.qual]) + '\n')
                nunmapped += 1
            else:
                if stranded == 'n': 
                    mapreads = [x for x in reads if not x.is_unmapped]
                elif stranded == 'y':
                    mapreads = [x for x in reads if (not x.is_unmapped) and (not x.is_reverse)] 
                if len(mapreads) >= 1: # at least one is mapped properly => bam file
                    rgtag = '_'.join(sorted([r.get_tag('RG').rsplit('.')[-1] if 'RG' in [t[0] for t in r.get_tags()] else '-' for r in mapreads]))
                    rgtag = rgtag.replace('-ribo','')
                    nmapped.update([rgtag])
                    r0 = mapreads[0]
                    new_tags = [t if t[0] != 'RG' else ('RG',rgtag) for t in r0.tags]
                    r0.tags = new_tags
                    bam_out.write(r0)
                else:
                    r0 = reads[0]
                    r0.qname = '@' + r0.qname
                    fout.write('\n'.join([r0.qname, r0.seq, '+', r0.qual]) + '\n')
                    nunmapped += 1
            reads = [r]
            nreads += 1

fout.close()
bam.close()
bam_out.close()

fout = open(output + '.ribo-map.log', 'w')
fout.write('Number of reads: ' + str(nreads)+'\n')
fout.write('Number of unmapped reads: '+str(nunmapped)+'\n')
fout.write('Number of mapped reads: '+'\n')
for tag in nmapped:
    fout.write('\t'+tag+': '+str(nmapped[tag])+'\n')
fout.close()
os.system('gzip '+output+'.nonRibo.fastq')
