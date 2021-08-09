#!/usr/bin/env python3
import sys, os
from pandas.io.parsers import read_csv
import numpy as np
import pandas as pd
from collections import Counter
import pickle
import gzip 
from itertools import islice
import multiprocess
import glob

try:
    cellFolder = sys.argv[1]
    output = sys.argv[2]
    protocol = sys.argv[3]
    cellidFROMfilename = sys.argv[4]
except:
    sys.exit("Please, provide:\n(1) input folder with all cell files; (2) output file (3) and protocol; (4) cellid from file name (f) or from read name (r)")

if cellidFROMfilename not in ['f','r']:
    sys.exit('parameter 4 needs to be f (from file) or r (from read)')

ncores = 8

def addCount(cnt, gene, umi, label):
    try:
        cnt[gene][umi].update([label])
    except:
        try:
            cnt[gene][umi] =  Counter([label])
        except:
            cnt[gene] = {umi: Counter([label])}
    return cnt

def get_UMI(name, protocol = 'vasa'):
    if protocol in ['vasa','10x','smartseq_UMI']:
        name_info = {x.rsplit(':')[0]: x.rsplit(':')[1] for x in name.rsplit(';')[1:]}
        umi = name_info['RX']
    elif protocol == 'ramda': 
        umi = 'A'
    elif protocol == 'smartseq_noUMI': 
        name_info = {x.rsplit(':')[0]: x.rsplit(':')[1] for x in name.rsplit(';')[1:]}
        umi = 'A'
    return umi

def gene_assignment(indf):
    rname, protocol, df = indf
    umi = get_UMI(rname, protocol)
    label = ''; gene = ''
 
    df = df[df.nMs == df.nMs.min()]
    if 'IN' in df.jSs.values: 
        df = df[df.jSs == 'IN']

    biotypeWsplicing = ['ProteinCoding','lncRNA','UnprocessedPseudogene','ProcessedPseudogene','TranscribedProcessedPseudogene','TranscribedUnprocessedPseudogene', 'TrCGene','PolymorphicPseudogene', 'UnitaryPseudogene', 'TranscribedUnitaryPseudogene']

    if any([b not in biotypeWsplicing for b in df['Biotype']]): # priority no non-spliced biotypes
        df = df.loc[[idx for idx in df.index if df.loc[idx,'Biotype'] not in biotypeWsplicing]]

    if len(df) == 1:
        gene, label = df[['Gene','Label']].iloc[0]
    else:
        if len(set(df['Gene']))==1:
            gene = df['Gene'].iloc[0]
            label = 'intron' if 'intron' in df['Label'].values else 'exon'
        else: # multiple genes
            # factors to consider: length, exon/intron, biotypes
            gdf = {g: df_g for g, df_g in df.groupby('Gene')}
            ldf = {g: 'intron' if 'intron' in gdf[g]['Label'].values else 'exon' for g in gdf}
            if len(set(ldf.values())) == 1: # if all multimappers are the same label (intron/exon)
                # here we could add a length/biotype layer
                gene = '-'.join(gdf.keys()); label = df.iloc[0]['Label']
            else: # if there are different labels, give priority to exons
                # here we could add length/biotype layer
                gene = '-'.join([l for l in ldf.keys() if ldf[l] == 'exon'])
                label = 'exon'
    return umi, label, gene

def get_cellDict(cell):
    cnt = {}; cellid = 'NaN'
    biotypeWsplicing = ['ProteinCoding','lncRNA','UnprocessedPseudogene','ProcessedPseudogene','TranscribedProcessedPseudogene','TranscribedUnprocessedPseudogene', 'TrCGene','PolymorphicPseudogene', 'UnitaryPseudogene', 'TranscribedUnitaryPseudogene']
    cellfiles = glob.glob(cell + '*_genes.bed.gz')
    cellfiles = [f for f in cellfiles if os.stat(f).st_size != 0] # this does not always work! Delete before hand empty bed.gz files as an alternative. 
    for cellfile in cellfiles:
        df = read_csv(cellfile, sep = '\t' , header = None, names = ['Chromosome','Start','End','Name','Strand','Gene','Info','Length','Cov'], low_memory = False)
        if len(df) == 0:
            continue
        if cellidFROMfilename == 'f':
            cellid = cellfile[:cellfile.index('_cbc')] 
        else:
            cellid = {x.rsplit(':')[0]: x.rsplit(':')[1] for x in df.iloc[0]['Name'].rsplit(';') if ':' in x}['SM']
        df['Gene'] = df.apply(lambda x: x['Gene'].replace('-','.') + '_tRNA' if 'tRNA' in x['Gene'] else x['Gene'], axis = 1)
        df['Label'] = df.apply(lambda x: x['Gene'].rsplit('_')[-1], axis = 1)
        df['Biotype'] = df.apply(lambda x: x['Gene'].rsplit('_')[-2], axis = 1)
        df['jSs'] = df.apply(lambda x: x['Info'].rsplit(';jS:')[-1] if ';jS:' in x['Info'] else 'unknown', axis = 1)
        df['nMs'] = df.apply(lambda x: int(x['Info'].rsplit(';nM:')[1].rsplit(';jS:')[0]), axis = 1)
        df = df[df.Biotype != 'TEC']
        df = df.loc[[idx for idx in df.index if np.invert(df.loc[idx,'Biotype'] not in biotypeWsplicing and df.loc[idx,'jSs'] != 'IN')]]
        df['Gene'] = df['Gene'].apply(lambda x: '_'.join(x.rsplit('_')[:-1]))

        gdf = df.groupby('Name')
        
        for g in gdf.groups: 
            umi, label, xgene = gene_assignment((g, protocol, df.loc[gdf.groups[g]]))
            if xgene != '':
                cnt = addCount(cnt, xgene, umi, label)
    return cellid,cnt
     
cells = glob.glob(cellFolder + '/*.singlemappers_genes.bed.gz')
cells = [c[:-len('.singlemappers_genes.bed.gz')] for c in cells]
pool = multiprocess.Pool(ncores)
gcnt = {}
for (cell,cnt) in pool.imap_unordered(get_cellDict, cells):
    print(cell)
    if len(cnt) > 0:
        gcnt[cell] = cnt

pickle.dump(gcnt, open(output + 'dict.pickle', 'wb'))

cntdf = pd.DataFrame(gcnt)
del gcnt

pickle.dump(cntdf, open(output + '.pickle', 'wb'))
os.system('gzip '+output + '.pickle')

sys.exit()



