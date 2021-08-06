#!/usr/bin/env python3
import sys, os
from pandas.io.parsers import read_csv
import numpy as np
import pandas as pd
from collections import Counter
import pickle
import gzip

try:
    inputpickle = sys.argv[1]
    output = sys.argv[2]
    protocol = sys.argv[3]
    filt_unigenes = sys.argv[4]
except:
    sys.exit("Please, provide:\n(1) input pickle.gz;\n(2) output file;\n(3) protocol;\n(4) filter uniq genes (y/n)")

cntdf = pickle.load(gzip.open(inputpickle, 'rb'))
cntdf = cntdf[sorted(cntdf.columns)]

#
def sumCounts(x1, x2):
    if type(x1) == dict and type(x2) == dict:
        x = x1
        for umi in x2:
            if umi not in x1:
                x[umi] = x2[umi]
            else:
                x[umi] = x[umi]+x2[umi]
    elif type(x1) == dict and type(x2) != dict:
        x = x1
    elif type(x1) != dict and type(x2) == dict:
        x = x2
    else:
        x = x1
    return x

def aggCounts(xs):
    ax = xs[0]
    if any([type(x) == dict for x in xs]):
        for x in xs[1:]:
            ax = sumCounts(ax,x)
    return ax

# count tables
# total reads
def countTotalReads(x, protocol = 'vasa'):
    if protocol in ['vasa','10x','smartseq_UMI']:
        y = sum([sum(x[u].values()) for u in x]) if type(x) == dict else 0
    elif protocol in ['ramda','smartseq_noUMI']: 
        umi = 'A'
        y = sum(x[umi].values()) if type(x) == dict else 0
    return y

def countTotalUMI(x):
    return len(x) if type(x) == dict else 0

# reads with no introns
def countExonReads(x, protocol = 'vasa'):
    if protocol in ['vasa','10x','smartseq_UMI']:
        y = sum([sum(x[u].values()) for u in x if 'intron' not in ['-'.join(set(k.rsplit('-'))) for k in x[u]]]) if type(x) == dict else 0
    elif protocol in ['ramda','smartseq_noUMI']:
        umi = 'A'
        y = 0
        if type(x) == dict:
            y = sum([x[umi][k] if 'exon' in k else 0 for k in x[umi]])
    return y

def countExonUMI(x):
    return len([u for u in x if 'intron' not in ['-'.join(set(k.rsplit('-'))) for k in x[u]]]) if type(x) == dict else 0

# unspliced reads
def countIntronReads(x, protocol = 'vasa'):
    if protocol in ['vasa','10x', 'smartseq_UMI']:
        y = sum([sum(x[u].values()) for u in x if 'intron' in ['-'.join(set(k.rsplit('-'))) for k in x[u]]]) if type(x) == dict else 0
    elif protocol in ['ramda','smartseq_noUMI']: 
        umi = 'A'
        y = 0
        if type(x) == dict:
            y = sum([x[umi][k] if 'intron' in k else 0 for k in x[umi]])
    return y

def countIntronUMI(x):
    return len([u for u in x if 'intron' in ['-'.join(set(k.rsplit('-'))) for k in x[u]]]) if type(x) == dict else 0

# estimated transcripts
umi = sorted([x for x in cntdf[cntdf.columns[0]] if type(x)==dict][0].keys())[0]
K = 4**len(umi)
def bc2trans(x):
    if x >= K:
        t = np.log(1.-(float(K)-1e-3)/K)/np.log(1.-1./K)
    elif x > 0 and x < K:
        t = np.log(1.-float(x)/K)/np.log(1.-1./K)
    elif x == 0:
        t = 0
    return  t

fout = open(output + '_mapStats.log', 'w')

total_reads_df = cntdf.applymap(lambda x: countTotalReads(x, protocol))
total_umi_df = cntdf.applymap(lambda x: countTotalUMI(x))
transcripts_total_df = total_umi_df.applymap(bc2trans)

fout.write('Total mapped reads:\t' + str(total_reads_df.sum().sum())  + '\n')
fout.write('Dimension of raw dataset:\t' + str(cntdf.shape) + '\n')

total_reads_df.to_csv(output + '_total.ReadCounts.tsv', sep = '\t')
total_umi_df.to_csv(output + '_total.UFICounts.tsv', sep = '\t')
transcripts_total_df.to_csv(output + '_total.TranscriptCounts.tsv', sep = '\t')

# select tRNA and genes indexes:
tRNAs = [idx for idx in cntdf.index if 'tRNA' in idx]
genes = [idx for idx in cntdf.index if idx not in tRNAs]

fout.write('Total reads assigned to tRNA:\t' + str(total_reads_df.loc[tRNAs].sum().sum())  + '\n')
fout.write('Number of uni/multi-tRNA detected:\t' + str(len(tRNAs))  + '\n')

# tRNA tables
cntdf_tRNA = cntdf.loc[tRNAs]
total_reads_tRNA = cntdf_tRNA.applymap(lambda x: countTotalReads(x, protocol))
total_UMI_tRNA = cntdf_tRNA.applymap(countTotalUMI)

total_reads_tRNA['type'] = ['-'.join(sorted(set([t.rsplit('.')[-1] for t in idx.rsplit('-')]))) for idx in total_reads_tRNA.index]
total_reads_tRNA = total_reads_tRNA.groupby('type').aggregate(sum)
total_UMI_tRNA['type'] = ['-'.join(sorted(set([t.rsplit('.')[-1] for t in idx.rsplit('-')]))) for idx in total_UMI_tRNA.index]
total_UMI_tRNA = total_UMI_tRNA.groupby('type').aggregate(sum)

total_reads_tRNA.to_csv(output + '_tRNA.ReadCounts.tsv', sep = '\t')
total_UMI_tRNA.to_csv(output + '_tRNA.UFICounts.tsv', sep = '\t')

fout.write('Number of tRNA detected after collapsing:\t' + str(len(total_reads_tRNA))  + '\n')

# gene tables
uni_genes = [g for g in genes if '-' not in g]
fout.write('Total reads assigned to genes:\t' + str(total_reads_df.loc[genes].sum().sum())  + '\n')
fout.write('Number of uni/multi-genes detected:\t' + str(len(genes))  + '\n')
fout.write('Total reads assigned to uni-genes:\t' + str(total_reads_df.loc[uni_genes].sum().sum())  + '\n')
fout.write('Number of uni-genes:\t' + str(len(uni_genes))  + '\n')

cntdf_genes = cntdf.loc[genes].copy()

def reduceGeneName(gene, uni_genes):
    rg = gene
    if gene.count('-') == 0:
        rg = gene
    else:
        bios = set([x.rsplit('_')[-1] for x in gene.rsplit('-')])
        shortlived = ['miRNA', 'tRNA','MtTrna']
        longstuff = ['lncRNA']
        shortstuff = ['snRNA','snoRNA','MiscRna','scaRNA']
        ribos = ['rRNA','ribozyme']
        if any([b in ribos for b in bios]):
            gene = '-'.join([g for g in gene.rsplit('-') if g.rsplit('_')[-1] in ribos])
            rg = gene
        if any([b not in shortlived for b in bios])  and any([b in shortlived for b in bios]):
            gene = '-'.join([g for g in gene.rsplit('-') if g.rsplit('_')[-1] not in shortlived])
            rg = gene
        if any([b in shortstuff for b in bios]) and any([b not in shortstuff for b in bios]): 
            gene = '-'.join([g for g in gene.rsplit('-') if g.rsplit('_')[-1] in shortstuff])
            rg = gene
        if sum([g in uni_genes for g in gene.rsplit('-')]) == 1:
            rg = [g for g in gene.rsplit('-') if g in uni_genes][0]
            gene = rg
        if gene.count('-') >= 1 and sum([g.rsplit('_')[1][:2]!="Gm" for g in gene.rsplit('-')]) == 1:
            rg = [g for g in gene.rsplit('-') if g.rsplit('_')[1][:2]!="Gm"][0]
    return rg

def fixGeneLabels(xdf):
    for idx in xdf.index:
        if idx != xdf.loc[idx,'new_gene']:
            i = np.argmax([x==xdf.loc[idx,'new_gene'] for x in idx.rsplit('-')])
            for cell in xdf.columns[:-1]:
                if type(xdf.loc[idx,cell]) == dict:
                    for umi in xdf.loc[idx,cell]:
                        try:
                            xdf.loc[idx,cell][umi] = Counter([k.rsplit('-')[i] for k in xdf.loc[idx,cell][umi].elements()])
                        except:
                            a = True
    return xdf

if filt_unigenes == 'y':
    ncells = max(5, round(0.01*len(cntdf.columns)))
    nreads = 1
    uni_genes_filt = np.array(uni_genes)[(total_reads_df.loc[uni_genes] >= nreads).sum(axis = 1) >= ncells]
else:
    uni_genes_filt = uni_genes
print(len(uni_genes_filt))
cntdf_genes['new_gene'] = [reduceGeneName(idx, uni_genes_filt) for idx in cntdf_genes.index]

total_reads_genes_df = total_reads_df.loc[genes]
total_reads_genes_df['new_gene'] = [reduceGeneName(idx, uni_genes_filt) for idx in cntdf_genes.index]

print(len(set(cntdf_genes['new_gene'])))
cntdf_genes = fixGeneLabels(cntdf_genes)
agg_cntdf_genes = cntdf_genes.groupby('new_gene').aggregate(aggCounts)

uni_genes = [g for g in agg_cntdf_genes.index if '-' not in g]
multi_genes = [g for g in agg_cntdf_genes.index if '-' in g]

total_aggreads_df = agg_cntdf_genes.applymap(lambda x: countTotalReads(x))

fout.write('Total reads after gene aggregation:\t' + str(total_aggreads_df.sum().sum())  + '\n')
fout.write('Dimension of gene-dataset after aggregation:\t' + str(total_aggreads_df.shape) + '\n')
fout.write('Total reads assigned to uni-gene after aggregation:\t' + str(total_aggreads_df.loc[uni_genes].sum().sum())  + '\n')
fout.write('Number of uni-genes after aggregation:\t' + str(len(uni_genes)) + '\n')
fout.write('Total reads assigned to multi-gene after aggregation:\t' + str(total_aggreads_df.loc[multi_genes].sum().sum())  + '\n')
fout.write('Number of multi-gene after aggregation:\t' + str(len(multi_genes))  + '\n')

multi_genes_singleLabel = []
multi_genes_multiLabel = []
multi_genes_multiTags = []
for g in multi_genes:
    ks = set(['&'.join(sorted(x[umi].keys())) for x in agg_cntdf_genes.loc[g] if type(x) == dict for umi in x])
    if len(ks) == 1:
        ks = set(list(ks)[0].rsplit('-'))
        if len(ks) == 1:
            multi_genes_singleLabel.append(g)
        else:
            multi_genes_multiLabel.append(g)
    else:
        multi_genes_multiTags.append(g)

fout.write('Total reads assigned to multi-genes after aggregation that have single labels (all exons/introns):\t' + str(total_aggreads_df.loc[multi_genes_singleLabel].sum().sum())  + '\n')
fout.write('Number of multi-gene after aggregation:\t' + str(len(multi_genes_singleLabel))  + '\n')
fout.write('Total reads assigned to multi-genes after aggregation that have multiple labels (exon-intron etc):\t' + str(total_aggreads_df.loc[multi_genes_multiLabel].sum().sum())  + '\n')
fout.write('Number of multi-gene after aggregation that have multiple labels:\t' + str(len(multi_genes_multiLabel))  + '\n')
fout.write('Total reads assigned to multi-genes after aggregation that have distinct tags (exon-intron and intron-intron):\t' + str(total_aggreads_df.loc[multi_genes_multiTags].sum().sum())  + '\n')
fout.write('Number of multi-genes after aggregation that have distinct tags (exon-intron and intron-intron):\t' + str(len(multi_genes_multiTags))  + '\n')

fout.close()
uni_genes += multi_genes_singleLabel

# unique gene tables
cntdf_unigenes = agg_cntdf_genes.loc[uni_genes]
cntdf_multigenes = agg_cntdf_genes.loc[[idx for idx in agg_cntdf_genes.index if idx not in uni_genes]]

total_reads_multigenes = cntdf_multigenes.applymap(lambda x: countTotalReads(x, protocol))
total_UMI_multigenes = cntdf_multigenes.applymap(countTotalUMI)
total_transcripts_multigenes = total_UMI_multigenes.applymap(bc2trans)

total_reads_unigenes = cntdf_unigenes.applymap(lambda x: countTotalReads(x, protocol))
spliced_reads_unigenes = cntdf_unigenes.applymap(lambda x: countExonReads(x, protocol))
unspliced_reads_unigenes = cntdf_unigenes.applymap(lambda x: countIntronReads(x, protocol))
total_UMI_unigenes = cntdf_unigenes.applymap(countTotalUMI)
spliced_UMI_unigenes = cntdf_unigenes.applymap(countExonUMI)
unspliced_UMI_unigenes = cntdf_unigenes.applymap(countIntronUMI)
total_transcripts_unigenes = total_UMI_unigenes.applymap(bc2trans)
spliced_transcripts_unigenes = spliced_UMI_unigenes.applymap(bc2trans)
unspliced_transcripts_unigenes = unspliced_UMI_unigenes.applymap(bc2trans)

total_reads_unigenes.to_csv(output + '_uniaggGenes_total.ReadCounts.tsv', sep = '\t')
total_UMI_unigenes.to_csv(output + '_uniaggGenes_total.UFICounts.tsv', sep = '\t')
total_transcripts_unigenes.to_csv(output + '_uniaggGenes_total.TranscriptCounts.tsv', sep = '\t')

total_reads_multigenes.to_csv(output + '_multiaggGenes_total.ReadCounts.tsv', sep = '\t')
total_UMI_multigenes.to_csv(output + '_multiaggGenes_total.UFICounts.tsv', sep = '\t')
total_transcripts_multigenes.to_csv(output + '_multiaggGenes_total.TranscriptCounts.tsv', sep = '\t')

unspliced_reads_unigenes.to_csv(output + '_uniaggGenes_unspliced.ReadCounts.tsv', sep = '\t')
unspliced_UMI_unigenes.to_csv(output + '_uniaggGenes_unspliced.UFICounts.tsv', sep = '\t')
unspliced_transcripts_unigenes.to_csv(output + '_uniaggGenes_unspliced.TranscriptCounts.tsv', sep = '\t')

spliced_reads_unigenes.to_csv(output + '_uniaggGenes_spliced.ReadCounts.tsv', sep = '\t')
spliced_UMI_unigenes.to_csv(output + '_uniaggGenes_spliced.UFICounts.tsv', sep = '\t')
spliced_transcripts_unigenes.to_csv(output + '_uniaggGenes_spliced.TranscriptCounts.tsv', sep = '\t')

def remove_ENSandGm(gene):
    rg = sorted(set(['_'.join(x.rsplit("_")[1:]) for x in gene.rsplit('-')]))
    xg = [g for g in rg if g[:2] != 'Gm']
    if len(xg) == 0:
        xg = rg
    xg = '-'.join(xg)
    return xg

total_transcripts_unigenes['new_gene2'] = [remove_ENSandGm(idx) for idx in total_transcripts_unigenes.index] 
unspliced_transcripts_unigenes['new_gene2'] = [remove_ENSandGm(idx) for idx in unspliced_transcripts_unigenes.index]
spliced_transcripts_unigenes['new_gene2'] = [remove_ENSandGm(idx) for idx in spliced_transcripts_unigenes.index]

agg_total_transcripts_unigenes = total_transcripts_unigenes.groupby('new_gene2').aggregate(sum)
agg_unspliced_transcripts_unigenes = unspliced_transcripts_unigenes.groupby('new_gene2').aggregate(sum)
agg_spliced_transcripts_unigenes = spliced_transcripts_unigenes.groupby('new_gene2').aggregate(sum)

agg_total_transcripts_unigenes.to_csv(output + '_shortGeneNames_uniaggGenes_total.TranscriptCounts.tsv', sep = '\t')
agg_unspliced_transcripts_unigenes.to_csv(output + '_shortGeneNames_uniaggGenes_unspliced.TranscriptCounts.tsv', sep = '\t')
agg_spliced_transcripts_unigenes.to_csv(output + '_shortGeneNames_uniaggGenes_spliced.TranscriptCounts.tsv', sep = '\t')


