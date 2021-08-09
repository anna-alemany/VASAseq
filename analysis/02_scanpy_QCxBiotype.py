#!/usr/bin/env python
# coding: utf-8

# Second round of QC: produces some plots and saves data to scanpy object.
# Creates 3 different scanpy objects: 1 for only spliced counts, another for only unspliced counts, and the 3rd for both spliced and unspliced counts.
# In the 3rd object, spliced and unspliced counts from the same gene are given different gene identifiers.
# Each scanpy object is further split into new objects, in which only counts from different biotypes are used (protein-coding, smallRNA, lncRNA, tRNA)
# For all objects, only genes present in at least 2 cells are kept. 

# ## Libraries
import sys, os
import numpy as np
import pandas as pd
from pandas.io.parsers import read_csv
from collections import Counter
import scanpy as sc
import bbknn
import plot_aautils as plaa
import sc_aautils as scaa
import pickle
import itertools as it
from scipy.stats import fisher_exact as fisher_exact
import glob
import matplotlib.pyplot as plt
from pandarallel import pandarallel
import scipy

try:
    timepoint = sys.argv[1] # E65, E75, E85, E95
    genebody = sys.argv[2] # all/high
except:
    sys.exit("Please, give:\n(1) timepoint (E65, E75, E85, E95);\n(2) gene body coverage (all/high)")

# path to input files for each timepoint (metadata and merged count tables)
inputmetadir = '/hpc/hub_oudenaarden/aalemany/vasaseq/NovaSeq/'+timepoint+'/mergedData/'
if timepoint == 'E65' and genebody == 'all':
    inputfeatherdir = '/hpc/hub_oudenaarden/aalemany/vasaseq/NovaSeq/E65/all_covs/count_tables_filters/'
elif timepoint == 'E75' and genebody == 'all':
    inputfeatherdir = '/hpc/hub_oudenaarden/aalemany/vasaseq/NovaSeq/E75/allCov/count_tables_with_filters/'
elif timepoint == 'E85' and genebody == 'all':: 
    inputfeatherdir = '/hpc/hub_oudenaarden/aalemany/vasaseq/NovaSeq/E85/allCov/count_tables_with_filter/'
elif timepoint == 'E95' and genebody == 'all'::
    inputfeatherdir = '/hpc/hub_oudenaarden/aalemany/vasaseq/NovaSeq/E95/allCov/count_Tables_with_filters/'
elif timepoint == 'E65' and genebody == 'high'
    inputfeatherdir = '/hpc/hub_oudenaarden/aalemany/vasaseq/NovaSeq/E65/high_cov/filtered_tables/'
elif timepoint == 'E75' and genebody == 'high'
    inputfeatherdir = '/hpc/hub_oudenaarden/aalemany/vasaseq/NovaSeq/E75/highCov/count_tables_with_filters/'
elif timepoint == 'E85' and genebody == 'high'
    inputfeatherdir = '/hpc/hub_oudenaarden/aalemany/vasaseq/NovaSeq/E85/highCov/count_tables_with_filters/'
elif timepoint == 'E95' and genebody == 'high'
    inputfeatherdir = '/hpc/hub_oudenaarden/aalemany/vasaseq/NovaSeq/E95/highCov/count_Tables_with_filters/'
    
# path for output file, depending to each timepoint
outdir = '../'+timepoint+'/res_scanpy_'+genebody+'_rawQC_2021'
os.system('mkdir -p '+outdir)

# ## Read cell metadata (output from script 01)
meta_df = read_csv(inputmetadir + '/cell_meta'+timepoint+'_VASA.tsv', sep = '\t', index_col = 0)
meta_df = meta_df[meta_df['tech_filter']&np.invert(meta_df['doublet'])]

# ## biotype information
p2TF = '/hpc/hub_oudenaarden/aalemany/vasaseq/ref_seqs/' # HPC
#p2TF = '/Users/anna/Desktop/vasaseq/ref_genomes/mouse_ensembl_99/TFs_AnimalTFDB3/' # LOCAL
tfs = read_csv(p2TF + '/Mus_musculus_TF.txt', sep = '\t')
cofs = read_csv(p2TF + '/Mus_musculus_TF_cofactors.txt', sep = '\t')

# ## Read transcriptome data
ufile = glob.glob(inputfeatherdir + '/*shortGeneNames_uniaggGenes_unspliced.TranscriptCounts.feather')
if len(ufile) != 1:
    sys.exit("multiple input unspliced files or not at all, "+str(len(ufile)))
udf = pd.read_feather(ufile[0])
sfile = glob.glob(inputfeatherdir + '/*shortGeneNames_uniaggGenes_spliced.TranscriptCounts.feather')
if len(sfile) > 1: 
    sys.exit("multiple input unspliced files or not at all, "+str(len(sfile)))
sdf = pd.read_feather(sfile[0])

udf = udf.set_index(udf.columns[0])
sdf = sdf.set_index(sdf.columns[0])

udf = udf.loc[[idx for idx in udf.index if type(idx)==str]]
sdf = sdf.loc[[idx for idx in sdf.index if type(idx)==str]]

udf = udf[meta_df.index]
sdf = sdf[meta_df.index]

# ## Prepare scanpy objects
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
sc.settings.figdir = outdir
sc.settings.set_figure_params(dpi=80)

adata_s = sc.AnnData(sdf.T) # only spliced transcriptome
adata_u = sc.AnnData(udf.T) # only unspliced transcriptome

udf.index = [idx + '_u' for idx in udf.index]
sdf.index = [idx + '_s' for idx in sdf.index]
usdf = pd.concat([sdf, udf])

adata_su = sc.AnnData(usdf.T) # both spliced and unspliced transcripts, separated identities

udf.index = ['_'.join(idx.rsplit('_')[:-1]) for idx in udf.index]
sdf.index = ['_'.join(idx.rsplit('_')[:-1]) for idx in sdf.index]

# ### Addition of cell metadata and organization of gene metadata
def add_metadata(adata):
    adata.obs = meta_df.loc[adata.obs.index]
    adata.var['id'] = ['-'.join([k.rsplit('_')[0] for k in idx.rsplit('-')]) for idx in adata.var.index]
    if adata.var.index[0].rsplit('_')[-1] in ['s','u']:
        adata.var['biotype'] = ['-'.join([k.rsplit('_')[-1] for k in idx[:-2].rsplit('-')]) for idx in adata.var.index]
    else:
        adata.var['biotype'] = ['-'.join([k.rsplit('_')[-1] for k in idx.rsplit('-')]) for idx in adata.var.index]
    adata.obs['n_counts'] = adata.X.sum(axis=1)
    adata.obs['n_genes'] = (adata.X>0).sum(axis=1)
    adata.var['n_cells'] = (adata.X>0).sum(axis=0)
    adata.obs['sample'] = ['-'.join(idx.rsplit('-')[:-1]) for idx in adata.obs.index]
    adata.var['ubiotype'] = ['-'.join(sorted(set(b.rsplit('-')))) for b in adata.var['biotype']]
    adata.var['mbiotype'] = ['ProteinCoding' if adata.var.loc[idx,'ubiotype']=='ProteinCoding' else (
                          'lncRNA' if adata.var.loc[idx,'ubiotype']=='lncRNA' else (
                          'sncRNA' if adata.var.loc[idx,'ubiotype'] in ["snRNA","snoRNA","MiscRna","scaRNA",'ribozyme','miRNA'] else (
                          'mixed' if '-' in  adata.var.loc[idx,'ubiotype'] else 'other'))) for idx in adata.var.index]
    adata.var['reg'] = adata.var['id'].apply(lambda x: 'TF' if x in list(tfs['Symbol']) else
                                             ('Cof' if x in list(cofs['Symbol']) else '-'))
    return adata

adata_s = add_metadata(adata_s)
adata_u = add_metadata(adata_u)
adata_su = add_metadata(adata_su)

# ## Biotype splitting
def biotype_split(adata):
    adata = {'All': adata}
    for t, bts in [('ProteinCoding',["ProteinCoding"]),('lncRNA',["lncRNA"]),('smallRNA',["snRNA","snoRNA","MiscRna","scaRNA",'ribozyme','miRNA']),('tRNA',['tRNA'])]:
        adata[t] = adata['All'][:, [g in bts for g in adata['All'].var['ubiotype']] ]
    for t, reg in [('TF', ["TF"]),('Cofactor',["Cof"])]:
        adata[t] = adata['All'][:, [g in reg for g in adata['All'].var['reg']] ]

    for k in adata:
        adata[k].obs['n_counts_'+k] = adata[k].X.sum(axis=1)
        adata[k].obs['n_genes_'+k] = (adata[k].X>0).sum(axis=1)

    return adata

adata_s = biotype_split(adata_s)
adata_u = biotype_split(adata_u)
adata_su = biotype_split(adata_su)

# ## Filter data
# To speed things up, keep genes that are detected in more than 2 cells
print('keep genes in more than 2 cells')
def filter_lowGenes(adata):
    for i, k in enumerate(adata):
        print(k)
        print(adata[k].shape)
        adata[k] = adata[k][:,adata[k].var['n_cells']>2]
        print(adata[k].shape)
        print()
    return adata

adata_s = filter_lowGenes(adata_s)
adata_u = filter_lowGenes(adata_u)
adata_su = filter_lowGenes(adata_su)
    
# ### Check transcript count and gene number per cell
def scatterplot_ncountsVSngenes(adatas, outputname):
    fig, maxs = plaa.template_plot(ncols = len(adatas[0]), nrows = len(adatas), figsize=(len(adatas[0])*3*1.6, 3*len(adatas)))
    for j, (axs, adata) in enumerate(zip(maxs, adatas)):
        for k, (ax, n) in enumerate(zip(axs, adata)):
            df = adata[n].obs
            for i, rep in enumerate(sorted(set(df['sample']))):
                cells = df[df['sample']==rep].index
                ax.scatter(df.loc[cells,'n_counts_'+n], df.loc[cells,'n_genes_'+n], s = 5, c = plaa.colors()[i], label = rep, alpha = 0.5)
            ax.grid(False); ax.grid(axis = 'both', ls = '--', c = 'silver', lw = 0.5)
            if j == 0:
                ax.set_title(n)
            if j == len(maxs)-1:
                ax.set_xlabel('number of transcripts')
            if k == 0:
                ax.set_ylabel('number of genes')
    maxs[-1][0].legend(loc = 2, bbox_to_anchor = (-0.1,-0.25), ncol = len(adata))
    fig.savefig(outdir + '/' + outputname, bbox_inches = 'tight')
    plt.close()
    return

scatterplot_ncountsVSngenes([adata_s, adata_u, adata_su],'scatter_s-u-su_ncountsVSngenes.pdf')

def histo_readsThresholds(adatas, outputname):
    fig, maxs = plaa.template_plot(ncols = len(adatas[sorted(adatas.keys())[0]]), nrows = len(adatas), figsize=(len(adatas[sorted(adatas.keys())[0]])*3*1.6, 3*len(adatas)))
    for k, (axs, bio) in enumerate(zip(maxs, adatas)):
        adata = adatas[bio]
        for j, (ax, n) in enumerate(zip(axs, adata)):
            df = adata[n].obs
            xra = np.log10(df['n_counts_'+n]+1)
            for i, rep in enumerate(sorted(set(df['sample']))):
                cells = df[df['sample']==rep].index
                x = np.log10(df.loc[cells,'n_counts_'+n]+1)
                ax.hist(x, bins = 100, range = (xra.min()*0.9, xra.max()*1.1), color = plaa.colors()[i], label = rep, alpha = 0.5)
                h = np.histogram(x, bins = 100, range = (xra.min()*0.9, xra.max()*1.1))
                if i == 0:
                    histdf = pd.DataFrame({rep: h[0]}, index = [h[1][i:i+1].mean() for i in range(100)])
                else:
                    histdf[rep] = h[0]
            ax.grid(False); ax.grid(axis = 'both', ls = '--', c = 'silver', lw = 0.5)
            if k == 0:
                ax.set_title(n)
            if k == len(maxs)-1:
                ax.set_xlabel('number of transcripts')
            histdf.to_csv(outdir + '/' + outputname + '_' + bio + '_' + n + '.hist', sep = '\t')
        axs[0].set_ylabel('frequency')
    maxs[-1][0].legend(loc = 2, bbox_to_anchor = (-0.1,-0.25), ncol = len(adata))
    fig.savefig(outdir + '/' + outputname + '.pdf', bbox_inches = 'tight')
    plt.close()
    return

histo_readsThresholds({'s': adata_s, 'u': adata_u, 'su': adata_su}, 'histogram_s-u-su_ncounts')

def histo_geneThresholds(adatas, outputname, log = False):
    fig, maxs = plaa.template_plot(ncols = len(adatas[sorted(adatas.keys())[0]]), nrows = len(adatas), figsize=(len(adatas[sorted(adatas.keys())[0]])*3*1.6, 3*len(adatas)))
    for k, (axs, bio) in enumerate(zip(maxs, adatas)):
        adata = adatas[bio]
        for j, (ax, n) in enumerate(zip(axs, adata)):
            df = adata[n].obs
            xra = np.log10(df['n_genes_'+n]+1) if log else df['n_genes_'+n]
            for i, rep in enumerate(sorted(set(df['sample']))):
                cells = df[df['sample']==rep].index
                x = np.log10(df.loc[cells,'n_genes_'+n]+1) if log else df.loc[cells,'n_genes_'+n]
                ax.hist(x, bins = 100, range = (xra.min()*0.9, xra.max()*1.1), color = plaa.colors()[i], label = rep, alpha = 0.5)
                h = np.histogram(x, bins = 100, range = (xra.min()*0.9, xra.max()*1.1))
                if i == 0:
                    histdf = pd.DataFrame({rep: h[0]}, index = [h[1][i:i+1].mean() for i in range(100)])
                else:
                    histdf[rep] = h[0]
            ax.grid(False); ax.grid(axis = 'both', ls = '--', c = 'silver', lw = 0.5)
            if k == 0:
                ax.set_title(n)
            if k == len(maxs)-1:
                if log:
                    ax.set_xlabel('number of log(genes)')
                else: 
                    ax.set_xlabel('number of genes')
            histdf.to_csv(outdir + '/' + outputname + '_' + bio + '_' + n + '.hist', sep = '\t')
        axs[0].set_ylabel('frequency')
    maxs[-1][0].legend(loc = 2, bbox_to_anchor = (-0.1,-0.25), ncol = len(adata))
    fig.savefig(outdir + '/' + outputname + '.pdf', bbox_inches = 'tight')
    plt.close()
    return

histo_geneThresholds({'s': adata_s, 'u': adata_u, 'su': adata_su}, 'histogram_s-u-su_ngenes', log = False)
histo_geneThresholds({'s': adata_s, 'u': adata_u, 'su': adata_su}, 'histogram_s-u-su_logngenes', log = True)

def histo_fracBiotype(adatas, outputname):
    fig, maxs = plaa.template_plot(ncols = len(adatas[sorted(adatas.keys())[0]]), nrows = len(adatas), figsize=(len(adatas[sorted(adatas.keys())[0]])*3*1.6, 3*len(adatas)))
    for k, (axs, bio) in enumerate(zip(maxs, adatas)):
        adata = adatas[bio]
        for j, (ax, n) in enumerate(zip(axs, adata)):
            df = adata[n].obs
            xra = df['n_counts_'+n]/df['n_counts']
            for i, rep in enumerate(sorted(set(df['sample']))):
                cells = df[df['sample']==rep].index
                x = df.loc[cells,'n_counts_'+n]/df.loc[cells,'n_counts']
                ax.hist(x, bins = 100, color = plaa.colors()[i], label = rep, alpha = 0.5, range = (0.9*xra.min(), 1.1*xra.max()))
                h = np.histogram(x, bins = 100, range = (0.9*xra.min(), 1.1*xra.max()))
                if i == 0:
                    histdf = pd.DataFrame({rep: h[0]}, index = [h[1][i:i+1].mean() for i in range(100)])
                else:
                    histdf[rep] = h[0]
            ax.grid(False); ax.grid(axis = 'both', ls = '--', c = 'silver', lw = 0.5)
            histdf.to_csv(outdir + '/' + outputname + '_' + bio + '_' + n + '.hist', sep = '\t')
            if k == 0:
                ax.set_title(n)
            if k == len(maxs)-1:
                ax.set_xlabel('fraction of transcripts')
        axs[0].set_ylabel('frequency')
    maxs[-1][0].legend(loc = 2, bbox_to_anchor = (-0.1,-0.25), ncol = len(adata))
    fig.savefig(outdir + '/' + outputname + '.pdf', bbox_inches = 'tight')
    plt.close()
    return

histo_fracBiotype({'s': adata_s, 'u': adata_u, 'su': adata_su}, 'histogram_s-u-su_fracBiotypeCounts')

for n in adata_s:
    adata_s[n].write(outdir + '/VASA_s_raw_'+n+'.h5ad', compression = 'gzip')
    adata_u[n].write(outdir + '/VASA_u_raw_'+n+'.h5ad', compression = 'gzip')
    adata_su[n].write(outdir + '/VASA_su_raw_'+n+'.h5ad', compression = 'gzip')

    adata_s[n].obs.to_csv(outdir + '/VASAobs_s_raw_'+n+'.tsv', sep = '\t')
    adata_u[n].obs.to_csv(outdir + '/VASAobs_u_raw_'+n+'.tsv', sep = '\t')
    adata_su[n].obs.to_csv(outdir + '/VASAobs_su_raw_'+n+'.tsv', sep = '\t')

