#!/usr/bin/env python3
# coding: utf-8
# **Description**
# 
# In this notebook we perform the quality checks (QC) of the VASA libraries (mouse embryo data: E6.5, E7.5, E8.5 and E9.5) and identify doublet cells with scrublet. 
# The notebook takes as an input the merged tables (in feather format) for both counts obtained using reads mapping in the whole gene body, 
# or counts obtained from reads that fall into the 80% side of the 5' end of the gene. 
# To run in, write in the terminal `01_qc_checks_HPCversion.py E65` for E6.5, or `01_qc_checks_HPCversion.py E75` for E7.5, etc. 

# ### Libraries
import os, sys
import pandas as pd
import numpy as np
from pandas.io.parsers import read_csv
from collections import Counter
import plot_aautils as plaa
import sc_aautils as scaa
import glob
from scipy.cluster import hierarchy
import scrublet as scr
import random
import matplotlib.pyplot as plt

random.seed(2020)

# input information
try:
    timepoint = sys.argv[1]
except:
    sys.exit("Please, give timepoint (E65, E75, E85, E95)")

# provide paths to all input files for each timepoint
if timepoint == 'E65': 
    outdir = '/hpc/hub_oudenaarden/aalemany/vasaseq/NovaSeq/E65/mergedData/'
    outdirfiles = '/hpc/hub_oudenaarden/aalemany/vasaseq/NovaSeq/E65/mergedData/'
    p2merged = '/hpc/hub_oudenaarden/aalemany/vasaseq/NovaSeq/E65/mergedData/'
    p2all = '/hpc/hub_oudenaarden/aalemany/vasaseq/NovaSeq/E65/all_covs/count_tables_filters/'
    prefFeatherAll = 'NovaSeq_E6.5_allCov_filtered'
    p2hi = '/hpc/hub_oudenaarden/aalemany/vasaseq/NovaSeq/E65/high_cov/filtered_tables/'
    prefFeatherHi = 'NovaSeq_E6.5_highCov_filtered'

elif timepoint == 'E75': 
    outdir = '/hpc/hub_oudenaarden/aalemany/vasaseq/NovaSeq/E75/mergedData/'
    outdirfiles = '/hpc/hub_oudenaarden/aalemany/vasaseq/NovaSeq/E75/mergedData/'
    p2merged = '/hpc/hub_oudenaarden/aalemany/vasaseq/NovaSeq/E75/mergedData/'
    p2all = '/hpc/hub_oudenaarden/aalemany/vasaseq/NovaSeq/E75/allCov/count_tables_with_filters/'
    prefFeatherAll = 'NovaSeq_E7.5_allCov_filtered'
    p2hi= '/hpc/hub_oudenaarden/aalemany/vasaseq/NovaSeq/E75/highCov/count_tables_with_filters/'
    prefFeatherHi = 'NovaSeq_E7.5_highCov_filtered'

elif timepoint == 'E85': 
    outdir = '/hpc/hub_oudenaarden/aalemany/vasaseq/NovaSeq/E85/mergedData/'
    outdirfiles = '/hpc/hub_oudenaarden/aalemany/vasaseq/NovaSeq/E85/mergedData/'
    p2merged = '/hpc/hub_oudenaarden/aalemany/vasaseq/NovaSeq/E85/mergedData/'
    p2all = '/hpc/hub_oudenaarden/aalemany/vasaseq/NovaSeq/E85/allCov/count_tables_with_filter/'
    prefFeatherAll = 'NovaSeq_E8.5_allCov_filtered'
    p2hi = '/hpc/hub_oudenaarden/aalemany/vasaseq/NovaSeq/E85/highCov/count_tables_with_filters/'
    prefFeatherHi = 'NovaSeq_E8.5_highCov_filtered'

elif timepoint == 'E95':
    outdir = '/hpc/hub_oudenaarden/aalemany/vasaseq/NovaSeq/E95/mergedData/'
    outdirfiles = '/hpc/hub_oudenaarden/aalemany/vasaseq/NovaSeq/E95/mergedData/'
    p2merged = '/hpc/hub_oudenaarden/aalemany/vasaseq/NovaSeq/E95/mergedData/'
    p2all = '/hpc/hub_oudenaarden/aalemany/vasaseq/NovaSeq/E95/allCov/count_Tables_with_filters/'
    prefFeatherAll = 'NovaSeq_E9.5_allCov_filtered'
    p2hi = '/hpc/hub_oudenaarden/aalemany/vasaseq/NovaSeq/E95/highCov/count_Tables_with_filters/'
    prefFeatherHi = 'NovaSeq_E9.5_highCov_filtered'

# common filtering parameters for all libraries
minReads = 10**(3.5)
minMap = 0.6
minMap3 = 0.2
maxrRNA = 0.4
maxGenes = 15000

# # Read data
# ### Raw read counts
cellBCfile = p2merged + '/VASA_' + timepoint + '_cellBC.tsv'
rawreads1_df = read_csv(cellBCfile, sep = ',', index_col = 0, usecols = [0,1])

fqreadfile = p2merged + '/VASA_' + timepoint + '_FQreadCounts.tsv'
rawreads2_df = read_csv(fqreadfile, sep = ',', index_col = 0)

riboreadfile = p2merged + '/VASA_' + timepoint + '_Ribo.ReadCounts.tsv'
riboreads_df = read_csv(riboreadfile, sep = ',', index_col = 0)

tmapfiles = p2all + '/' + prefFeatherAll + '_total.ReadCounts.feather'
tmapreads_df = pd.read_feather(tmapfiles)
tmapreads_df = tmapreads_df.set_index('index')

tRNAfiles = p2all + '/' + prefFeatherAll + '_tRNA.ReadCounts.feather'
tRNAreads_df = pd.read_feather(tRNAfiles)
tRNAreads_df = tRNAreads_df.set_index(tRNAreads_df.columns[0])

sreadfiles = p2all + '/' + prefFeatherAll + '_uniaggGenes_spliced.ReadCounts.feather'
sreads_df = pd.read_feather(sreadfiles)
sreads_df = sreads_df.set_index(sreads_df.columns[0])

ureadfiles = p2all + '/' + prefFeatherAll + '_uniaggGenes_unspliced.ReadCounts.feather'
ureads_df = pd.read_feather(ureadfiles)
ureads_df = ureads_df.set_index(ureads_df.columns[0])

tmapfiles = p2hi + '/' + prefFeatherHi + '_total.ReadCounts.feather'
tmapreads3_df = pd.read_feather(tmapfiles)
tmapreads3_df = tmapreads3_df.set_index(tmapreads3_df.columns[0])

tRNAfiles = p2hi + '/' + prefFeatherHi + '_tRNA.ReadCounts.feather'
tRNAreads3_df = pd.read_feather(tRNAfiles)
tRNAreads3_df = tRNAreads3_df.set_index(tRNAreads3_df.columns[0])

sreadfiles = p2hi + '/' + prefFeatherHi + '_uniaggGenes_spliced.ReadCounts.feather'
sreads3_df = pd.read_feather(sreadfiles)
sreads3_df = sreads3_df.set_index(sreads3_df.columns[0])

ureadfiles = p2hi + '/' + prefFeatherHi + '_uniaggGenes_unspliced.ReadCounts.feather'
ureads3_df = pd.read_feather(ureadfiles)
ureads3_df = ureads3_df.set_index(ureads3_df.columns[0])

tRNAfiles = p2all + '/' + prefFeatherAll + '_tRNA.UFICounts.feather'
tRNAtrans_df = pd.read_feather(tRNAfiles)
tRNAtrans_df = tRNAtrans_df.set_index(tRNAtrans_df.columns[0])

sreadfiles = p2all + '/' + prefFeatherAll + '_uniaggGenes_spliced.TranscriptCounts.feather'
strans_df = pd.read_feather(sreadfiles)
strans_df = strans_df.set_index(strans_df.columns[0])
strans_df = strans_df.loc[[idx for idx in strans_df.index if type(idx)==str]]

ureadfiles = p2all + '/' + prefFeatherAll + '_uniaggGenes_unspliced.TranscriptCounts.feather'
utrans_df = pd.read_feather(ureadfiles)
utrans_df = utrans_df.set_index(utrans_df.columns[0])

tRNAfiles = p2hi + '/' + prefFeatherHi + '_tRNA.UFICounts.feather'
tRNAtrans3_df = pd.read_feather(tRNAfiles)
tRNAtrans3_df = tRNAtrans3_df.set_index(tRNAtrans3_df.columns[0])

sreadfiles = p2hi + '/' + prefFeatherHi + '_uniaggGenes_spliced.TranscriptCounts.feather'
strans3_df = pd.read_feather(sreadfiles)
strans3_df = strans3_df.set_index(strans3_df.columns[0])

ureadfiles = p2hi + '/' + prefFeatherHi + '_uniaggGenes_unspliced.TranscriptCounts.feather'
utrans3_df = pd.read_feather(ureadfiles)
utrans3_df = utrans3_df.set_index(utrans3_df.columns[0])

print(rawreads1_df.shape, rawreads2_df.shape, utrans_df.shape, utrans3_df.shape)
if rawreads1_df.shape[0] > rawreads2_df.shape[0]:
#    rawreads1_df = rawreads1_df.loc[[idx for idx in rawreads2_df.index if idx in rawreads1_df.index]]
#    c = rawreads1_df.columns
#    rawreads1_df = rawreads1_df.merge(rawreads2_df, how = 'left', left_index = True, right_index = True).fillna(0)
#    rawreads1_df = rawreads1_df[c]
    rawreads1_df = rawreads1_df.loc[rawreads2_df.index]


# # Prepare cell metadata
mdf = {
    'raw_reads': rawreads1_df,
    'trim_reads': rawreads2_df,
    'dep.ribo_reads': riboreads_df, 
    'mapped_reads': tmapreads_df.sum(),
    'genes': (tmapreads_df>0).sum(),
    'mapped3_reads': tmapreads3_df.sum(),
    'genes3': (tmapreads3_df>0).sum(),
    'tRNA_reads': tRNAtrans_df.sum(),
    'spliced_reads': sreads_df.sum(),
    'unspliced_reads': ureads_df.sum(),
    'spliced_trans': strans_df.sum(),
    'unspliced_trans': utrans_df.sum(),
    'tRNA3_reads': tRNAreads3_df.sum(),
    'spliced3_reads': sreads3_df.sum(),
    'unspliced3_reads': ureads3_df.sum(),
    'spliced3_trans': strans3_df.sum(),
    'unspliced3_trans': utrans3_df.sum(),
    }


for k in mdf:
    if type(mdf[k]) == pd.core.frame.DataFrame:
        mdf[k].columns = [k]
    else:
        mdf[k] = pd.DataFrame(mdf[k], columns = [k]) 

mdf = pd.concat([mdf[k] for k in mdf], axis = 1, join = 'outer')
mdf = mdf.fillna(0)
mdf['sample'] = ['-'.join(idx.rsplit('-')[:-1]) for idx in mdf.index]

print(mdf.head())
mdf.to_csv(outdirfiles + '/cell_meta'+timepoint+'_VASA.tsv', sep = '\t')

# # Exploration of cell qualities
gmdf = {s: df_s for s, df_s in mdf.groupby('sample')}

nc = 4; nr = round(len(gmdf)/nc)
fig, myaxs = plaa.template_plot(ncols = nc, nrows = nr, figsize=(nc*3*1.6, nr*3))
axs = myaxs.reshape(nc*nr)
for i, (ax, k) in enumerate(zip(axs, sorted(gmdf))):
    df = gmdf[k]
    ax.hist(np.log10(df['trim_reads']), color = 'blue', bins = 50, density = False, label = k)
    ax.set_xlabel('log10(Trimmed reads)')
    ax.set_ylabel('Frequency')
    ax.legend()
for j in range(i+1, len(axs)):
    axs[j].axis(False)
plt.show()
fig.savefig(outdir + '/histo_totalCounts_lib.pdf', bbox_inches = 'tight')

nc = 4; nr = round(len(gmdf)/nc)
fig, myaxs = plaa.template_plot(ncols = nc, nrows = nr, figsize=(nc*3*1.6, nr*3))
axs = myaxs.reshape(nc*nr)
for i, (ax, k) in enumerate(zip(axs, sorted(gmdf))):
    df = gmdf[k]
    ax.hist(df['mapped_reads']/df['trim_reads'], color = 'blue', bins = 100, density = False, label = k)
    ax.set_xlabel('Mappability')
    ax.set_ylabel('Frequency')
    ax.legend()
for j in range(i+1, len(axs)):
    axs[j].axis(False)
plt.show()
fig.savefig(outdir + '/histo_mappabilities_lib.pdf', bbox_inches = 'tight')

nc = 4; nr = round(len(gmdf)/nc)
fig, myaxs = plaa.template_plot(ncols = nc, nrows = nr, figsize=(nc*3*1.6, nr*3))
axs = myaxs.reshape(nc*nr)
for i, (ax, k) in enumerate(zip(axs, sorted(gmdf))):
    df = gmdf[k]
    ax.hist(df['mapped3_reads']/df['trim_reads'], color = 'blue', bins = 100, density = False, label = k)
    ax.set_xlabel('Mappability (3'' end)')
    ax.set_ylabel('Frequency')
    ax.legend()
for j in range(i+1, len(axs)):
    axs[j].axis(False)
plt.show()
fig.savefig(outdir + '/histo_mappabilities3_lib.pdf', bbox_inches = 'tight')

nc = 4; nr = round(len(gmdf)/nc)
fig, myaxs = plaa.template_plot(ncols = nc, nrows = nr, figsize=(nc*3*1.6, nr*3))
axs = myaxs.reshape(nc*nr)
for i, (ax, k) in enumerate(zip(axs, sorted(gmdf))):
    df = gmdf[k]
    ax.scatter(df['mapped3_reads']/df['trim_reads'], df['mapped_reads']/df['trim_reads'], s = 5, rasterized = True, label = k, c = 'b')
    ax.set_xlabel('Mappability (3'' end)')
    ax.set_ylabel('Total mappability')
    ax.legend()
for j in range(i+1, len(axs)):
    axs[j].axis(False)
plt.show()
fig.savefig(outdir + '/scatter_mapVSmap3_lib.pdf', bbox_inches = 'tight')

fig, ax = plaa.template_plot()
for i, k in enumerate(sorted(gmdf)):
    df = gmdf[k]
    ax.scatter(df['mapped3_reads']/df['trim_reads'], df['mapped_reads']/df['trim_reads'], s = 5, rasterized = True, label = k, alpha = 0.5)
    ax.set_xlabel('Mappability (3'' end)')
    ax.set_ylabel('Total mappability')
ax.legend(loc = 2, bbox_to_anchor = (1,0.95), ncol = 1)
plt.show()
fig.savefig(outdir + '/scatter_mapVSmap3.pdf', bbox_inches = 'tight')

nc = 4; nr = round(len(gmdf)/nc)
fig, axs = plaa.template_plot(ncols = 2, figsize=(2*3*1.6, 3))
ax = axs[0]
for i, k in enumerate(sorted(gmdf)):
    df = gmdf[k]
    ax.scatter(df['trim_reads'], df['mapped_reads']/df['trim_reads'], s = 5, rasterized = True, label = k, alpha = 0.5)
    ax.set_xlabel('raw reads')
    ax.set_ylabel('Total mappability')
ax.set_xscale('log'); 

ax = axs[1]
for i, k in enumerate(sorted(gmdf)):
    df = gmdf[k]
    ax.scatter(df['trim_reads'], df['mapped3_reads']/df['trim_reads'], s = 5, rasterized = True, label = k, alpha = 0.5)
    ax.set_xlabel('raw reads')
    ax.set_ylabel('3'' mappability')
ax.set_xscale('log');
ax.legend(loc = 2, bbox_to_anchor = (1,1), ncol = 1)
plt.show()
fig.savefig(outdir + '/scatter_mapsVcounts.pdf', bbox_inches = 'tight')

nc = 4; nr = round(len(gmdf)/nc)
fig, myaxs = plaa.template_plot(ncols = nc, nrows = nr, figsize=(nc*3*1.6, nr*3))
axs = myaxs.reshape(nc*nr)
for i, (ax, k) in enumerate(zip(axs, sorted(gmdf))):
    df = gmdf[k]
    ax.hist(df['dep.ribo_reads']/df['trim_reads'], color = 'blue', bins = 100, density = False, label = k)
    ax.set_xlabel('Fraction of rRNA reads')
    ax.set_ylabel('Frequency')
    ax.legend()
for j in range(i+1, len(axs)):
    axs[j].axis(False)
plt.show()
fig.savefig(outdir + '/histo_rRNA_lib.pdf', bbox_inches = 'tight')

nc = 4; nr = round(len(gmdf)/nc)
fig, myaxs = plaa.template_plot(ncols = nc, nrows = nr, figsize=(nc*3*1.6, nr*3))
axs = myaxs.reshape(nc*nr)
for i, (ax, k) in enumerate(zip(axs, sorted(gmdf))):
    df = gmdf[k]
    ax.hist(df['genes'], color = 'blue', bins = 40, density = False, label = k)
    ax.set_xlabel('Number of genes')
    ax.set_ylabel('Frequency')
    ax.legend()
for j in range(i+1, len(axs)):
    axs[j].axis(False)
plt.show()
fig.savefig(outdir + '/histo_genes_lib.pdf', bbox_inches = 'tight')

nc = 4; nr = round(len(gmdf)/nc)
fig, myaxs = plaa.template_plot(ncols = nc, nrows = nr, figsize=(nc*3*1.6, nr*3))
axs = myaxs.reshape(nc*nr)
for i, (ax, k) in enumerate(zip(axs, sorted(gmdf))):
    df = gmdf[k]
    ax.hist(df['genes3'], color = 'blue', bins = 40, density = False, label = k)
    ax.set_xlabel('Number of genes')
    ax.set_ylabel('Frequency')
    ax.legend()
for j in range(i+1, len(axs)):
    axs[j].axis(False)
plt.show()
fig.savefig(outdir + '/histo_genes3_lib.pdf', bbox_inches = 'tight')

nc = 4; nr = round(len(gmdf)/nc)
fig, myaxs = plaa.template_plot(ncols = nc, nrows = nr, figsize=(nc*3*1.6, nr*3))
axs = myaxs.reshape(nc*nr)
for i, (ax, k) in enumerate(zip(axs, sorted(gmdf))):
    df = gmdf[k]
    ax.scatter(df['genes'], df['genes3'], s = 5, rasterized = True, label = k, c = 'b')
    ax.set_xlabel('Number of genes')
    ax.set_ylabel('Number of genes (3''end)')
    ax.legend()
for j in range(i+1, len(axs)):
    axs[j].axis(False)
plt.show()
fig.savefig(outdir + '/scatter_geneVgene3_lib.pdf', bbox_inches = 'tight')

nc = 4; nr = round(len(gmdf)/nc)
fig, ax = plaa.template_plot()
for i, k in enumerate(sorted(gmdf)):
    df = gmdf[k]
    ax.scatter(df['genes'], df['genes3'], s = 5, rasterized = True, label = k, alpha = 0.5)
    ax.set_xlabel('Number of genes')
    ax.set_ylabel('Number of genes (3'' end)')
ax.legend(loc = 2, bbox_to_anchor = (1,1), ncol = 1)
plt.show()
fig.savefig(outdir + '/scatter_geneVSgene3.pdf', bbox_inches = 'tight')

nc = 4; nr = round(len(gmdf)/nc)
fig, axs = plaa.template_plot(ncols = 2, figsize=(2*3*1.6, 3))
ax = axs[0]
for i, k in enumerate(sorted(gmdf)):
    df = gmdf[k]
    ax.scatter(df['mapped_reads'], df['genes'], s = 5, rasterized = True, label = k, alpha = 0.5)
    ax.set_xlabel('Mapped reads')
    ax.set_ylabel('Number of genes')
ax.set_xscale('log'); 

ax = axs[1]
for i, k in enumerate(sorted(gmdf)):
    df = gmdf[k]
    ax.scatter(df['mapped3_reads'], df['genes3'], s = 5, rasterized = True, label = k, alpha = 0.5)
    ax.set_xlabel('Mapped reads (3'' end)')
    ax.set_ylabel('Number of genes (3'' end)')
ax.set_xscale('log');
ax.legend(loc = 2, bbox_to_anchor = (1,1), ncol = 1)
plt.show()
fig.savefig(outdir + '/scatter_genesVmappedreads_lib.pdf', bbox_inches = 'tight')

nc = 4; nr = round(len(gmdf)/nc)
fig, axs = plaa.template_plot(ncols = 2, figsize=(2*3*1.6, 3))
ax = axs[0]
for i, k in enumerate(sorted(gmdf)):
    df = gmdf[k]
    ax.scatter(df['trim_reads'], df['genes'], s = 5, rasterized = True, label = k, alpha = 0.5)
    ax.set_xlabel('Raw reads')
    ax.set_ylabel('Number of genes')
ax.set_xscale('log'); 

ax = axs[1]
for i, k in enumerate(sorted(gmdf)):
    df = gmdf[k]
    ax.scatter(df['trim_reads'], df['genes3'], s = 5, rasterized = True, label = k, alpha = 0.5)
    ax.set_xlabel('Mapped reads (3'' end)')
    ax.set_ylabel('Number of genes (3'' end)')
ax.set_xscale('log');
ax.legend(loc = 2, bbox_to_anchor = (1,1), ncol = 1)
plt.show()
fig.savefig(outdir + '/scatter_genesVrawreads_lib.pdf', bbox_inches = 'tight')

mdf['tech_filter'] = (mdf['trim_reads']>minReads) & (mdf['genes'] < maxGenes) & (mdf['mapped_reads']/mdf['trim_reads'] >= minMap) & (mdf['mapped3_reads']/mdf['trim_reads'] >= minMap3) & (mdf['dep.ribo_reads']/mdf['trim_reads'] <= maxrRNA)

mdf.to_csv(outdirfiles + '/cell_meta'+timepoint+'_VASA.tsv', sep = '\t')

# # Doublets
#mdf = read_csv(outdirfiles + '/cell_meta'+timepoint+'_VASA.tsv', sep = '\t', index_col = 0)
print(strans_df.head())
print(strans_df.head().sum(axis=1))

gmdf = {g: mdf_g for g, mdf_g in mdf.groupby('sample')}

scrubs = {k: [] for k in gmdf}
for k in gmdf:
    sdf = strans_df[gmdf[k].index[gmdf[k]['tech_filter']]]
    sdf = sdf.loc[[idx for idx in sdf.index if type(idx)==str]]
    sdf.index = [idx + '_s' for idx in sdf.index]
    udf = utrans_df[gmdf[k].index[gmdf[k]['tech_filter']]]
    udf = udf.loc[[idx for idx in udf.index if type(idx)==str]]
    udf.index = [idx + '_u' for idx in udf.index]
    tdf = pd.concat([sdf, udf])
    scrubs[k] = scr.Scrublet(tdf.T, sim_doublet_ratio = 20, random_state = 52720)
    print(k, tdf.shape)
    del tdf


doubletScoresPredicted = {lib: scrubs[lib].scrub_doublets() for lib in scrubs}
scrubScores = {lib: pd.DataFrame({'cells': gmdf[lib][gmdf[lib]['tech_filter']].index, 
                                'scrublet_scores': scrubs[lib].doublet_scores_obs_}).set_index('cells') for lib in scrubs}

score_th = {}
for lib in sorted(scrubs):
    scrub = scrubs[lib]
    fig, axs = plaa.template_plot(ncols = 2, nrows = 1, figsize=(2*3*1.6,3))
    axs[0].hist(scrub.doublet_scores_obs_, bins = 80, density = True)
    h = axs[1].hist(scrub.doublet_scores_sim_, bins = 80, density = True)
    hdf = pd.DataFrame({'y': h[0], 'x': [h[1][i:i+2].mean() for i in range(len(h[1])-1)]})
    hidx = hdf.index[hdf['y']==hdf[(hdf['x']<0.5)&(hdf['x']>0)]['y'].min()]
    score = hdf.loc[hidx,'x'].values[0]
    score_th[lib] = score
    for ax,t in zip(axs, ['observed','simulated']):
        ax.set_yscale('log')
        ax.set_xlabel('doublet score')
        ax.set_ylabel('frequency')
        ax.set_title(lib + ', ' + t)
        ax.grid(ls = '--')
        ax.axvline(score, c = 'r')
    plt.savefig(outdir + '/scrub_'+lib+'.pdf', bbox_to_inches = 'tight')

mdf = mdf.merge(pd.DataFrame(pd.concat([scrubScores[k] for k in scrubScores])), how = 'left', right_index = True, left_index = True)

mdf['doublet'] = mdf['scrublet_scores']>0.2

mdf.to_csv(outdirfiles + '/cell_meta'+timepoint+'_VASA.tsv', sep = '\t')

