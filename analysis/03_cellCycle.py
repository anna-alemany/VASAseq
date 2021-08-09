#!/usr/bin/env python
# coding: utf-8

# Here, we identify cells in S-phase for all timepoints, separatedly and combined, based on histone expression (using only spliced read counts). 
# We use only count tables produced with all gene counts (falling into any position along the gene body)
# We filter the data based on the parameters in filterParams.
# Performs differential gene expression analysis between S-phase and non-S-phase cells, for both spliced and unspliced counts (separately)

# Run script as: 
# ```./03_cellCycle.py ```

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
import lmfit

from filterParams import filterParams_allCov

# ## output folder
outdir = '../cellCycle'
os.system('mkdir -p '+outdir)

# ## Scanpy settings
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
sc.settings.figdir = outdir
sc.settings.set_figure_params(dpi=80)

# ## Read h5ad raw data
def read_h5ad_zipped(adatafile):
    if os.path.isfile(adatafile + '.gz'):
        os.system('gunzip '+adatafile + '.gz')
    adata = sc.read_h5ad(adatafile)
    os.system('gzip '+adatafile)
    return adata

timepoints = ['E65','E75','E85','E95']
folder = 'res_scanpy_all_rawQC_2021' # this is the folder were I keep the raw scanpy objects
bio = 'All'

adatasu = {t: read_h5ad_zipped('../' + t + '/' + folder + '/VASA_su_raw_' + bio + '.h5ad') for t in timepoints}
adatas = {t: read_h5ad_zipped('../' + t + '/' + folder + '/VASA_s_raw_' + bio + '.h5ad') for t in timepoints}
adatau = {t: read_h5ad_zipped('../' + t + '/' + folder + '/VASA_u_raw_' + bio + '.h5ad') for t in timepoints}

# ## Filter data
def filter_BiotypeThresholds(adata, reads_th, fracs_th, sample_out):
    for n in adata:
        adata[n] = adata[n][[s not in sample_out for s in adata[n].obs['sample']],:]
        if n in reads_th:
            th_min = reads_th[n][0]; th_max = reads_th[n][1]
            adata[n] = adata[n][(th_min < adata[n].obs['n_counts'])&(adata[n].obs['n_counts'] < th_max),:]
        if n in fracs_th: 
            th_min = fracs_th[n][0]; th_max = fracs_th[n][1]
            adata[n] = adata[n][(th_min < adata[n].obs['n_counts_'+n]/adata[n].obs['n_counts'])&(adata[n].obs['n_counts_'+n]/adata[n].obs['n_counts'] < th_max),:]
    return adata

for timepoint in timepoints:
    sample_out, reads_su_th, fracs_su_th, n_pca, resolution_s, resolution_su, min_disp = filterParams_allCov(timepoint)
    print(timepoint)
    print('Raw:', adatasu[timepoint].shape)
    adatasu[timepoint] = filter_BiotypeThresholds({bio: adatasu[timepoint]}, reads_su_th, fracs_su_th, sample_out)[bio]
    print('First filter:', adatasu[timepoint].shape)
    print()

    adatas[timepoint] = adatas[timepoint][[idx in adatasu[timepoint].obs.index for idx in adatas[timepoint].obs.index],:]
    adatau[timepoint] = adatau[timepoint][[idx in adatasu[timepoint].obs.index for idx in adatau[timepoint].obs.index],:]

# ## Normalization of scanpy objects and gene selection
def normalize_data(adata):
    for n in adata:
        sc.pp.normalize_per_cell(adata[n], counts_per_cell_after=1e4) 
        sc.pp.log1p(adata[n]) 
        adata[n].raw = adata[n]
    return adata

adatas = normalize_data(adatas)

# ## merge all datasets
mdf = adatas[timepoints[0]]
for i, t in enumerate(timepoints[1:]):
    mdf = scaa.concatSCobj(mdf, adatas[t], names = [timepoints[i], timepoints[i+1]])

adatas['merged'] = mdf

# ## Histone and Mito genes
histones = read_csv('HistoneGenes.tsv', sep = '\t', index_col = 0, header = None)
def tag_genes(adata, genelist, colname = 'histone', multis = 'n'):
    if multis == 'n':
        if adata.var.index[0][-2:] in ['_s','_u']:
            adata.var[colname] = ['_'.join(idx.rsplit('_')[:-1]) in genelist for idx in adata.var.index]
            adata.raw.var[colname] = ['_'.join(idx.rsplit('_')[:-1]) in genelist for idx in adata.raw.var.index]
            adata.raw.var[colname] = ['_'.join(idx.rsplit('_')[:-1]) in genelist for idx in adata.raw.var.index]
        else:
            adata.var[colname] = [idx in genelist for idx in adata.var.index]
            adata.raw.var[colname] = [idx in genelist for idx in adata.raw.var.index]
            adata.raw.var[colname] = [idx in genelist for idx in adata.raw.var.index]
        adata.obs[colname + '_logfraction'] = adata[:,adata.var[colname]].X.sum(axis=1)/adata.X.sum(axis=1)
        adata.obs[colname + '_fraction'] = (np.exp(adata[:,adata.var[colname]].X)-1).sum(axis=1)/(np.exp(adata.X)-1).sum(axis=1)
    elif multis == 'y':
        if adata.var.index[0][-2:] in ['_s','_u']:
            adata.var[colname] = [any([g in genelist for g in idx[:-2].rsplit('-')]) for idx in adata.var.index]
            adata.raw.var[colname] = [any([g in genelist for g in idx[:-2].rsplit('-')]) for idx in adata.raw.var.index]
        else:
            adata.var[colname] = [any([g in genelist for g in idx.rsplit('-')]) for idx in adata.var.index]
            adata.raw.var[colname] = [any([g in genelist for g in idx.rsplit('-')]) for idx in adata.raw.var.index]
        adata.obs[colname + '_logfraction'] = adata[:,adata.var[colname]].X.sum(axis=1)/adata.X.sum(axis=1)
        adata.obs[colname + '_fraction'] = (np.exp(adata[:,adata.var[colname]].X)-1).sum(axis=1)/(np.exp(adata.X)-1).sum(axis=1)
    return adata

def mito_genes(adata):
    adata.var['mito'] = ['mt.' in idx for idx in adata.var.index]
    adata.obs['mito_logfraction'] = adata[:,adata.var['mito']].X.sum(axis=1)/adata.X.sum(axis=1)
    return adata

for t in adatas: 
    adatas[t] = tag_genes(adatas[t], list(histones[1]))
    adatas[t] = mito_genes(adatas[t])

# ## Classify cells according to Histone-gene fraction (S-phase True/False)
def bimodal(x, A, s1, s2, mu1, mu2):
    f = A*np.exp(-0.5*(x-mu1)**2/s1)/np.sqrt(2*np.pi*s1) + (1.-A)*np.exp(-0.5*(x-mu2)**2/s2)/np.sqrt(2*np.pi*s2)
    return f

fig, axs = plaa.template_plot(ncols = len(adatas), figsize = (len(adatas)*3*1.6, 3))
for ax, t in zip(axs, adatas.keys()):
    adata = adatas[t]
    ax.grid(False)
    ax.set_title(t)
    v = adata.obs['histone_logfraction']
    hist = ax.hist(v, bins = 100, density = True)
    ax.set_xlabel('Histone total log-fraction')
    # find threshold using fits
    histdf = pd.DataFrame({'y': hist[0], 'x': [hist[1][i:i+1].mean() for i in range(100)]})
    histdf.iloc[1:]
    bimod = lmfit.Model(bimodal)
    bimod.set_param_hint('s1', value = v.var(), min = 1e-5)
    bimod.set_param_hint('s2', value = v.var(), min = 1e-5)
    bimod.set_param_hint('A', value = 0.5, min = 0, max = 1)
    bimod.set_param_hint('mu1', value = 0.5*(v.mean()+v.max()), min = 1e-5, max = v.max())
    bimod.set_param_hint('mu2', value = 0.5*(v.mean()+v.min()), min = 1e-5, max = v.min())
    bimod.make_params()
    histdf = histdf[histdf['y']>0]
    bimod_fit = bimod.fit(histdf['y'], x = histdf['x'], mu1 = 0.5*(v.mean()+v.max()), mu2 = 0.5*(v.mean()+v.min()), s1 = v.var(), s2 = v.var(), A = 0.5)
    histdf['best_fit'] = bimod_fit.best_fit
    ax.plot(histdf['x'], bimod_fit.best_fit, c = 'k', ls = '--')
    ax.axvline(bimod_fit.best_values['mu1'], c = 'r'); ax.axvline(bimod_fit.best_values['mu2'], c = 'r');
    histdf = histdf[histdf['x']<bimod_fit.best_values['mu1']]
    th = histdf.iloc[np.where(histdf['best_fit']==histdf['best_fit'].min())]['x'].values[0]
    adata.uns['threshold_histone_logfraction'] = th
    ax.axvline(th, lw = 2, c = 'indigo', ls = '--')
    ax.text(th, 40, '%.4f' % th)
    # annotate S-phase cells
    adata.obs['S-phase'] = adata.obs['histone_logfraction']>th
    adata.obs['S-phase'] = adata.obs['S-phase'].astype(str)
    adata.obs['S-phase'] = adata.obs['S-phase'].astype('category')
axs[0].set_ylabel('density')
fig.savefig(outdir + '/all_histos_HistoneFraction.pdf', bbox_inches = 'tight')

fig, axs = plaa.template_plot(ncols = len(adatas), figsize = (len(adatas)*3*1.6, 3))
for ax, t in zip(axs, adatas):
    adata = adatas[t]
    ax.grid(False)
    ax.set_title(t)
    ax.hist(adata.obs['mito_logfraction'], bins = 100)
    ax.set_xlabel('MT fraction');
axs[0].set_ylabel('frequency')
fig.savefig(outdir + '/histos_MT-Fraction.pdf', bbox_inches = 'tight')

# ## reclassify merged, because the fit does not work well there
for t in timepoints: 
    for idx in adatas[t].obs.index:
        adatas['merged'].obs.loc[idx,'S-phase'] = adatas[t].obs.loc[idx,'S-phase']

# ## Differential gene expression analysis in S-phase cells
dexs_cycling = {t: scaa.difGeneExpr(adatas[t], ['S-phase'])['S-phase'] for t in adatas}
with pd.ExcelWriter(outdir + '/dex_cyclingGenes.xlsx') as writer:
    for t in dexs_cycling:
        dexs_cycling[t]['True'].sort_values(by='logfoldchanges', ascending=False).to_excel(writer, sheet_name = str(t))

# ## Final list of genes
cycling_genes = []
lfc = 0.25
for t in dexs_cycling:
    print(t)
    genes = dexs_cycling[t]['True'][dexs_cycling[t]['True']['logfoldchanges']>lfc].index
    for idx in genes:
        for g in idx.rsplit('-'):
            if len(g) > 3:
                c = [x for x in adatas[t].var.index if g in x and '-' not in x]
                if len(c) > 0:
                    cycling_genes += c
                else:
                    cycling_genes += [idx]
                if len(c) != 1: 
                    print(g,c,idx)
cycling_genes = sorted(set(cycling_genes))
cgdf = pd.DataFrame(pd.Series(cycling_genes, name = 'S-phase_genes'))
cgdf.to_csv(outdir + '/Sphase_genes.tsv', sep = '\t')

# ## Check if there are interesting unspliced genes

for t in timepoints:
    adatau[t].obs['S-phase'] = adatas[t].obs.loc[adatau[t].obs.index, 'S-phase']

mdf = adatau[timepoints[0]]
for i, t in enumerate(timepoints[1:]):
    mdf = scaa.concatSCobj(mdf, adatau[t], names = [timepoints[i], timepoints[i+1]])
adatau['merged'] = mdf

dexs_cycling_u = {t: scaa.difGeneExpr(adatau[t], ['S-phase'])['S-phase'] for t in adatau}

with pd.ExcelWriter(outdir + '/dex_cyclingGenes_unspliced.xlsx') as writer:
    for t in dexs_cycling_u:
        dexs_cycling_u[t]['True'].sort_values(by='logfoldchanges', ascending=False).to_excel(writer, sheet_name = str(t))

cycling_genes_u = []
lfc = 0.25
for t in dexs_cycling_u:
    print(t)
    genes = dexs_cycling_u[t]['True'][dexs_cycling_u[t]['True']['logfoldchanges']>lfc].index
    for idx in genes:
        for g in idx.rsplit('-'):
            c = [x for x in adatas[t].var.index if g in x and '-' not in x]
            if len(c) > 0:
                cycling_genes_u += c
            else:
                cycling_genes_u += [idx]
                print(g,c,idx)

cycling_genes_u = sorted(set(cycling_genes_u))
cgudf = pd.DataFrame(pd.Series(cycling_genes_u, name = 'S-phase_unspliced-genes'))
cgudf.to_csv(outdir + '/Sphase_unsplicedgenes.tsv', sep = '\t')

pd.DataFrame(adatas['merged'].obs['S-phase']).to_csv(outdir + '/Sphase_cells.tsv', sep = '\t')

sys.exit()


