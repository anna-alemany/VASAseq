#!/usr/bin/env python
# coding: utf-8
import sys, os
import numpy as np
import pandas as pd
from pandas.io.parsers import read_csv
from collections import Counter
import scanpy as sc
import plot_aautils as plaa
import sc_aautils as scaa
import pickle
import itertools as it
from scipy.stats import fisher_exact as fisher_exact
import glob
import matplotlib.pyplot as plt
from pandarallel import pandarallel
import bbknn
import scipy
from scipy.stats import norm
#import plotly

plt.ion()

try:
    timepoint = sys.argv[1]
except:
    sys.exit("Please, provide timepoint (E65,E75,E85)")

outdir = '../'+timepoint+'/res_scanpy_merged_2021_final'
os.system('mkdir -p '+outdir)

# folders with input h5ad objects
scanpyPJall = 'res_scanpy_all_2021_v2'
scanpyVASAall = 'res_scanpy_all_2021_noHistones'
scanpyPJhigh = 'res_scanpy_high_2021_v2'
scanpyVASAhigh = 'res_scanpy_high_2021_noHistones'

if timepoint == 'E65':
    n_pca = 20
else:
    n_pca = 40

# some functions
def plUmapsDiscrete(adata, metric = 'correlation', columns = ['leiden_correlation'], size = 5, factorsize = 3):
    if metric == None:
        udf = pd.DataFrame(adata.obsm['X_umap'], columns = ['u1','u2'], index = adata.obs.index)        
    else:
        udf = pd.DataFrame(adata.obsm['X_umap_'+metric], columns = ['u1','u2'], index = adata.obs.index)
    if len(columns) == 1:
        column = columns[0]
        fig, axs = plaa.template_plot(figsize = (factorsize*1.6,factorsize))
        if 'nan' in set(adata.obs[column]):
            cl = 'nan'
            cells = adata.obs[adata.obs[column]==cl].index
            axs.scatter(udf.loc[cells,'u1'], udf.loc[cells,'u2'], s = size, label = cl, c = 'gray')
        if 'not-assigned' in set(adata.obs[column]):
            cl = 'not-assigned'
            cells = adata.obs[adata.obs[column]==cl].index
            axs.scatter(udf.loc[cells,'u1'], udf.loc[cells,'u2'], s = size, label = cl, c = 'silver')
        for i, cl in enumerate(sorted(set(adata.obs[column]))):
            cells = adata.obs[adata.obs[column]==cl].index
            if cl not in ['nan','not-assigned'] : 
                axs.scatter(udf.loc[cells,'u1'], udf.loc[cells,'u2'], s = size, label = cl, c = plaa.colors()[i])                              
        axs.grid(False); axs.set_xticks([]); axs.set_yticks([])
    else:
        fig, axs = plaa.template_plot(ncols=len(columns), figsize=(factorsize*1.6*len(columns), factorsize))
        for ax, column in zip(axs, columns):
            if 'nan' in set(adata.obs[column]):
                cl = 'nan'
                cells = adata.obs[adata.obs[column]==cl].index
                ax.scatter(udf.loc[cells,'u1'], udf.loc[cells,'u2'], s = size, label = cl, c = 'gray')
            if 'not-assigned' in set(adata.obs[column]):
                cl = 'not-assigned'
                cells = adata.obs[adata.obs[column]==cl].index
                ax.scatter(udf.loc[cells,'u1'], udf.loc[cells,'u2'], s = size, label = cl, c = 'silver')
            for i, cl in enumerate(sorted(set(adata.obs[column]))):
                cells = adata.obs[adata.obs[column]==cl].index
                if cl not in ['nan','not-assigned']: 
                    ax.scatter(udf.loc[cells,'u1'], udf.loc[cells,'u2'], s = size, label = cl, c = plaa.colors()[i])
            ax.grid(False); ax.set_xticks([]); ax.set_yticks([])
    return fig, axs


# # 10X libraries
xasu = scaa.read_h5ad_zipped('../../pijuansala/'+timepoint+'/'+scanpyPJall+'/PJ_su_All_'+timepoint+'.h5ad')
x3s = scaa.read_h5ad_zipped('../../pijuansala/'+timepoint+'/'+scanpyPJhigh+'/PJ_s_All_'+timepoint+'.h5ad')
x3su = scaa.read_h5ad_zipped('../../pijuansala/'+timepoint+'/'+scanpyPJhigh+'/PJ_su_All_'+timepoint+'.h5ad')

print(x3s.shape, x3su.shape)
print('x3s:', Counter(x3s.obs['PJ_celltype']))
print('x3su:', Counter(x3su.obs['PJ_celltype']))

# ## VASA libraries
vasu = scaa.read_h5ad_zipped('../'+timepoint+'/'+scanpyVASAall+'/VASA_su_All_'+timepoint+'.h5ad')
v3s = scaa.read_h5ad_zipped('../'+timepoint+'/'+scanpyVASAhigh+'/VASA_s_All_'+timepoint+'.h5ad')
v3su = scaa.read_h5ad_zipped('../'+timepoint+'/'+scanpyVASAhigh+'/VASA_su_All_'+timepoint+'.h5ad')

print(v3s.shape, v3su.shape)

# ## PJ Celltypes
fig, axs = plaa.template_plot(ncols = 2, figsize = (2*3*1.6, 3))
for ax, adata, t in zip(axs, [xasu, x3su], ['all reads', '3 reads']):
    u10 = pd.DataFrame(adata.obsm['X_umap'], index = adata.obs.index, columns = ['u1','u2'])
    for i, cl in enumerate(sorted(set(xasu.obs['PJ_celltype']))):
        c10 = adata.obs[adata.obs['PJ_celltype'] == cl].index
        ax.scatter(u10.loc[c10, 'u1'], u10.loc[c10,'u2'], c = plaa.colors()[i], s = 3, label = '-'.join([str(i+1),cl]))
        ax.text(u10.loc[c10, 'u1'].mean(), u10.loc[c10,'u2'].mean(), str(i+1), va = 'center', ha = 'center', fontsize = 8)
    ax.grid(False); ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(t)
lgn = axs[0].legend(loc = 2, bbox_to_anchor = (0,0), ncol = 3)
for handle in lgn.legendHandles:
    handle.set_sizes([20.0])
plt.savefig(outdir + '/umapS_allreads_10x_PJcelltypes_'+timepoint+'.pdf', bbox_inches= 'tight')
plt.close()

# select clusters from vasu and xasu to transfer
clalgo = 'leiden_clean_manhattan'
xasu.obs['asu_leiden'] = xasu.obs[clalgo]
vasu.obs['asu_leiden'] = vasu.obs[clalgo]

x3s.obs['asu_leiden'] = [xasu.obs.loc[idx,'asu_leiden'] if idx in xasu.obs.index else -1 for idx in x3s.obs.index]
x3su.obs['asu_leiden'] = [xasu.obs.loc[idx,'asu_leiden'] if idx in xasu.obs.index else -1 for idx in x3su.obs.index]
v3s.obs['asu_leiden'] = [vasu.obs.loc[idx,'asu_leiden'] if idx in vasu.obs.index else -1 for idx in v3s.obs.index]
v3su.obs['asu_leiden'] = [vasu.obs.loc[idx,'asu_leiden'] if idx in vasu.obs.index else -1 for idx in v3su.obs.index]

x3s = x3s[x3s.obs['asu_leiden']!=-1,:]
x3su = x3su[x3su.obs['asu_leiden']!=-1,:]
v3s = v3s[v3s.obs['asu_leiden']!=-1,:]
v3su = v3su[v3su.obs['asu_leiden']!=-1,:]

for adata in [x3s, x3su, v3s, v3su]:
    adata.obs['asu_leiden'] = adata.obs['asu_leiden'].astype('category')

for adatas, technology in zip([(xasu,x3su), (vasu, v3su)], ['10x','vasa']):
    fig, axs = plaa.template_plot(ncols = 2, figsize = (2*3*1.6, 3))
    for ax, adata, t in zip(axs, adatas, ['all reads', '3 reads']):
        u10 = pd.DataFrame(adata.obsm['X_umap'], index = adata.obs.index, columns = ['u1','u2'])
        for i, cl in enumerate(sorted(set(adata.obs['asu_leiden']))):
            c10 = adata.obs[adata.obs['asu_leiden'] == cl].index
            ax.scatter(u10.loc[c10, 'u1'], u10.loc[c10,'u2'], c = plaa.colors()[i], s = 3, label = cl)
            ax.text(u10.loc[c10, 'u1'].mean(), u10.loc[c10,'u2'].mean(), i, va = 'center', ha = 'center', fontsize = 8)
        ax.grid(False); ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(t)
    lgn = axs[0].legend(loc = 2, bbox_to_anchor = (0,0), ncol = 3)
    for handle in lgn.legendHandles:
        handle.set_sizes([20.0])
    plt.savefig(outdir + '/umapS_'+technology+'_asuLeiden_'+timepoint+'.pdf', bbox_inches= 'tight')
    plt.close()

# # Merge 10X with VASA
# ## Extract first raw data (no filtered by noisy genes etc)
rv3s = scaa.getRawData(v3s)
rv3su = scaa.getRawData(v3su)
rvasu = scaa.getRawData(vasu)

rx3s = scaa.getRawData(x3s)
rx3su = scaa.getRawData(x3su)
rxasu = scaa.getRawData(xasu)

# ## Concatenate vasa and 10x raw data
m3s = scaa.concatSCobj(rv3s, rx3s, ['vasa','10x'])
m3su = scaa.concatSCobj(rv3su, rx3su, ['vasa','10x'])
masu = scaa.concatSCobj(rvasu, rxasu, ['vasa','10x'])

# ## Filter merged datasets
def basicFilter(mdf, f = True):
    mdf.var['n_counts'] = (np.exp(mdf.X)-1).sum(axis=0)
    mdf.var['n_cells'] = (mdf.X>0).sum(axis=0)
    mdf = mdf[:, ['-' not in idx for idx in mdf.var.index]]
    if f: 
        mdf = mdf[:, mdf.var['n_cells']>10] # keeping genes present in at least 10 cells
    return mdf

m3s = basicFilter(m3s, f = False)
m3su = basicFilter(m3su, f = False)
masu = basicFilter(masu, f = False)

# ## annotate genes when detected in vasa or in 10x
def geneInBatch(mdata, b1data, b2data, name_batch1, name_batch2):
    """Annotates batch in gene metadata"""
    mdata.var['in-'+name_batch1] = [idx in b1data.var.index for idx in mdata.var.index]
    mdata.var['in-'+name_batch2] = [idx in b2data.var.index for idx in mdata.var.index]
    mdata.var['n_counts_'+name_batch1] = mdata[mdata.obs['batch']==name_batch1,].X.sum(axis=0)
    mdata.var['n_counts_'+name_batch2] = mdata[mdata.obs['batch']==name_batch2,].X.sum(axis=0)
    mdata.var['n_cells_'+name_batch1] = (mdata[mdata.obs['batch']==name_batch1,].X>0).sum(axis=0)
    mdata.var['n_cells_'+name_batch2] = (mdata[mdata.obs['batch']==name_batch2,].X>0).sum(axis=0)
    return mdata

m3s = geneInBatch(m3s, rv3s, rx3s, 'vasa','10x')
m3su = geneInBatch(m3su, rv3su, rx3su, 'vasa','10x')
masu = geneInBatch(masu, rvasu, rxasu, 'vasa','10x')

# ## pie charts showing common genes
fig, axs = plt.subplots(ncols = 3, figsize = (3*3*1.5,3))
i = 0
for ax, adata, t in zip(axs, [m3s, m3su, masu], ['3'', only spliced', '3'', spliced and unspliced', 'all, spliced and unspliced']):
    cnt = pd.crosstab(adata.var['in-10x'], adata.var['in-vasa']).values.reshape(4)
    ax.pie(cnt[1:], labels = ['only vasa','only 10x', 'both'], autopct='%1.1f%%')
    ax.set_title(t)
fig.savefig(outdir + '/piechart_detectedGenes.pdf', bbox_inches = 'tight')
#os.system('cp '+outdir + '/piechart_detectedGenes.pdf /Users/anna/Dropbox/vasa/fig3/panel_b/piechart_detectedGenes_'+timepoint+'.pdf')
plt.close()

fig, axs = plt.subplots(ncols = 3, figsize = (3*3*1.5,3))
i = 0
for ax, adata, t in zip(axs, [m3s, m3su, masu], ['3'', only spliced', '3'', spliced and unspliced', 'all, spliced and unspliced']):
    cnt = [adata.var[(adata.var['in-10x'])&(np.invert(adata.var['in-vasa']))][['n_counts']].sum().values[0], adata.var[(adata.var['in-vasa'])&(np.invert(adata.var['in-10x']))][['n_counts']].sum().values[0], adata.var[(adata.var['in-10x'])&(adata.var['in-vasa'])][['n_counts']].sum().values[0]]
    ax.pie(cnt, labels = ['only 10x','only vasa', 'both'], autopct='%1.1f%%')
    ax.set_title(t)
fig.savefig(outdir + '/piechart_detectedReads.pdf', bbox_inches = 'tight')
#os.system('cp '+outdir + '/piechart_detectedReads.pdf /Users/anna/Dropbox/vasa/fig3/panel_b/piechart_detectedReads_'+timepoint+'.pdf')
plt.close()


m3su.var[['in-10x','in-vasa','n_counts','n_cells', 'n_counts_vasa','n_cells_vasa','n_counts_10x','n_cells_10x']].to_csv(outdir + '/piechart_detectedGenes_3su_'+timepoint+'.tsv', sep = '\t')
masu.var[['in-10x','in-vasa', 'n_counts','n_cells','n_counts_vasa','n_cells_vasa','n_counts_10x','n_cells_10x']].to_csv(outdir + '/piechart_detectedGenes_asu_'+timepoint+'.tsv', sep = '\t')

# I do the transfer using only 3' reads
m3s.raw = m3s
m3su.raw = m3su

m3s = m3s[:, (m3s.var['in-vasa'])&(m3s.var['in-10x'])]
m3su = m3su[:, (m3su.var['in-vasa'])&(m3su.var['in-10x'])]

sc.pp.highly_variable_genes(m3s, min_mean=0.0125, max_mean=5, min_disp=0.5, n_bins = 30)
sc.pp.highly_variable_genes(m3su, min_mean=0.0125, max_mean=5, min_disp=0.5, n_bins = 30)

m3s = m3s[:, m3s.var['highly_variable']] 
m3su = m3su[:, m3su.var['highly_variable']]

# ## combat to correct for batch effects between the 10X and VASA
sc.pp.combat(m3s, key='batch')
sc.pp.combat(m3su, key='batch')

# ## PCA
def pcaplease(mdf):
    mdf.raw = mdf
    sc.pp.regress_out(mdf, ['batch'])
    sc.pp.regress_out(mdf, ['n_counts'])
    sc.pp.scale(mdf, max_value=10)
    sc.tl.pca(mdf, svd_solver='arpack', n_comps = min(150, mdf.shape[1]-50))
    return mdf

m3s = pcaplease(m3s)
m3su = pcaplease(m3su)

def plotPCA(mas):
    pdf = pd.DataFrame(mas.obsm['X_pca'], index = mas.obs.index)
    fig, axs = plaa.template_plot(ncols = 2, figsize=(2*3*1.6,3))
    for b in set(mas.obs['batch']):
        cells = mas.obs[mas.obs['batch']==b].index
        axs[0].scatter(pdf.loc[cells,0], pdf.loc[cells,1], s = 1, label = b)
    axs[0].set_xlabel('PCA 1'); axs[0].set_ylabel('PCA 2')
    axs[1].scatter(range(1,len(mas.uns['pca']['variance_ratio'])+1), mas.uns['pca']['variance_ratio'], s = 10)
    axs[1].set_yscale('log')
    axs[1].set_xlim(0,50)
    axs[1].set_xlabel('PCA component'); axs[1].set_ylabel('variance')
    return fig, axs

for adata, t in zip([m3s, m3su],['m3s','m3su']):
    fig, ax = plotPCA(adata)
    plt.savefig(outdir + '/pca_'+t+'.pdf', bbox_inches = 'tight')
    plt.close()

# this generates UMAPs with bbknn distances
def pleaseUmap(adata, n_pca = n_pca, metric = 'manhattan'):
    #sc.pp.neighbors(mdf, n_neighbors = new_k, n_pcs = n_pca, random_state = 1971723, metric = metric)
    bbknn.bbknn(adata, batch_key = 'batch', n_pcs = n_pca, metric = metric)#, neighbors_within_batch = 10)
    sc.tl.umap(adata, n_components = 2, random_state = 921225, min_dist = 0.5, spread = 1)
    adata.obsm['X_umap_'+metric] = adata.obsm['X_umap']
    return adata

m3s = pleaseUmap(m3s)
m3su = pleaseUmap(m3su)

# select cells that are also present in vasu and xasu, and transfer those clusters

# interbatch distances
def gaussian(x, mu, s2):
    return (1./np.sqrt(2*np.pi*s2))*np.exp(-0.5*(x-mu)**2/s2)

fig, aaxs = plaa.template_plot(ncols = 4, nrows = 2, figsize = (4*3*1.6, 3*2))
for axs, t, adata in zip(aaxs, ['m3s','m3su'], [m3s, m3su]):
    print(t)
    clalgo = 'asu_leiden'
    distdf = scaa.get_distDF(adata, metric= '')
    adata.uns['knn_th'] = {'vasa=>10x': {}, '10x=>vasa': {}}
    for ax, (b1, b2) in zip(axs, it.product(['vasa','10x'],['vasa','10x'])):
        print(b1, b2)
        subdistdf = distdf.loc[adata.obs['batch']==b1, adata.obs['batch']==b2]
        ds = np.array(subdistdf)[np.where(subdistdf>0)] # This array will have all annotated neigbors
#        ds = np.array(((subdistdf==0)*1e3 + subdistdf).min(axis=1)) # this will have only first nearest neighbor for each cell
        m = np.mean(ds); std = np.std(ds)
        ax.hist(ds, bins = 100, label = '-'.join([b1,b2]))
        ax.axvline(m, c = 'r'); ax.text(m, 50, 'mean')
        ax.axvline(m+1.5*std, c = 'r', ls = '--'); ax.text(m+1.5*std, 40, 'mean+1.5*std')
        ax.legend()
        ax.set_xlabel('distances')
        ax.set_ylabel('frequencies')

        N = len(set(adata.obs.loc[subdistdf.index, clalgo]))
        fig2, ax2 = plaa.template_plot(nrows = N, figsize = (3*1.6,3*N))
        for aaxx22, ct in zip(ax2, sorted(set(adata.obs.loc[subdistdf.index, clalgo]))):
            cells = subdistdf.index[adata.obs.loc[subdistdf.index, clalgo] == ct]
            dmins = np.array(subdistdf.loc[cells])[np.where(subdistdf.loc[cells]>0)] # this will have all annotated neigbors
#            dmins = subdistdf.loc[cells][subdistdf.loc[cells]>0].min(axis=1) # this will have only first nearest neighbor
            aaxx22.hist(dmins, bins = 100, label = ct, density = True)
            aaxx22.legend()
            xra = np.linspace(dmins.min()*0.8, dmins.max()*1.1, 100)
            aaxx22.plot(xra, gaussian(xra, dmins.mean(), dmins.var()), c = 'red')
            print(ct, b1, b2, dmins.mean())
            if b2 != b1:
                if len(dmins) > 1:
                    adata.uns['knn_th'][b1 + '=>' + b2][ct] = {'mu': dmins.mean(), 's2': dmins.var()}
                else:
                    adata.uns['knn_th'][b1 + '=>' + b2][ct] = {'mu': dmins.mean(), 's2': 1e-5}
        fig2.savefig(outdir + '/bbknn_distances_'+t+'_'+'-'.join([b1, b2])+'.pdf', bbox_inches = 'tight')

    axs[0].set_title(t)
    adata.obsp['clean_bbknn_distances'] = scipy.sparse.csr_matrix(distdf.values)
fig.savefig(outdir + '/bbknn_distances.pdf', bbox_inches = 'tight')
plt.close('all')

# cell type assignment
m3s_dist = scaa.get_distDF(m3s, metric = '')
m3su_dist = scaa.get_distDF(m3su, metric = '')

def find_compatible_cluster(dist_v, adata, ref_batch = '10x', mapped_batch = 'vasa'):
    clalgo = 'asu_leiden'
    cell = dist_v.name
    k = ref_batch + '=>' + mapped_batch
    if adata.obs.loc[cell,'batch'] == ref_batch:
        ct = adata.obs.loc[cell,clalgo]
        score = 1
    else:
        knn = dist_v[dist_v>0].index # find kNN (should be 5 + itself)
        knn_ref = knn[adata.obs.loc[knn,'batch']==ref_batch] # select kNN that belong to reference batch (should be 3)
        if len(knn_ref) > 0:
            d = pd.DataFrame(dist_v.loc[knn_ref])
            d.columns = ['dist']
            d['celltype'] = [adata.obs.loc[idx,clalgo] for idx in d.index]
            d['score'] = d.apply(lambda x: 1-norm(loc = adata.uns['knn_th'][k][x['celltype']]['mu'], scale = np.sqrt(adata.uns['knn_th'][k][x['celltype']]['s2'])).cdf(x['dist']), axis = 1)
            ct = d[d['score']==d['score'].max()]['celltype'].values[0]
            score = d[d['score']==d['score'].max()]['score'].values[0]
        else: 
            ct = 'not-assigned'
            score = 0
    return pd.Series({'celltype': ct, 'score': score})

def new_function(dist_v, adata):
    cell = dist_v.name
    clalgo = 'asu_leiden'
    celltype = adata.obs.loc[cell,'asu_leiden']
    knn = dist_v[dist_v>0].index
    knn_ref = knn[adata.obs.loc[knn,'batch']!=adata.obs.loc[cell,'batch']]
    d = pd.DataFrame(dist_v.loc[knn_ref])
    d.columns = ['dist']
    d['celltype'] = [adata.obs.loc[idx,clalgo] for idx in d.index]
    k = adata.obs.loc[d.index[0],'batch'] + '=>' + adata.obs.loc[cell,'batch'] 
    d['score'] = d.apply(lambda x: 1-norm(loc = adata.uns['knn_th'][k][x['celltype']]['mu'], scale = np.sqrt(adata.uns['knn_th'][k][x['celltype']]['s2'])).cdf(x['dist']), axis = 1)
    ct = d[d['score']==d['score'].max()]['celltype'].values[0]
    score = d[d['score']==d['score'].max()]['score'].values[0]
    dist =  d[d['score']==d['score'].max()]['dist'].values[0]
    return pd.Series({'celltype': celltype, 'assigned': ct, 'score': score, 'batch': adata.obs.loc[cell,'batch'], 'dist': dist})

pandarallel.initialize(nb_workers = 8)
m3su_celltypes_assignment = m3su_dist.parallel_apply(lambda x: new_function(x, m3su), axis = 1)

pandarallel.initialize(nb_workers = 8)
m3s_celltypes_assignment = m3s_dist.parallel_apply(lambda x: new_function(x, m3s), axis = 1)

m3su.obs['celltype_assigned'] = m3su_celltypes_assignment.loc[m3su.obs.index, 'assigned']
m3su.obs['assignment_score'] =  m3su_celltypes_assignment.loc[m3su.obs.index, 'score']
m3su.obs['assignment_distance'] =  m3su_celltypes_assignment.loc[m3su.obs.index, 'dist']

m3s.obs['celltype_assigned'] = m3s_celltypes_assignment.loc[m3s.obs.index, 'assigned']
m3s.obs['assignment_score'] =  m3s_celltypes_assignment.loc[m3s.obs.index, 'score']
m3s.obs['assignment_distance'] =  m3s_celltypes_assignment.loc[m3s.obs.index, 'dist']

# ## 
vasu.obs['x2v_assigned_3su'] = [m3su_celltypes_assignment.loc[idx,'assigned'] if idx in m3su_celltypes_assignment.index else -1 for idx in vasu.obs.index] 
vasu.obs['x2v_score_3su'] = [m3su_celltypes_assignment.loc[idx,'score'] if idx in m3su_celltypes_assignment.index else 0 for idx in vasu.obs.index]
vasu.obs['x2v_dist_3su'] = [m3su_celltypes_assignment.loc[idx,'dist'] if idx in m3su_celltypes_assignment.index else 0 for idx in vasu.obs.index]

xasu.obs['v2x_assigned_3su'] = [m3su_celltypes_assignment.loc[idx,'assigned'] if idx in m3su_celltypes_assignment.index else -1 for idx in xasu.obs.index] 
xasu.obs['v2x_score_3su'] = [m3su_celltypes_assignment.loc[idx,'score'] if idx in m3su_celltypes_assignment.index else 0 for idx in xasu.obs.index]
xasu.obs['v2x_dist_3su'] = [m3su_celltypes_assignment.loc[idx,'dist'] if idx in m3su_celltypes_assignment.index else 0 for idx in xasu.obs.index]

# plot results
fig, ax = plaa.template_plot()
ax.hist(vasu.obs['x2v_score_3su'], bins = 100, alpha = 0.5, label = 'vasa cells with 10x label', density = True)
ax.hist(xasu.obs['v2x_score_3su'], bins = 100, alpha = 0.5, label = '10x cells with vasa labels', density = True)
ax.set_xlabel('score')
ax.set_ylabel('density')
ax.legend()
plt.savefig(outdir + '/hist_scores.pdf', bbox_inches= 'tight')

# labels from 10x to vasa
mudf = pd.DataFrame(m3su.obsm['X_umap'], columns = ['u1','u2'], index = m3su.obs.index)
xudf = pd.DataFrame(xasu.obsm['X_umap'], columns = ['u1','u2'], index = xasu.obs.index)
vudf = pd.DataFrame(vasu.obsm['X_umap'], columns = ['u1','u2'], index = vasu.obs.index)

N = len(set(xasu.obs['asu_leiden']))
fig, maxs = plaa.template_plot(ncols = 4, nrows = N, figsize = (4*3*1.6, N*3))
for axs, cl in zip(maxs, sorted(set(xasu.obs['asu_leiden']))):
    for ax, udf in zip(axs, [mudf, mudf, xudf, vudf]):
        ax.scatter(udf['u1'], udf['u2'], s = 5, c = 'silver', rasterized = True)
        ax.grid(False); ax.set_xticks([]); ax.set_yticks([])
    xcells = m3su.obs[(m3su.obs['batch']=='10x')&(m3su.obs['asu_leiden']==cl)].index
    axs[0].scatter(mudf.loc[xcells,'u1'], mudf.loc[xcells,'u2'], s = 5, c = 'red', label = cl); axs[0].legend()
    axs[2].scatter(xudf.loc[xcells,'u1'], xudf.loc[xcells,'u2'], s = 5, c = 'red')
    vcells = m3su.obs[(m3su.obs['batch']=='vasa')&(m3su.obs['celltype_assigned']==cl)].index
    vcells = m3su.obs.loc[vcells].sort_values(by='assignment_score').index
    axs[1].scatter(mudf.loc[vcells,'u1'], mudf.loc[vcells,'u2'], s = 5, c = m3su.obs.loc[vcells,'assignment_score'], cmap = 'gnuplot', vmin = 0, vmax = 1)
    im = axs[3].scatter(vudf.loc[vcells,'u1'], vudf.loc[vcells,'u2'], s = 5, c = vasu.obs.loc[vcells,'x2v_score_3su'], cmap = 'gnuplot', vmin = 0, vmax = 1)
    plt.colorbar(im, ax = axs[3], label = 'score')
maxs[0][0].set_title('merged, 10x cells'); maxs[0][1].set_title('merged, vasa cells'); maxs[0][2].set_title('10x'); maxs[0][3].set_title('vasa')
plt.savefig(outdir + '/umapS_m3su_10xtoVASAassignments.pdf', bbox_inches = 'tight')

# labels from vasa to 10x
N = len(set(vasu.obs['asu_leiden']))
fig, maxs = plaa.template_plot(ncols = 4, nrows = N, figsize = (4*3*1.6, N*3))
for axs, cl in zip(maxs, sorted(set(vasu.obs['asu_leiden']))):
    for ax, udf in zip(axs, [mudf, mudf, vudf, xudf]):
        ax.scatter(udf['u1'], udf['u2'], s = 5, c = 'silver', rasterized = True)
        ax.grid(False); ax.set_xticks([]); ax.set_yticks([])
    vcells = m3su.obs[(m3su.obs['batch']=='vasa')&(m3su.obs['asu_leiden']==cl)].index
    axs[0].scatter(mudf.loc[vcells,'u1'], mudf.loc[vcells,'u2'], s = 5, c = 'red', label = cl); axs[0].legend()
    axs[2].scatter(vudf.loc[vcells,'u1'], vudf.loc[vcells,'u2'], s = 5, c = 'red')
    xcells = m3su.obs[(m3su.obs['batch']=='10x')&(m3su.obs['celltype_assigned']==cl)].index
    xcells = m3su.obs.loc[xcells].sort_values(by='assignment_score').index
    axs[1].scatter(mudf.loc[xcells,'u1'], mudf.loc[xcells,'u2'], s = 5, c = m3su.obs.loc[xcells,'assignment_score'], cmap = 'gnuplot', vmin = 0, vmax = 1)
    im = axs[3].scatter(xudf.loc[xcells,'u1'], xudf.loc[xcells,'u2'], s = 5, c = xasu.obs.loc[xcells,'v2x_score_3su'], cmap = 'gnuplot', vmin = 0, vmax = 1)
    plt.colorbar(im, ax = axs[3], label = 'score')
maxs[0][0].set_title('merged, vasa cells'); maxs[0][1].set_title('merged, 10x cells'); maxs[0][2].set_title('vasa'); maxs[0][3].set_title('10x')
plt.savefig(outdir + '/umapS_m3su_VASAto10xassignments.pdf', bbox_inches = 'tight')


# output dataframes
uv = pd.DataFrame(vasu.obsm['X_umap_clean_manhattan'], columns = ['umap1', 'umap2'], index = vasu.obs.index).merge(vasu.obs, how = 'inner', left_index = True, right_index = True)
ux = pd.DataFrame(xasu.obsm['X_umap_clean_manhattan'], columns = ['umap1', 'umap2'], index = xasu.obs.index).merge(xasu.obs, how = 'inner', left_index = True, right_index = True)

uv.to_csv(outdir + '/VASA_su_All_UmapMeta_'+timepoint+'.tsv', sep = '\t')
#uv.to_csv('/Users/anna/Dropbox/vasa/fig3/panel_c/umaps/VASA_su_All_UmapMeta_'+timepoint+'.tsv', sep = '\t')
ux.to_csv(outdir + '/PJ_su_All_UmapMeta_'+timepoint+'.tsv', sep = '\t')
#ux.to_csv('/Users/anna/Dropbox/vasa/fig3/panel_c/umaps/PJ_su_All_UmapMeta_'+timepoint+'.tsv', sep = '\t')

sys.exit()

# aluvial plots?
def genSankey(df,cat_cols=[],value_cols='',title='Sankey Diagram'):
    """from https://medium.com/kenlok/how-to-create-sankey-diagrams-from-dataframes-in-python-e221c1b4d6b0"""
    # maximum of 6 value cols -> 6 colors
    colorPalette = ['#4B8BBE','#306998','#FFE873','#FFD43B','#646464']
    labelList = []; colorNumList = []
    for catCol in cat_cols:
        labelListTemp =  list(set(df[catCol].values))
        colorNumList.append(len(labelListTemp))
        labelList = labelList + labelListTemp

    # remove duplicates from labelList
    labelList = list(dict.fromkeys(labelList))

    # define colors based on number of levels
    colorList = []
    for idx, colorNum in enumerate(colorNumList):
        colorList = colorList + [colorPalette[idx]]*colorNum

    # transform df into a source-target pair
    for i in range(len(cat_cols)-1):
        if i==0:
            sourceTargetDf = df[[cat_cols[i],cat_cols[i+1],value_cols]]
            sourceTargetDf.columns = ['source','target','count']
        else:
            tempDf = df[[cat_cols[i],cat_cols[i+1],value_cols]]
            tempDf.columns = ['source','target','count']
            sourceTargetDf = pd.concat([sourceTargetDf,tempDf])
        sourceTargetDf = sourceTargetDf.groupby(['source','target']).agg({'count':'sum'}).reset_index()

    # add index for source-target pair
    sourceTargetDf['sourceID'] = sourceTargetDf['source'].apply(lambda x: labelList.index(x))
    sourceTargetDf['targetID'] = sourceTargetDf['target'].apply(lambda x: labelList.index(x))
    
    # creating the sankey diagram
    data = dict(
        type='sankey',
        node = dict(
          pad = 15, thickness = 20, line = dict(color = "black", width = 0.5),
          label = labelList, color = colorList
        ),
        link = dict(
          source = sourceTargetDf['sourceID'], target = sourceTargetDf['targetID'], value = sourceTargetDf['count']
        )
      )
    
    layout =  dict(
        title = title, font = dict(size = 10)
    )
       
    fig = dict(data=[data], layout=layout)
    return fig

#df = pd.DataFrame({'v1': ['A']*8+['B']*3, 'v2': ['AP']*3+['AC']*3+['AB']*2+['BE','BR','BA'], 'v3': ['APP','APE', 'APA', 'ACT', 'ACC', 'ACE', 'ABL', 'ABO', 'BET',' BRE', 'BAK'], 'count': [5,2,3,8,2,10,1,3,4,6,3]})
#fig = genSankey(df,cat_cols=['v1','v2','v3'],value_cols='count',title='Word Etymology')
#plotly.offline.plot(fig, validate=False)

c1 = 'leiden_clean_manhattan'; c2 = 'PJ_celltype'
df = ux[[c1, c2]]
for c in df.columns:
    df[c] = np.array(df[c])
cdf = pd.DataFrame({'count': df.groupby(list(df.columns)).agg(len)})
for i in range(len(cdf.index.names)):
    cdf[cdf.index.names[i]] = np.array([idx[i] for idx in cdf.index])
fig = genSankey(cdf,cat_cols=[c1, c2],value_cols='count',title='')
plotly.offline.plot(fig, output_type = 'file', image = 'png', image_filename = outdir + '/aluvial_xasu_Leiden2celltypes_' + timepoint )

c1 = 'leiden_clean_manhattan_3'; c2 = 'PJ_celltype'
df = ux[[c1, c2]]
df = df[df[c1] != -1]
for c in df.columns:
    df[c] = np.array(df[c])
cdf = pd.DataFrame({'count': df.groupby(list(df.columns)).agg(len)})
for i in range(len(cdf.index.names)):
    cdf[cdf.index.names[i]] = np.array([idx[i] for idx in cdf.index])
fig = genSankey(cdf,cat_cols=[c1, c2],value_cols='count',title='')
plotly.offline.plot(fig, output_type = 'file', image = 'png', image_filename = outdir + '/aluvial_x3su_Leiden2celltypes_' + timepoint )

c1 = 'x2v_assigned_3su'; c2 = 'leiden_clean_manhattan'
df = uv[[c1, c2, 'x2v_score_3su']]
df = df[df['x2v_score_3su']>0.2][[c1,c2]]
for c in df.columns:
    df[c] = np.array(df[c])
for l, c in zip(['10x', 'vasa'], [c1, c2]): 
    df[c] = ['-'.join([str(int(x)),l]) for x in df[c]]
cdf = pd.DataFrame({'count': df.groupby(list(df.columns)).agg(len)})
for i in range(len(cdf.index.names)):
    cdf[cdf.index.names[i]] = np.array([idx[i] for idx in cdf.index])
cdf = cdf[cdf['count'] != 1]
fig = genSankey(cdf,cat_cols=[c1, c2],value_cols='count',title='')
plotly.offline.plot(fig, output_type = 'file', image = 'png', image_filename = outdir + '/aluvial_x2v_HQ_assignment_' + timepoint )

sys.exit()

# Differential gene expression 
def difGeneExpr(adata, clusters =  ['leiden_manhattan','louvain_manhattan']):
    dex = {}
    for cluster in clusters:
        print('*'+cluster+'*')
        try:
            adata.obs[cluster] = adata.obs[cluster].astype(str)
            sc.tl.rank_genes_groups(adata, cluster, method='t-test', n_genes=200)
            dex[cluster] = {}
            for cl in set(adata.obs[cluster]):
                try:
                    keys = ['names', 'logfoldchanges', 'pvals', 'pvals_adj', 'scores']
                    dex[cluster][cl] = pd.DataFrame({k: [adata.uns['rank_genes_groups'][k][i][cl] for i in range(200)] for k in keys})
                    dex[cluster][cl] = dex[cluster][cl].set_index('names')
                except:
                    continue
        except:
            continue
    return dex

cnt = Counter(vasu.obs['leiden_assigned_3su'])
for cl in cnt:
    if cnt[cl]<=1:
        for cell in vasu.obs[vasu.obs['leiden_assigned_3su']==cl].index:
            vasu.obs.loc[cell,'leiden_assigned_3su'] = -1
dex = difGeneExpr(vasu, clusters = ['leiden_assigned_3su'])

# 
udf = pd.DataFrame(vasu.obsm['X_umap'], columns = ['u1','u2'], index = vasu.obs.index).merge(vasu.obs, how = 'inner', left_index = True, right_index = True)
udf.to_csv('/Users/anna/Dropbox/vasa/umaps/E65/umap-metadata-'+timepoint+'.tsv', sep = '\t')


