#!/usr/bin/env python
# coding: utf-8
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
    timepoint = sys.argv[1]
except:
    sys.exit("Please, give timepoint (E65, E75, E85, E95)")

metric = 'manhattan'

if timepoint == 'E65':
    sample_out = []
    reads_su_th = {'All': [10**3.5,10**6]}; reads_s_th = {'All': [10**3.5,10**5]}
    genes_su_th = {}; genes_s_th = {'All': [2500,10000]}
    fracs_su_th = {'ProteinCoding': [.7,.9], 'lncRNA': [.007,.025], 'smallRNA': [.05,0.35]}; fracs_s_th = {}
    n_pca = 20
    resolution_s = {'All': 0.7, 'ProteinCoding': 0.65, 'lncRNA': 0.75, 'smallRNA': 0.5, 'TF': 0.7, 'Cofactor': 0.7}
    resolution_su = {'All': 0.8, 'ProteinCoding': 0.65, 'lncRNA': 0.75, 'smallRNA': 0.5, 'TF': 0.7, 'Cofactor': 0.7}
    min_disp = 1
elif timepoint == 'E75':
    sample_out = []
    reads_su_th = {'All': [10**3,10**6]}; reads_s_th = {'All': [10**3,10**5]}
    genes_su_th = {}; genes_s_th = {'All': [1200,8000]}
    fracs_su_th = {'ProteinCoding': [.7,.9], 'lncRNA': [.007,.02], 'smallRNA': [.05,0.3]}; fracs_s_th = {}
    n_pca = 40 
    resolution_s = {'All': 1, 'ProteinCoding': 1, 'lncRNA': 1, 'smallRNA': 0.75, 'TF': 1, 'Cofactor': 0.75}
    resolution_su =  {'All': 1, 'ProteinCoding': 1, 'lncRNA': 1, 'smallRNA': 0.75, 'TF': 1, 'Cofactor': 0.75}
    min_disp = 0.5
elif timepoint == 'E85': 
    sample_out = ['E8.5-8_i21']
    reads_su_th = {'All': [10**3,10**6]}; reads_s_th = {'All': [10**3, 10**5]}
    genes_su_th = {}; genes_s_th = {'All': [1250,8000]}
    fracs_su_th = {'ProteinCoding': [.7,.9], 'lncRNA': [.007,.02], 'smallRNA': [.07,0.3]}; fracs_s_th = {}
    resolution_s = {'All': 1, 'ProteinCoding': 1, 'lncRNA': 1, 'smallRNA': 0.75, 'TF': 1, 'Cofactor': 0.75}
    resolution_su =  {'All': 1, 'ProteinCoding': 1, 'lncRNA': 1, 'smallRNA': 0.75, 'TF': 1, 'Cofactor': 0.75}
    n_pca = 40 
    min_disp = 0.5
elif timepoint == 'E95':
    sample_out = []
    reads_su_th = {'All': [10**3,10**6]}; reads_s_th = {'All': [10**3, 10**5]}
    genes_su_th = {}; genes_s_th = {'All': [1000,7000]}
    fracs_su_th = {'ProteinCoding': [.7,.9], 'lncRNA': [.007,.02], 'smallRNA': [.07,0.3]}; fracs_s_th = {}
    n_pca = 40
    resolution_s = {'All': 1, 'ProteinCoding': 1, 'lncRNA': 1, 'smallRNA': 0.75, 'TF': 1, 'Cofactor': 0.75}
    resolution_su = {'All': 1, 'ProteinCoding': 1, 'lncRNA': 1, 'smallRNA': 0.75, 'TF': 1, 'Cofactor': 0.75}
    min_disp = 0.5

outdir = '../'+timepoint+'/res_scanpy_high_2021_v2'
os.system('mkdir -p '+outdir)

# ## Read raw scanpy objects
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
sc.settings.figdir = outdir
sc.settings.set_figure_params(dpi=80)

bios = ['All','ProteinCoding','lncRNA','smallRNA','TF','Cofactor']
folder = 'res_scanpy_high_rawQC_2021'
adata_s = {bio: sc.read_h5ad('../' + timepoint + '/' + folder + '/VASA_s_raw_' + bio + '.h5ad') for bio in bios}
adata_u = {bio: sc.read_h5ad('../' + timepoint + '/' + folder + '/VASA_u_raw_' + bio + '.h5ad') for bio in bios}
adata_su = {bio: sc.read_h5ad('../' + timepoint + '/' + folder + '/VASA_su_raw_' + bio + '.h5ad') for bio in bios}

# ## Filter data
def filter_BiotypeThresholds(adata, reads_th, genes_th, fracs_th, sample_out):
    for n in adata:
        print(n)
        print(adata[n].shape)
        adata[n] = adata[n][[s not in sample_out for s in adata[n].obs['sample']],:]
        if n in reads_th:
            th_min = reads_th[n][0]; th_max = reads_th[n][1]
            adata[n] = adata[n][(th_min < adata[n].obs['n_counts'])&(adata[n].obs['n_counts'] < th_max),:]
        if n in genes_th:
            th_min = genes_th[n][0]; th_max = genes_th[n][1]
            adata[n] = adata[n][(th_min < adata[n].obs['n_genes'])&(adata[n].obs['n_genes'] < th_max),:]
        if n in fracs_th: 
            th_min = fracs_th[n][0]; th_max = fracs_th[n][1]
            adata[n] = adata[n][(th_min < adata[n].obs['n_counts_'+n]/adata[n].obs['n_counts'])&(adata[n].obs['n_counts_'+n]/adata[n].obs['n_counts'] < th_max),:]
        print(adata[n].shape)
        print('--')
    return adata

data_su = filter_BiotypeThresholds(adata_su, reads_su_th, genes_su_th, fracs_su_th, sample_out)
adata_s = filter_BiotypeThresholds(adata_s, reads_s_th, genes_s_th, fracs_s_th, sample_out)

good_cells = []
for n in adata_s:
    good_cells += list(adata_s[n].obs.index)
    good_cells += list(adata_su[n].obs.index)

cnt_good_cells = Counter(good_cells)
selected_cells = [cell for cell in cnt_good_cells if cnt_good_cells[cell]==len(adata_su)+len(adata_s)]
print('Cell survival: '+ str(len(cnt_good_cells))+', '+str(len(selected_cells)))

for n in adata_s:
    adata_s[n] = adata_s[n][[idx in selected_cells for idx in adata_s[n].obs.index],:]
    adata_u[n] = adata_u[n][[idx in selected_cells for idx in adata_u[n].obs.index],:]
    adata_su[n] = adata_su[n][[idx in selected_cells for idx in adata_su[n].obs.index],:]

del adata_u['smallRNA']

# ## Normalization of scanpy objects and gene selection
def normalize_data(adata):
    for n in adata:
        sc.pp.normalize_per_cell(adata[n], counts_per_cell_after=1e4) 
        sc.pp.log1p(adata[n]) 
        adata[n].raw = adata[n]
    return adata

adata_s = normalize_data(adata_s)
adata_u = normalize_data(adata_u)
adata_su = normalize_data(adata_su)

# ## Neighbors check
def get_distDF(adata, metric = 'manhattan'):
    if len(metric) > 0:
        dist_df = pd.DataFrame.sparse.from_spmatrix(adata.obsp[metric + '_distances'])
    else:
        dist_df = pd.DataFrame.sparse.from_spmatrix(adata.obsp['distances'])
    dist_df.columns = adata.obs.index
    dist_df.index = adata.obs.index
    return dist_df

def getRawData(scdata):
    rsc = sc.AnnData(scdata.raw.X)
    rsc.var = scdata.raw.var
    rsc.obs = scdata.obs
    return rsc

# correct distances and filter out cells, but only for su[All]
adata = adata_su['All'].copy()
n_iter = 0; N = 1000
while N > 0:
    n_iter += 1; Ns = Counter()
    adata.var['n_cells'] = (adata.X>0).sum(axis=0)
    adata = adata[:,adata.var['n_cells']>2]
    sc.pp.highly_variable_genes(adata, min_mean = 0.0125, max_mean = 5, min_disp = min_disp, n_bins = 30)
    adata = adata[:, adata.var['highly_variable']]
    sc.pp.regress_out(adata, ['n_counts'])
    sc.pp.scale(adata, max_value=10)
    ncomps = min(150, abs(adata.shape[1]-50))
    sc.tl.pca(adata, svd_solver='arpack', n_comps = ncomps)
    sc.pp.neighbors(adata, n_neighbors = 10, n_pcs = n_pca, random_state = 1971723, metric = metric, key_added = metric) # regression
    distdf = get_distDF(adata, metric = metric)
    ds = np.array(distdf)[np.where(distdf>0)]
    ds_th = np.mean(ds)+1.5*np.std(ds)
    distdf =  ((distdf < ds_th)*distdf).astype(float)
    adata.obsp['clean_'+metric+'_distances'] = scipy.sparse.csr_matrix(distdf.values)
    adata.obs['kNN_'+metric] = (distdf>0).sum(axis=1)
    Ns = list(adata.obs[adata.obs['kNN_'+metric]==0].index)

    adata = getRawData(adata)
    adata = adata[ [idx not in Ns for idx in adata.obs.index] ,:]
    adata.raw = adata

    N = len(Ns)     
    print('=====>', n_iter, N, adata.shape)

final_cells = adata.obs.index

#
def correct_kNNdistances(adata, cells):
    for n in adata: 
        adata[n] = getRawData(adata[n])
        adata[n] = adata[n][ [idx in cells for idx in adata[n].obs.index] ,:]
        adata[n].raw = adata[n]
        adata[n].var['n_cells'] = (adata[n].X>0).sum(axis=0)
        sc.pp.highly_variable_genes(adata[n], min_mean = 0.0125, max_mean = 5, min_disp = min_disp, n_bins = 30)
        adata[n] = adata[n][:, adata[n].var['highly_variable']]
        sc.pp.regress_out(adata[n], ['n_counts'])
        sc.pp.scale(adata[n], max_value=10)
        ncomps = min(150, abs(adata[n].shape[1]-50))
        sc.tl.pca(adata[n], svd_solver='arpack', n_comps = ncomps)
        sc.pp.neighbors(adata[n], n_neighbors = 10, n_pcs = n_pca, random_state = 1971723, metric = metric, key_added = metric)
        distdf = get_distDF(adata[n], metric = metric)
        ds = np.array(distdf)[np.where(distdf>0)]
        ds_th = np.mean(ds)+1.5*np.std(ds)
        distdf =  ((distdf < ds_th)*distdf).astype(float)
        adata[n].obsp['clean_'+metric+'_distances'] = scipy.sparse.csr_matrix(distdf.values)
        adata[n].obs['kNN_'+metric] = (distdf>0).sum(axis=1)
        adata[n].obs['pointers_'+metric] = (distdf>0).sum()
    return adata

adata_s = correct_kNNdistances(adata_s, final_cells)
adata_su = correct_kNNdistances(adata_su, final_cells)
adata_u = correct_kNNdistances(adata_u, final_cells)

# plot PCAs
def scatter_pca_varianceRatio(adata, outputname):
    fig, maxs = plaa.template_plot(ncols = 2, nrows = len(adata), figsize=(2*3*1.6,3*len(adata)))
    for axs, n in zip(maxs, adata):
        df = adata[n]

        ax = axs[0]
        pca_df = pd.DataFrame(df.obsm['X_pca'])
        pca_df.index = df.obs.index
        for i, p in enumerate(sorted(set(df.obs['sample']))):
            cells = df.obs[df.obs['sample']==p].index
            ax.scatter(pca_df.loc[cells, 0], pca_df.loc[cells,1], s = 1, c = plaa.colors()[i], label = p)
        ax.set_xlabel(f"PCA1 ({df.uns['pca']['variance_ratio'][0]*100:4.2f}%)")
        ax.set_ylabel(f"PCA2 ({df.uns['pca']['variance_ratio'][1]*100:4.2f}%)")
        ax.legend(loc = 2, bbox_to_anchor =(0,-0.2))
        ax.grid(False); ax.grid(c = 'silver', lw = 0.5, ls = '--')
        ax.set_title(n)

        ax = axs[1]
        ax.scatter(range(len(df.uns['pca']['variance_ratio'])), df.uns['pca']['variance_ratio'], s = 5)
        ax.set_xlim(-1,80)
        ax.set_xlabel('ranking'); ax.set_ylabel('variance ratio')
        ax.grid(False); ax.grid(c = 'silver', lw = 0.5, ls = '--')
        ax.set_yscale('log')
        ax.set_title(n)
    fig.savefig(outdir + '/'+ outputname, bbox_inches = 'tight')
    plt.close()

scatter_pca_varianceRatio(adata_s, 'pca_s.pdf')
scatter_pca_varianceRatio(adata_u, 'pca_u.pdf')
scatter_pca_varianceRatio(adata_su, 'pca_su.pdf')

# UMAPs
def findSigma(cell, new_dist_v, rho, new_k):
    if len(new_dist_v) <= 1:
        sigma = 1
    else:
        eps = min([1e-2, (new_dist_v.sort_values()[-1]-new_dist_v.sort_values()[-2])/100])
        sigma1 = np.array(new_dist_v).min()/100;  sigma2 = new_k
        exps = np.exp(-(new_dist_v-rho))
        log2k = np.log2(new_k)
        f1 = (exps**(1./sigma1)).sum()-log2k; f2 = (exps**(1./sigma2)).sum()-log2k
        i = 0; j = 0
        while abs(sigma2-sigma1)>eps:
            if f1*f2 > 0 and f1 > 0:
                sigma1 /= 10; f1 = (exps**(1./sigma1)).sum()-log2k
                j += 1
                if j == 5:
                    sigma1 = sigma2 = sigma = 1e-5
            elif f2*f2 > 0 and f2 < 0:
                sigma2 *= 10; f2 = (exps**(1./sigma2)).sum()-log2k
                i += 1
                if i == 5:
                    sigma1 = sigma2 = sigma = 1000
            elif f1*f2 < 0:
                sigma = 0.5*(sigma1+sigma2)
                f =  (exps**(1./sigma)).sum()-log2k
                if f > 0:
                    sigma2 = sigma; f2 = f
                elif f < 0:
                    sigma1 = sigma; f1 = f
    return sigma

def directed_weigths(new_dist_cell, rho_cell, sigma_cell):
    knn = new_dist_cell[new_dist_cell>0].index
    cv = pd.Series(0, index = new_dist_cell.index)
    for c in knn:
        cv.loc[c] = np.exp(-(new_dist_cell.loc[c]-rho_cell)/sigma_cell)
    return cv

# clean one last time for disconnected cells
def remove_disconnectedCells(adata):
    for n in adata:
        rm = adata[n].obs[(adata[n].obs['kNN_'+metric]<=1)&(adata[n].obs['pointers_'+metric]<=1)].index
        adata[n] = adata[n][[idx not in rm for idx in adata[n].obs.index],:]
        distdf = get_distDF(adata[n], metric = 'clean_manhattan')
        adata[n].obs['kNN_'+metric+'_2'] = (distdf>0).sum(axis=1)
        adata[n].obs['pointers_'+metric+'_2'] = (distdf>0).sum()
        rm  = adata[n].obs[(adata[n].obs['kNN_'+metric+'_2']==0)&(adata[n].obs['pointers_'+metric+'_2']==0)].index
        adata[n] = adata[n][[idx not in rm for idx in adata[n].obs.index],:]
    return adata

adata_s = remove_disconnectedCells(adata_s)
adata_u = remove_disconnectedCells(adata_u)
adata_su = remove_disconnectedCells(adata_su)

# obtain umaps
def easy_umap(adata):
    for n in adata:
        sc.tl.umap(adata[n], n_components = 2, random_state = 921225, min_dist = 0.5, spread = 1, neighbors_key = metric)
        adata[n].obsm['X_umap_'+metric] = adata[n].obsm['X_umap']
    return adata

adata_s = easy_umap(adata_s)
adata_u = easy_umap(adata_u)
adata_su = easy_umap(adata_su)

def filtered_umap(adata):
    distdf = get_distDF(adata, metric = 'clean_manhattan')

    pandarallel.initialize(nb_workers = 8)
    rho = distdf.parallel_apply(lambda x: x[x>0].min(), axis = 1)

    pandarallel.initialize(nb_workers = 8)
    sigma = distdf.parallel_apply(lambda x: findSigma(x.name, distdf.loc[x.name,distdf.columns[distdf.loc[x.name]>0]], rho.loc[x.name], 10), axis = 1)

    pandarallel.initialize(nb_workers = 8)
    wdf = distdf.parallel_apply(lambda x: directed_weigths(x, rho.loc[x.name], sigma.loc[x.name]), axis = 1) 

    bdf = wdf + wdf.T - wdf * wdf.T
    adata.obsp['clean_'+metric+'_connectivities'] = scipy.sparse.csr_matrix(bdf.values)

    adata.uns['clean_'+metric] = {'connectivities_key': 'clean_'+metric+'_connectivities',
                       'distances_key': 'clean_'+metric+'_distances',
                       'params': {'method': 'mine','metric': metric, 'n_neighbors': 10, 'n_pca': 20}}

    sc.tl.umap(adata, n_components = 2, random_state = 921225, min_dist = 0.5, spread = 1, neighbors_key = 'clean_'+metric)
    adata.obsm['X_umap_clean_'+metric] = adata.obsm['X_umap']
    return adata

adata_su['All'] = filtered_umap(adata_su['All'])

### 
def plotUMAPS(adatas = adata_s, col = 'sample', metric = 'clean_manhattan'):
    if len(adatas) > 1:
        fig, axs = plaa.template_plot(ncols = len(adatas), nrows = 1, figsize=(len(adatas)*3*1.6,3))
        for n, (ax, k) in enumerate(zip(axs, adatas)):
            df = adatas[k]
            umap_df = pd.DataFrame(df.obsm['X_umap_'+metric], columns = ['u1','u2'], index = df.obs.index)
            for i, p in enumerate(sorted(set(df.obs[col]))):
                cells = df.obs[df.obs[col]==p].index
                ax.scatter(umap_df.loc[cells, 'u1'], umap_df.loc[cells,'u2'], s = 5, c = plaa.colors()[i], label = p)
                ax.grid(False); ax.set_title(k); ax.set_xticks([]); ax.set_yticks([])
    else:
        fig, ax = plaa.template_plot(ncols = 1, nrows = 1, figsize=(3*1.6,3))
        for k in adatas:
            df = adatas[k]
            umap_df = pd.DataFrame(df.obsm['X_umap_'+metric], columns = ['u1','u2'], index = df.obs.index)
            for i, p in enumerate(sorted(set(df.obs[col]))):
                cells = df.obs[df.obs[col]==p].index
                ax.scatter(umap_df.loc[cells, 'u1'], umap_df.loc[cells,'u2'], s = 5, c = plaa.colors()[i], label = p)
                ax.grid(False); ax.set_title(k); ax.set_xticks([]); ax.set_yticks([])
        axs = ax
    return fig, axs

fig, axs = plotUMAPS(adata_s, col = 'sample', metric = metric)
plt.savefig(outdir + '/umap_s_'+'-'.join(adata_s)+'_'+metric+'_samples.pdf', bbox_inches = 'tight')

fig, axs = plotUMAPS(adata_u, col = 'sample', metric = metric)
plt.savefig(outdir + '/umap_u_'+'-'.join(adata_u)+'_'+metric+'_samples.pdf', bbox_inches = 'tight')

fig, axs = plotUMAPS(adata_su, col = 'sample', metric = metric)
plt.savefig(outdir + '/umap_su_'+'-'.join(adata_su)+'_'+metric+'_samples.pdf', bbox_inches = 'tight')

fig, ax = plotUMAPS({'All': adata_su['All']}, col = 'sample', metric = 'clean_'+metric)
plt.savefig(outdir + '/umap_su_All_clean'+metric+'_samples.pdf', bbox_inches = 'tight')

plt.close('all')

# # Clustering
def leidenCluster(adata, leiden_resolution, metric):
    for n in adata:
        sc.tl.leiden(adata[n], resolution = leiden_resolution[n], neighbors_key = metric, key_added = 'leiden_'+metric)
        adata[n].obs['leiden_'+metric] = adata[n].obs['leiden_'+metric].astype(int).astype('category')
    return adata

adata_s =  leidenCluster(adata_s, resolution_s, metric)
adata_su =  leidenCluster(adata_su, resolution_su, metric)
adata_u =  leidenCluster(adata_u, resolution_s, metric)

adata =  leidenCluster({'All': adata_su['All']}, {'All': resolution_su['All']}, 'clean_'+metric)
adata_su['All'] = adata['All']

fig, axs = plotUMAPS(adata_s, col = 'leiden_'+metric, metric = metric)
plt.savefig(outdir + '/umap_s_'+'-'.join(adata_s)+'_leiden'+metric+'_samples.pdf', bbox_inches = 'tight')

fig, axs = plotUMAPS(adata_u, col = 'leiden_'+metric, metric = metric)
plt.savefig(outdir + '/umap_u_'+'-'.join(adata_u)+'_leiden'+metric+'_samples.pdf', bbox_inches = 'tight')

fig, axs = plotUMAPS(adata_su, col = 'leiden_'+metric, metric = metric)
plt.savefig(outdir + '/umap_su_'+'-'.join(adata_su)+'_leiden'+metric+'_samples.pdf', bbox_inches = 'tight')

fig, ax = plotUMAPS({'All': adata_su['All']}, col = 'leiden_clean_'+metric, metric = 'clean_'+metric)
plt.savefig(outdir + '/umap_su_All_leidenclean'+metric+'_samples.pdf', bbox_inches = 'tight')

# # Differential gene expression analysis
def difGeneExpr(adata, clusters =  ['leiden_manhattan','louvain_manhattan']):
    dex = {}
    for cluster in clusters:
        print('*'+cluster+'*')
        try:
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

dex_s = {bio: difGeneExpr(adata_s[bio], clusters = ['leiden_' + metric]) for bio in adata_s}
dex_u = {bio: difGeneExpr(adata_u[bio], clusters = ['leiden_' + metric]) for bio in adata_u}
dex_su = {bio: difGeneExpr(adata_su[bio], clusters = ['leiden_' + metric]) for bio in adata_su}
dex_su_clean = difGeneExpr(adata_su['All'], clusters = ['leiden_clean_' + metric])
dex_su['All-clean'] = dex_su_clean

# print data
for n in adata_s:
    adata_s[n].write(outdir + '/VASA_s_' + n + '_' + timepoint + '.h5ad')
    pickle.dump(dex_s[n], open(outdir + '/dex_VASA_s_' + n + '_' + timepoint + '.pickle', 'wb'))

for n in adata_u:
    adata_u[n].write(outdir + '/VASA_u_' + n + '_' + timepoint + '.h5ad')
    pickle.dump(dex_u[n], open(outdir + '/dex_VASA_u_' + n + '_' + timepoint + '.pickle', 'wb'))

for n in adata_su:
    adata_su[n].write(outdir + '/VASA_su_' + n + '_' + timepoint + '.h5ad')
    pickle.dump(dex_su[n], open(outdir + '/dex_VASA_su_' + n + '_' + timepoint + '.pickle', 'wb'))






