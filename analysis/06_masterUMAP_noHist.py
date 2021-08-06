#!/usr/bin/env python
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
import scipy
from pandarallel import pandarallel
import multiprocessing
import scipy.spatial.distance as scdist
from scipy.optimize import root
from umap.umap_ import smooth_knn_dist # returns sigmas and ros

plt.ion()

old_k = 80
new_k = 30
Ncells = 50000
onlySpliced = True
markers_or_hgv = 'hgv'

# option 1: seqPast
metric = 'manhattan'
tNN = {'E6.5': ['E6.5'],
        'E7.5': ['E7.5','E6.5'],
        'E8.5': ['E8.5','E7.5'],
        'E9.5': ['E9.5','E8.5']}
output = 'test_seqPast_'+metric+'_s_noHist_RegresHistOut'
# option 2: seqFuture
#metric = 'manhattan'
#tNN = {'E6.5': ['E6.5','E7.5'], 
#        'E7.5': ['E7.5','E8.5'], 
#        'E8.5': ['E8.5','E9.5'], 
#        'E9.5': ['E9.5']}
#output = 'test_seqFuture_manhattan_onlySpliced'
# option 3: seqSeq
#metric = 'manhattan'
#tNN = {'E6.5': ['E6.5','E7.5'],
#        'E7.5': ['E7.5','E8.5','E6.5'],
#        'E8.5': ['E8.5','E9.5','E7.5'],
#        'E9.5': ['E9.5','E8.5']}
#output = 'test_seqSeq_manhattan_onlySpliced'

# input meta
#scanpyDir = 'res_scanpy_all_2021_v2/'
scanpyDir = 'res_scanpy_all_histone'
inputScanpyE65 = '../E65/' + scanpyDir
inputScanpyE75 = '../E75/' + scanpyDir
inputScanpyE85 = '../E85/' + scanpyDir
inputScanpyE95 = '../E95/' + scanpyDir

# timepoints
timepoints = ['E65','E75','E85','E95']
scanpyinputs = [inputScanpyE65,inputScanpyE75,inputScanpyE85,inputScanpyE95]

# ============================================================================ #

# Identify markers
def findMarker(dex, clusterAlgorithm = 'leiden', lfc_th = 1, pv_th = 1e-2):
    markers = []
    for cl in dex[clusterAlgorithm]:
        m = dex[clusterAlgorithm][cl][(dex[clusterAlgorithm][cl]['logfoldchanges'] > lfc_th) & (dex[clusterAlgorithm][cl]['pvals_adj'] < pv_th)].index
        markers += list(m)
    markers = list(set(markers))
    return markers

#markers = {'manhattan': {}, 'correlation': {}}
#for timepoint, p2scanpy in zip(timepoints, scanpyinputs):
#    print('markers for ', timepoint)
#    if onlySpliced:
#        os.system('gunzip '+p2scanpy + '/dex_VASA_s_All_'+timepoint+'.pickle.gz')
#        dexfile = glob.glob(p2scanpy + '/dex_VASA_s_All_'+timepoint+'.pickle')
#        dex = pickle.load(open(dexfile[0], 'rb'))
#        os.system('gzip '+p2scanpy + '/dex_VASA_s_All_'+timepoint+'.pickle')
#    else:
#        os.system('gunzip '+p2scanpy + '/dex_VASA_su_All_'+timepoint+'.pickle.gz')
#        dexfile = glob.glob(p2scanpy + '/dex_VASA_su_All_'+timepoint+'.pickle')
#        dex = pickle.load(open(dexfile[0], 'rb'))
#        os.system('gzip '+p2scanpy + '/dex_VASA_su_All_'+timepoint+'.pickle')
#    markers['manhattan'][timepoint] = findMarker(dex, clusterAlgorithm = 'leiden_manhattan', lfc_th = 1.1, pv_th = 1e-3)
#    markers['correlation'][timepoint] = findMarker(dex, clusterAlgorithm = 'leiden_correlation')

# merge all markers into a list
#print([len(markers[metric][t]) for t in markers[metric]])
all_markers = []
#for t in markers[metric]:
#    all_markers += markers[metric][t]
#all_markers = list(set(all_markers))
#print(len(all_markers))

# get raw data
def getRawData(scdata):
    rsc = sc.AnnData(scdata.raw.X)
    rsc.var = scdata.raw.var
    rsc.obs = scdata.obs
    return rsc

def getPCA(scdata):
    pca = pd.DataFrame(scdata.varm['PCs'], index = scdata.var.index)
    n_pca = scdata.uns[metric]['params']['n_pcs']
    pca = pca[pca.columns[:n_pca]]
    return pca

def get_distDF(adata):
    dist_df = pd.DataFrame.sparse.from_spmatrix(adata.obsp['distances'])
    dist_df.columns = adata.obs.index
    dist_df.index = adata.obs.index
    return dist_df

# read all data ad save it in a scanpy object
adatas = {}
pcas = {}
unss = {}
for timepoint, p2scanpy in zip(timepoints, scanpyinputs):
    print('reading counttables for ',timepoint)
    if onlySpliced:
        if os.path.isfile(p2scanpy + '/VASA_s_All_'+timepoint+'.h5ad.gz'):
            os.system('gunzip '+p2scanpy + '/VASA_s_All_'+timepoint+'.h5ad.gz')
        adata = sc.read_h5ad(p2scanpy + '/VASA_s_All_'+timepoint+'.h5ad')
        os.system('gzip '+p2scanpy + '/VASA_s_All_'+timepoint+'.h5ad')
    else:
        if os.path.isfile(p2scanpy + '/VASA_su_All_'+timepoint+'.h5ad.gz'):
            os.system('gunzip '+p2scanpy + '/VASA_su_All_'+timepoint+'.h5ad.gz')
        adata = sc.read_h5ad(p2scanpy + '/VASA_su_All_'+timepoint+'.h5ad')
        os.system('gzip '+p2scanpy + '/VASA_su_All_'+timepoint+'.h5ad')
    pcas[timepoint] = getPCA(adata)
    unss[timepoint] = adata.uns
    adata = getRawData(adata)
    print(timepoint, adata.shape)
    adatas[timepoint] = adata[[i<Ncells for i,c in enumerate(adata.obs.index)],:]

# merge all timeponints
def concatSCobj(sc1, sc2, names = ['1','2']):
    msc = sc1.concatenate(sc2, join = 'outer', batch_categories = None, index_unique = None)
    msc.X = np.nan_to_num(msc.X, nan = 0.0, posinf = 0, neginf = 0)
    return msc

mdf = adatas[timepoints[0]]
for i, t in enumerate(timepoints[1:]):
    mdf = concatSCobj(mdf, adatas[t], names = [timepoints[i], timepoints[i+1]])

print('merged', mdf.shape)
mdf.obs['batch'] = [idx.rsplit('-')[0] for idx in mdf.obs.index]

# Analyze
mdf.var['n_counts'] = mdf.X.sum(axis=0)
mdf.var['n_cells'] = (mdf.X>0).sum(axis=0)
#mdf = mdf[:, mdf.var['n_cells']>2]

sc.pp.highly_variable_genes(mdf, min_mean=0.0125, max_mean=5, min_disp=0.5, n_bins = 30)
mdf.var['marker'] = [g in all_markers for g in mdf.var.index]
 
for x in ['id','polyA','reg', 'ubiotype','mbiotype','biotype','n_counts','n_cells']:
    if x+'-1' in mdf.var.columns:
        mdf.var[x] = mdf.var[x+'-1']
        idcols = [c for c in mdf.var.columns if x+'-' in c]
        for c in idcols:
           del mdf.var[c]
     
mdf.raw = mdf
if markers_or_hgv == 'hgv':
    mdf = mdf[:, mdf.var['highly_variable']]
elif markers_or_hgv == 'markers':
    mdf = mdf[:, mdf.var['marker']]

# remove histone genes
mdf = mdf[:, np.invert(mdf.var['histone'])]

# s-phase cells
sc.pp.regress_out(mdf, ['n_counts'])
sc.pp.regress_out(mdf, ['S-phase'])
sc.pp.scale(mdf, max_value=10)
sc.tl.pca(mdf, svd_solver='arpack', n_comps = min(150, mdf.shape[1]-50))

#sc.pl.pca(mdf, color = 'batch')
#sc.pl.pca_variance_ratio(mdf, log = True)

n_pca = 40

# normal umap
sc.pp.neighbors(mdf, n_neighbors = new_k, n_pcs = n_pca, random_state = 1971723, metric = metric)
sc.tl.umap(mdf, n_components = 2, random_state = 921225, min_dist = 0.5, spread = 1)
mdf.obsm['X_umap_standard'] = mdf.obsm['X_umap']

#sc.pl.umap(mdf, color = 'batch')#, save = 'standard.pdf')

### some checks about distances and pca's
#pca_genes = pd.DataFrame(mdf.varm['PCs'], index = mdf.var.index)
#pca_cells = pd.DataFrame(mdf.obsm['X_pca'], index = mdf.obs.index)
#dist_df = get_distDF(mdf)

#c0 = dist_df.index[0]
#c1 = dist_df.columns[dist_df.loc[c0]>0][0]
#dist_df.loc[c0,c1]

#expr_v = pd.DataFrame(mdf[[idx in [c0,c1] for idx in mdf.obs.index],:].X, columns = mdf.var.index, index = [idx for idx in mdf.obs.index if idx in [c0,c1]])

#expr2pca = pd.DataFrame(np.dot(expr_v, pca_genes), index = expr_v.index)
#pca_cells.loc[[c0,c1]] # expr2pca and pca_cells are practically the same

#np.abs(pca_cells.loc[c0]-pca_cells.loc[c1]).sum()
#np.abs(pca_cells.loc[c0,:40]-pca_cells.loc[c1,:40]).sum() # <= this is it
#np.abs(expr2pca.loc[c0,:40]-expr2pca.loc[c1,:40]).sum() # <= this is also it
#import scipy.spatial.distance as scdist
#scdist.pdist(expr2pca[range(40)], metric = 'cityblock')
###

# new way to compute distances
def wagner_distances(x): # cell, mdf = mdf, pcas = pcas, old_k = old_k, new_k = new_k):
    cell, mdf, pcas, old_k, new_k = x
    tref = mdf.obs.loc[cell,'batch']
    pca_gene = pcas[tref.replace('.','')]
    cells = mdf.obs.index[[x in tNN[tref] for x in mdf.obs['batch']]]
    expr_v = pd.DataFrame(mdf.raw[[idx in cells for idx in mdf.obs.index],:].X, columns = mdf.raw.var.index, index = [idx for idx in mdf.obs.index if idx in cells]) 
    expr_v = expr_v[pca_gene.index]
    expr2pca = pd.DataFrame(np.dot(expr_v, pca_gene), index = expr_v.index)
    print(cell)
    vd = expr2pca.apply(lambda x: scdist.pdist([x, expr2pca.loc[cell]], metric = 'cityblock')[0], axis = 1)
    vd = vd.sort_values()[:old_k]
    d_min = vd[vd>0].min()#; d_std = vd[vd>0].std()
    vd = vd[vd < 3*d_min]
    vd = vd.sort_values()[:new_k]
#    distances = pd.Series({x: vd.loc[x] if x in vd.index else 0 for x in mdf.obs.index})
    return vd # distances

def wagner_distances_v2(tref, mdf = mdf, pcas = pcas, old_k = old_k, new_k = new_k):
    pca_gene = pcas[tref.replace('.','')]
    cells = mdf.obs.index[[x in tNN[tref] for x in mdf.obs['batch']]]
    expr_v = pd.DataFrame(mdf.raw[[idx in cells for idx in mdf.obs.index],:].X, columns = mdf.raw.var.index, index = [idx for idx in mdf.obs.index if idx in cells])
    expr_v = expr_v[pca_gene.index]
    expr2pca = pd.DataFrame(np.dot(expr_v, pca_gene), index = expr_v.index)
    scsub = sc.AnnData(expr_v)
    scsub.obs = mdf.obs.loc[scsub.obs.index]
    scsub.varm['PCs'] = np.array(pca_gene)
    scsub.obsm['X_pca'] = np.array(expr2pca)
    scsub.uns = unss[tref.replace('.','')]
    sc.pp.neighbors(scsub, n_neighbors = new_k, n_pcs = len(pca_gene.columns)-1, random_state = 1971723, metric = metric)
    distdf = get_distDF(scsub)
    cells = mdf.obs.index[[x in tref for x in mdf.obs['batch']]]
    return distdf.loc[cells]


#pandarallel.initialize(nb_workers = 8)
#new_dist_df = mdf.obs.head().parallel_apply(lambda x: wagner_distances(x.name), axis = 1) # does not work
#new_dist_df = mdf.obs.apply(lambda x: wagner_distances(x.name), axis = 1) # super slow
#pool = multiprocessing.Pool(8) #threads
#for idx, (vd) in enumerate(pool.imap_unordered(wagner_distances,[(idx, mdf, pcas, old_k, new_k) for idx in mdf.obs.index[:10]])): # did not work
#    print(vd[0])
#new_dist_df = mdf.obs.apply(lambda x: wagner_distances((x.name, mdf, pcas, old_k, new_k)), axis = 1) # at the end this is decent, but still slow

distdfs = {t: wagner_distances_v2(t) for t in tNN}
new_dist_df = pd.concat([distdfs['E6.5'], distdfs['E7.5'], distdfs['E8.5'], distdfs['E9.5']], sort = False)
new_dist_df = new_dist_df.fillna(0)

mdf.obsp['distances'] = scipy.sparse.csr_matrix(new_dist_df.values)

###############
# another aproach (no pca space correction)
# new connections
#sc.pp.neighbors(mdf, n_neighbors = old_k, n_pcs = n_pca, random_state = 1971723, metric = metric)

# now, modify distances
#dist_df = get_distDF(mdf)

#def reshape_distance_v(dist_df_v, new_k = new_k, tNN = tNN):
#    cell = dist_df_v.name
#    print(cell)
#    dist_df_v = pd.Series(np.array(dist_df_v), index = dist_df_v.index)
#    t_ref = cell.rsplit('-')[0] # timepoint of reference cell
#    vd = dist_df_v[dist_df_v>0].copy()# vector kNN distances
#    vt = pd.Series({idx: idx.rsplit('-')[0] for idx in vd.index}) # vector times
#    vtc = vt.apply(lambda x: x in tNN[t_ref])
#    # 1) clean neighbors that are not connected in a sequential temporal appraoch
#    vd = vd*vtc
#    # 2) keep distances if dij < 3*min(dik) 
#    if (vd>0).sum() > 1: 
#        d_min = vd[vd>0].min() 
#        d_std = vd[vd>0].std()
##       vd = (vd <= 3*d_min)*vd
#        vd = (vd <= d_min + 3.5*d_std)*vd
#    # 3) some extra conditions that for now i dont implement
#    # 4) reduce number of neighbors to new_k
#    if (vd>0).sum() > 0:
#        d_th = (vd[vd>0]).sort_values().iloc[min([new_k-1, sum(vd>0)-1])]
#        vd = (vd<=d_th)*vd
#    # reassign to dist_df_v
#    for cell in vd.index:
#        dist_df_v.loc[cell] = vd.loc[cell]
#    return dist_df_v

#new_dist_df = dist_df.apply(lambda x: reshape_distance_v(x), axis = 1)

#from pandarallel import pandarallel
#pandarallel.initialize(nb_workers = 8)
#new_dist_df = dist_df.parallel_apply(lambda x: reshape_distance_v(x), axis = 1)

#mdf.obsp['distances'] = scipy.sparse.csr_matrix(new_dist_df.values)

new_dist_df = new_dist_df.astype(float)

pandarallel.initialize(nb_workers = 8)
rho = new_dist_df.parallel_apply(lambda x: x[x>0].min(), axis = 1) # i can also do this in the dictionary...

###
def findSigma_v1(cell, new_dist_v, rho, new_k):
    if len(new_dist_v) <= 1:
        sigma = 1
    else:
        log2k = np.log2(new_k)
        eps = 1e-3 
        sigma1 = np.array(new_dist_v).min()/100; sigma2 = new_k
        exps = np.exp(-(new_dist_v-rho))
        f1 = (exps**(1./sigma1)).sum()-log2k; f2 = (exps**(1./sigma2)).sum()-log2k
        i = 0
        while abs(sigma2-sigma1)>eps:
            if f1*f2 > 0 and f1 > 0:
                print('ei')
                break
            elif f2*f2 > 0 and f2 < 0:
                sigma2 *= 10
                i += 1
                if i == 5:
                    sigma1 = sigma2 = sigma = 1000
                    break
            elif f1*f2 < 0:
                sigma = 0.5*(sigma1+sigma2)
                f =  (exps**(1./sigma)).sum()-log2k
                if f > 0:
                    sigma2 = sigma; f2 = f
                elif f < 0:
                    sigma1 = sigma; f1 = f
    return sigma

def fun_sigma(x, a, b, c): 
    return np.exp(-(a-b)/x).sum()-c

def findSigma_v2(cell, new_dist_v, rho, new_k):
    if len(new_dist_v) <= 1: 
        sigma = 1
    else: 
        sigma = root(lambda x: fun_sigma(x, new_dist_v, rho, np.log2(new_k)), x0 = 0.5*new_k).x[0]
    return sigma

def findSigma_v3(cell, new_dist_v, rho, new_k):
    if len(new_dist_v) <= 1:
        sigma = 1
    else:
        n_iter = 100
        eps = 1e-3 # min([1e-2, (new_dist_df.loc[cell,knn].sort_values()[-1]-new_dist_df.loc[cell,knn].sort_values()[-2])/100])
        lo = 0; hi = new_k
        sigma = 0.5*(lo + hi)
        log2k = np.log2(new_k)
        exps = np.exp(-(new_dist_v-rho))
        f = (exps**(1./sigma)).sum()-log2k
        for n in range(n_iter):
#        while np.abs(f) > eps:
            if np.abs(f) < eps:
                break
            if f > 0:
                hi = sigma
            else:
                lo = sigma
            sigma = 0.5*(lo + hi)
            f = (exps**(1./sigma)).sum()-log2k
    return sigma

import time
cell = np.random.choice(new_dist_df.index)
print('speed test: cell, sigma, computing time')
for fun in [findSigma_v1, findSigma_v2, findSigma_v3]:
    t0 = time.time(); 
    s = fun(cell, new_dist_df.loc[cell,new_dist_df.columns[new_dist_df.loc[cell]>0]], rho.loc[cell], new_k)
    t1 = time.time()
    print(cell, s, t1-t0)

pandarallel.initialize(nb_workers = 8)
sigma = new_dist_df.parallel_apply(lambda x: findSigma_v1(x.name, new_dist_df.loc[x.name,new_dist_df.columns[new_dist_df.loc[x.name]>0]], rho.loc[x.name], new_k), axis = 1)

# try umap function
#sigmaumap, rhoumap = smooth_knn_dist(np.array(new_dist_df), new_k)
#rho = pd.Series(rhoumap, index = new_dist_df.index)
#sigma = pd.Series(sigmaumap, index = new_dist_df.index)

## 
def directed_weigths(new_dist_cell, rho_cell, sigma_cell):
    knn = new_dist_cell[new_dist_cell>0].index
    cv = pd.Series(0, index = new_dist_cell.index)
    for c in knn:
        cv.loc[c] = np.exp(-(new_dist_cell.loc[c]-rho_cell)/sigma_cell)
    return cv

pandarallel.initialize(nb_workers = 8)
wdf = new_dist_df.parallel_apply(lambda x: directed_weigths(x, rho.loc[x.name], sigma.loc[x.name]), axis = 1)

bdf = wdf + wdf.T - wdf * wdf.T

mdf.obsp['connectivities'] = scipy.sparse.csr_matrix(bdf.values)
mdf.uns['neighbors']['params']['n_neighbors'] = new_k
#cdf = np.exp(-max([0,new_dist_df.T-rho])/sigma)

sc.tl.umap(mdf, n_components = 2, random_state = 921225, min_dist = 0.5, spread = 1)
#sc.pl.umap(mdf, color = 'batch')

#mdf.write(output + '.h5ad')

fig, axs = plt.subplots(ncols = 2, figsize = (2*3*1.6,3))
udfs = pd.DataFrame(mdf.obsm['X_umap_standard'], columns = ['u1','u2'], index = mdf.obs.index)
udfm = pd.DataFrame(mdf.obsm['X_umap'], columns = ['u1','u2'], index = mdf.obs.index)
for ax, udf, title in zip(axs, [udfs, udfm], ['standard','sequential']):
    for t in sorted(set(mdf.obs['batch'])):
        cells = mdf.obs[mdf.obs['batch']==t].index
        ax.scatter(udf.loc[cells,'u1'], udf.loc[cells,'u2'], s = 0.8, label = t, rasterized = True, alpha = 0.8)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(title)
lgn = axs[-1].legend(loc = 2, bbox_to_anchor = (1,1))
for handle in lgn.legendHandles:
    handle.set_sizes([12.0])

fig.savefig(output + '_umap_comparisons.pdf', bbox_inches = 'tight')
mdf.write(output + '_umap_comparisons.h5ad')

sys.exit()




