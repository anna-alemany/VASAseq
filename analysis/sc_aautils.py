import sys, os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
from scipy.cluster import hierarchy
from scipy.cluster.hierarchy import dendrogram, linkage
from pandas.io.parsers import read_csv
from matplotlib.patches import Rectangle, Circle
import itertools as it
import scipy.stats as stats
from collections import Counter

#### hierarchical clustering
def hierarchicalClustering(df, cth = 100, plot = False, method = 'ward', metric = 'euclidean', nolabels = 'True', leaf_colors = []):
    """performs hierarchical clustering using linkage and dendogram functions from scipy.cluster.hierarchy package"""
    if len(leaf_colors) > 0:
        hierarchy.set_link_color_palette(leaf_colors)
    Z = linkage(df, method=method, metric = metric)
    dg = dendrogram(Z, no_labels=nolabels, color_threshold=cth, no_plot = np.invert(plot))
    plt.show()
    return Z, dg

def getClusterByColor(dg, labels):
    """given a dendogram and labels, it groups labels by colors in the dendogram (i.e. clusters)"""
    kk = []
    ii = 0
    cluster = 0
    color = dg['color_list'][0]
    clusters = {cluster: []}
    for i in range(len(dg['icoord'])):
        v = dg['icoord'][i]
        for j in [0,2]:
            vj = int(round((v[j]-5.)/10))
            if (v[j]-5.)/10 == vj and vj not in kk:
                kk.append(vj)
                if dg['color_list'][i] == color:
                    clusters[cluster].append(labels[dg['leaves'][vj]])
                else:
                    color = dg['color_list'][i]
                    cluster += 1
                    clusters[cluster] = [labels[dg['leaves'][vj]]]
    return clusters

#### noise and coefficient of variation
def distance_point_line(xp, yp, a, b):
    """point coordinates are (xp,yp); straight line goes like y=a+b*x"""
    x = (xp+b*yp-a*b)/(1+b*b); y = a+b*x
    d2 = (x-xp)**2 + (y-yp)**2
    return np.sqrt(d2)

def coef_variation_df(ndf):
    """returns a dataframe with mu, SD, variance, coefficient of variance, and distance from Poissan behavior for each index"""
    cvdf = pd.DataFrame({'mu': ndf.mean(axis=1), 'SD': ndf.std(axis=1), 'var': ndf.var(axis=1)})
    cvdf['CV'] = cvdf.apply(lambda x: x['SD']/x['mu'], axis = 1)
    cvdf['distance'] = cvdf.apply(lambda x: distance_point_line(np.log10(x['mu']), np.log10(x['CV']), 0., -0.5), axis = 1)
    return cvdf

#### data treatment 
def zscore(df):
    df = df.loc[df.index[df.sum(axis=1)>0]]
    zdf = df.T
    zdf = (zdf-zdf.mean())/zdf.std()
    zdf = zdf.T
    return zdf

def scalezscore(zdf):
    df  = ((zdf.T>=0)*zdf.T/zdf.max(axis=1) + (zdf.T<0)*zdf.T/abs(zdf.min(axis=1))).T
    return df

def mergeDF(df1, df2, label1, label2):
    df1.columns = [c + '-' + label1 for c in df1.columns]
    df2.columns = [c + '-' + label2 for c in df2.columns]
    mdf = df1.merge(df2, how = 'outer', left_index = True, right_index = True)
    return mdf

#### kmedoids
def kmedoids(dist, numclusters):
    func_min = lambda x, medoids: medoids[dist.loc[medoids,x] == dist.loc[medoids,x].min()]
    findMedoid = lambda idxs, dists: pd.Series({idx: dist.loc[idx, idxs].sum() for idx in idxs}).sort_values().index[0]

    vj = pd.Series()
    for col in dist.columns:
        c = dist.loc[col].sum()
        vj.loc[col] = dist[col].sum()/c

    medoids = np.array(vj.sort_values(ascending=True).index[:numclusters])
    clusters = pd.Series({col: func_min(col, medoids)[0] for col in dist.columns})
    D0 = sum([dist.loc[col, clusters[col]] for col in clusters.index])

    while 2>1:
        medoids_new = []
        for m in medoids:
            idxs = clusters[clusters == m].index
            medoids_new.append(findMedoid(idxs, dist))
        medoids_new = np.array(medoids_new)
        clusters_new = pd.Series({col: func_min(col, medoids_new)[0] for col in dist.columns})
        
        D1 = sum([dist.loc[col, clusters_new[col]] for col in clusters_new.index])
        
        medoids = medoids_new
        clusters = clusters_new
        print(D0,D1)
        if D1 < D0:
            D0 = D1
        else:
            break
    return medoids, clusters

# count combinations
def count_cooccurrences(df, col1, col2):
    cnt_df = pd.DataFrame({v1: Counter(df[(df[col1]==v1)][col2]) for v1 in set(df[col1])}).fillna(0).astype(int)
    return cnt_df

def fisher_enrichments_coocurrence(cnt_df):
    pval_less_df = pd.DataFrame(columns = cnt_df.columns, index = cnt_df.index)
    pval_greater_df = pd.DataFrame(columns = cnt_df.columns, index = cnt_df.index)
    for col in cnt_df.columns:
        for idx in cnt_df.index:
            other_cols = [ocol for ocol in cnt_df.columns if ocol != col]
            other_idxs = [oidx for oidx in cnt_df.index if oidx != idx]
            contingency = pd.DataFrame({
                col: {idx: cnt_df.loc[idx, col], 'other_idxs': cnt_df.loc[idx, other_cols].sum()}, 
                'other_cols': {idx: cnt_df.loc[idx, other_cols].sum(), 'other_idxs': cnt_df.loc[other_idxs, other_cols].sum().sum()}
                })
            oddsratio, pvalue_less = stats.fisher_exact(contingency, alternative = 'less')
            oddsratio, pvalue_greater = stats.fisher_exact(contingency, alternative = 'greater')
            pval_less_df.loc[idx, col] = pvalue_less
            pval_greater_df.loc[idx, col] = pvalue_greater
    return pval_less_df, pval_greater_df

# merge scanpy objects
def getRawData(scdata):
    """Returns scanpy object with the raw data"""
    rsc = sc.AnnData(scdata.raw.X)
    rsc.var = scdata.raw.var
    rsc.obs = scdata.obs
    return rsc

def concatSCobj(sc1, sc2, names = ['1','2']):
    """Concatenation of two scanpy objects"""
    msc = sc1.concatenate(sc2, join = 'outer', batch_categories = names, index_unique = None)
    msc.X = np.nan_to_num(msc.X, nan = 0.0, posinf = 0, neginf = 0)
    return msc

# differential gene expression
def difGeneExpr(adata, clusters =  ['leiden_manhattan','louvain_manhattan'], n_genes = 200):
    dex = {}
    for cluster in clusters:
        adata.obs[cluster] = adata.obs[cluster].astype('category')
        print('*'+cluster+'*')
        cnt = Counter(adata.obs[cluster])
        rm_category = [k for k in cnt if cnt[k] == 1]
        if len(rm_category) > 0:
            print(rm_category)
            print(adata.shape)
            adata = adata[[cl not in rm_category for cl in adata.obs[cluster]],:]
            print(adata.shape)
        try:
            sc.tl.rank_genes_groups(adata, cluster, method='t-test', n_genes=n_genes)
            dex[cluster] = {}
            for cl in set(adata.obs[cluster]):
                try:
                    keys = ['names', 'logfoldchanges', 'pvals', 'pvals_adj', 'scores']
                    dex[cluster][cl] = pd.DataFrame({k: [adata.uns['rank_genes_groups'][k][i][cl] for i in range(200)] for k in keys})
                    dex[cluster][cl] = dex[cluster][cl].set_index('names')
                except:
                    continue
        except:
            print('=>it did not work')
            continue
    return dex


def get_distDF(adata, metric = 'manhattan'):
    if len(metric) > 0:
        dist_df = pd.DataFrame.sparse.from_spmatrix(adata.obsp[metric + '_distances'])
    else:
        dist_df = pd.DataFrame.sparse.from_spmatrix(adata.obsp['distances'])
    dist_df = pd.DataFrame(dist_df.values)
    dist_df.columns = adata.obs.index
    dist_df.index = adata.obs.index
    return dist_df

def get_connectDF(adata, metric = 'manhattan'):
    if len(metric) > 0:
        dist_df = pd.DataFrame.sparse.from_spmatrix(adata.obsp[metric + '_connectivities'])
    else:
        dist_df = pd.DataFrame.sparse.from_spmatrix(adata.obsp['connectivities'])
    dist_df = pd.DataFrame(dist_df.values)
    dist_df.columns = adata.obs.index
    dist_df.index = adata.obs.index
    return dist_df

# RNA velocity functions

def get_velocities(adata):
    velo_df =  pd.DataFrame(adata.layers['velocity'])
    velo_df.columns = adata.var.index
    velo_df.index = adata.obs.index
    return velo_df

def spliced(u, alpha, beta, gamma, s0, u0):
    s = (s0-alpha/gamma+(alpha-beta*u0)/(gamma-beta))*(((beta*u-alpha)/(beta*u0-alpha))**(gamma/beta)) + alpha/gamma + (beta*u-alpha)/(gamma-beta)
    return s

def spliced_t(t, alpha, beta, gamma, s0, u0):
    s = 3
    return s

def unspliced_t(t, alpha, beta, gamma, s0, u0):
    u = 3
    return u

def time_u(u, alpha, beta, gamma, s0, u0):
    t = 3
    return t

### read files
def read_h5ad_zipped(adatafile):
    if os.path.isfile(adatafile + '.gz'):
        os.system('gunzip '+adatafile + '.gz')
    adata = sc.read_h5ad(adatafile)
    os.system('gzip '+adatafile)
    return adata

