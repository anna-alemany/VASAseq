#!/usr/bin/env python
# coding: utf-8
import sys, os
import pandas as pd
import numpy as np
from pandas.io.parsers import read_csv
from collections import Counter
import plot_aautils as plaa
import scanpy as sc
import sc_aautils as scaa
from sklearn.metrics import silhouette_score
import scipy.spatial.distance as distance
from pandarallel import pandarallel
import matplotlib.pyplot as plt
plt.ion()

try: 
    timepoint = sys.argv[1]
except:
    sys.exit("Please, provide timepoint: (E65, E75, E85)")

if timepoint == 'E65':
    nCellsClusters = 4
elif timepoint == 'E75':
    nCellsClusters = 9
elif timepoint == 'E85':
    nCellsClusters = 9
lfc_th = 1; pv_th=1e-2
th_score_10x = 0.4; th_score_vasa = 0.2

inputDirUMAPs = '../'+timepoint+'/res_scanpy_merged_2021_final/'
scanpyDirVASA = 'res_scanpy_all_2021_noHistones'
scanpyDirPJ = 'res_scanpy_all_2021_v2'

# # Read UMAPS with all metadata
ux = read_csv(inputDirUMAPs + '/PJ_su_All_UmapMeta_'+timepoint+'.tsv', sep = '\t', index_col = 0)
uv = read_csv(inputDirUMAPs + '/VASA_su_All_UmapMeta_'+timepoint+'.tsv', sep = '\t', index_col = 0)
#plt.hist(uv['x2v_score_3su'], bins = 100, alpha = 0.5, label = 'vasa')
#plt.hist(ux['v2x_score_3su'], bins = 100, alpha = 0.5, label = '10x')

# # Read scanpy objects
vsc = {}
for biotype in ['All','ProteinCoding','lncRNA']:
    for readtype in ['su','s','u']: 
        print(biotype, readtype)
        bio = biotype + '-' + readtype
        vsc[bio] = scaa.read_h5ad_zipped('../'+timepoint+'/'+scanpyDirVASA+'/VASA_'+readtype+'_'+biotype+'_'+timepoint+'.h5ad')
        vsc[bio] = vsc[bio][:, ['-' not in idx for idx in vsc[bio].var.index] ]

xsc = {}
for biotype in ['All','ProteinCoding','lncRNA']:
    for readtype in ['su','s','u']:
        bio = biotype + '-' + readtype
        xsc[bio] = scaa.read_h5ad_zipped('../../pijuansala/'+timepoint+'/'+scanpyDirPJ+'/PJ_'+readtype+'_'+biotype+'_'+timepoint+'.h5ad')
        xsc[bio] = xsc[bio][:, ['-' not in idx for idx in xsc[bio].var.index] ]

# # plot umaps
def plotTransferX2V(axs, ux, uv, th_score = th_score_vasa):
    axs[0].scatter(ux[ 'umap1'], ux['umap2'], c = 'silver', s = 3, label = '_nolabel_', rasterized = True)
    axs[1].scatter(uv[ 'umap1'], uv['umap2'], c = 'silver', s = 3, label = '_nolabel_', rasterized = True)
    for i, cl in enumerate(sorted(set(uv['x2v_assigned_3su']))):
        if cl != -1:
            cx = ux[ux['leiden_clean_manhattan']==cl].index
            cv = uv[(uv['x2v_assigned_3su']==cl)&(uv['x2v_score_3su']>th_score)].index
            axs[0].scatter(ux.loc[cx, 'umap1'], ux.loc[cx,'umap2'], c = plaa.colors()[i], s = 3, label = cl)
            axs[1].scatter(uv.loc[cv, 'umap1'], uv.loc[cv,'umap2'], c = plaa.colors()[i], s = 3, label = cl)
            axs[0].text(ux.loc[cx, 'umap1'].mean(), ux.loc[cx,'umap2'].mean(), int(cl), va = 'center', ha = 'center')
            axs[1].text(uv.loc[cv, 'umap1'].mean(), uv.loc[cv,'umap2'].mean(), int(cl), va = 'center', ha = 'center')
    for ax in axs:
        ax.grid(False); ax.set_xticks([]); ax.set_yticks([])
    axs[0].set_title('10x'); axs[1].set_title('vasa')
    return axs

def plotTransferV2X(axs, ux, uv, th_score = th_score_10x):
    axs[1].scatter(ux[ 'umap1'], ux['umap2'], c = 'silver', s = 3, label = '_nolabel_', rasterized = True)
    axs[0].scatter(uv[ 'umap1'], uv['umap2'], c = 'silver', s = 3, label = '_nolabel_', rasterized = True)
    for i, cl in enumerate(sorted(set(ux['v2x_assigned_3su']))):
        if cl != -1:
            cv = uv[uv['leiden_clean_manhattan']==cl].index
            cx = ux[(ux['v2x_assigned_3su']==cl)&(ux['v2x_score_3su']>th_score)].index
            axs[0].scatter(uv.loc[cv, 'umap1'], uv.loc[cv,'umap2'], c = plaa.colors()[i], s = 3, label = cl)
            axs[1].scatter(ux.loc[cx, 'umap1'], ux.loc[cx,'umap2'], c = plaa.colors()[i], s = 3, label = cl)
            axs[1].text(ux.loc[cx, 'umap1'].mean(), ux.loc[cx,'umap2'].mean(), int(cl), va = 'center', ha = 'center')
            axs[0].text(uv.loc[cv, 'umap1'].mean(), uv.loc[cv,'umap2'].mean(), int(cl), va = 'center', ha = 'center')
    for ax in axs:
        ax.grid(False); ax.set_xticks([]); ax.set_yticks([])
    axs[1].set_title('10x'); axs[0].set_title('vasa')
    return axs

fig, axs = plaa.template_plot(ncols = 2, nrows = 2, figsize = (1.5*2*3*1.6, 1.5*2*3))
plotTransferX2V(axs[0], ux, uv)
plotTransferV2X(axs[1], ux, uv)
plt.savefig(inputDirUMAPs + '/umap_labeltransfers.pdf', bbox_inches = 'tight')
plt.close()

# # Creation of single clusters
fux = ux[ux['v2x_score_3su']>th_score_10x].copy()
fuv = uv[uv['x2v_score_3su']>th_score_vasa].copy()
print(ux.shape, uv.shape, '=> dimensions before score filtering (10x, vasa)')
print(fux.shape, fuv.shape, '=> dimensions after score filtering (10x, vasa)')

def plotTransferMergedV2X(axs, ux, uv, fux, fuv, minCells = 0):  
    axs[1].scatter(ux[ 'umap1'], ux['umap2'], c = 'silver', s = 3, label = '_nolabel_', rasterized = True)
    axs[0].scatter(uv[ 'umap1'], uv['umap2'], c = 'silver', s = 3, label = '_nolabel_', rasterized = True)
    j = 0
    for i, cl in enumerate(sorted(set(list(fux['merged_cluster'])+list(fuv['merged_cluster'])))):
        if cl != -1:
            cv = fuv[fuv['merged_cluster']==cl].index
            cx = fux[fux['merged_cluster']==cl].index
            if len(cv) > minCells and len(cx) > minCells:
                j += 1
                axs[0].scatter(fuv.loc[cv, 'umap1'], fuv.loc[cv,'umap2'], c = plaa.colors()[j], s = 3, label = cl)
                axs[1].scatter(fux.loc[cx, 'umap1'], fux.loc[cx,'umap2'], c = plaa.colors()[j], s = 3, label = cl)
    for ax in axs:
        ax.grid(False); ax.set_xticks([]); ax.set_yticks([])
    axs[1].set_title('10x-'+str(minCells)); axs[0].set_title('vasa-'+str(minCells))
    return axs

def cellCounter(fux, fuv, minCells = 0):  
    jx = 0; jv = 0
    for i, cl in enumerate(sorted(set(list(fux['merged_cluster'])+list(fuv['merged_cluster'])))):
        if cl != -1:
            cv = fuv[fuv['merged_cluster']==cl].index
            cx = fux[fux['merged_cluster']==cl].index
            if len(cv) > minCells and len(cx) > minCells:
                jx += len(cx); jv += len(cv)
    return jx, jv

def cellSelector(fux, fuv, minCells = 0):  
    jx = []; jv = []
    for i, cl in enumerate(sorted(set(list(fux['merged_cluster'])+list(fuv['merged_cluster'])))):
        if cl != -1:
            cv = fuv[fuv['merged_cluster']==cl].index
            cx = fux[fux['merged_cluster']==cl].index
            if len(cv) > minCells and len(cx) > minCells:
                jx += list(cx); jv += list(cv)
    return jx, jv

fux['merged_cluster'] = ['-'.join([str(fux.loc[idx,'leiden_clean_manhattan']),str(int(fux.loc[idx,'v2x_assigned_3su']))]) for idx in fux.index]
fuv['merged_cluster'] = ['-'.join([str(int(fuv.loc[idx,'x2v_assigned_3su'])),str(int(fuv.loc[idx,'leiden_clean_manhattan']))]) for idx in fuv.index]

print(len(set(fux['merged_cluster'])), len(set(fuv['merged_cluster'])),'total merged clusters (10x and vasa, resp.)')

fux.to_csv(inputDirUMAPs + '/PJ_su_All_UmapMeta_MERGEDCLUSTERS_'+timepoint+'.tsv', sep = '\t')
fuv.to_csv(inputDirUMAPs + '/VASA_su_All_UmapMeta_MERGEDCLUSTERS_'+timepoint+'.tsv', sep = '\t')

# SHOULD I TRY TO MAKE CLUSTERS SMALLER HERE, USING THE SCANPY TABLES?

# check the effect of filter by joint-cluster size
fig, axs = plaa.template_plot(ncols = 2, nrows = 2*nCellsClusters, figsize = (1.5*2*3*1.6, 1.5*3*2*nCellsClusters))
for ax, nc in zip(axs, range(1,2*nCellsClusters+1)):
    plotTransferMergedV2X(ax, ux, uv, fux, fuv, minCells=nc)
    ax[1].legend(loc = 2, bbox_to_anchor = (1,1), ncol = 2)
plt.savefig(inputDirUMAPs + '/umap_labeltransfer_ncellsFilters.pdf', bbox_inches = 'tight')
plt.close('all')

cdf = pd.DataFrame({cl: {'10x': Counter(fux['merged_cluster'])[cl], 'vasa': Counter(fuv['merged_cluster'])[cl]} for cl in set(list(set(fux['merged_cluster']))+list(set(fuv['merged_cluster'])))}).T
cdf = cdf.loc[cdf.sum(axis=1).sort_values(ascending=False).index]
cdf.to_csv(inputDirUMAPs + '/ncells_merged-cluster.tsv', sep = '\t')
ncdf = pd.DataFrame({i: cdf.loc[cdf.index[(cdf>i).sum(axis=1)==2]].sum() for i in range(25)}).T
ncdf['clusters'] = [len(cdf.loc[cdf.index[(cdf>i).sum(axis=1)==2]]) for i in ncdf.index]
fig, ax = plaa.template_plot()
im = ax.scatter(ncdf.index, ncdf['10x'], label = '10x', marker = 'o', c = ncdf['clusters'])
ax.scatter(ncdf.index, ncdf['vasa'], label = 'vasa', marker = 'v', c = ncdf['clusters'])
ax.set_xlabel('cell number threshold')
ax.set_ylabel('surviving cells')
plt.colorbar(im, ax = ax, label = 'number of clusters')
ax.legend()
plt.savefig(inputDirUMAPs + '/scatter_survivingCells_bilabel.pdf', bbox_inches = 'tight')
plt.close()

print(cellCounter(fux, fuv, minCells = nCellsClusters))

xcells, vcells = cellSelector(fux, fuv, minCells = nCellsClusters)
ffux = fux.loc[xcells]; ffuv = fuv.loc[vcells]
ncells_surv_df = pd.DataFrame({'10x': Counter(ffux['merged_cluster']), 'vasa': Counter(ffuv['merged_cluster'])})
print(ncells_surv_df)
xcells = []; vcells = []
for cl in ncells_surv_df.index:
    n = ncells_surv_df.loc[cl].min()
    xcells += list(ffux[ffux['merged_cluster'] == cl].sort_values(by='v2x_score_3su', ascending=False).index[:n])
    vcells += list(ffuv[ffuv['merged_cluster'] == cl].sort_values(by='x2v_score_3su', ascending=False).index[:n])
ffux = fux.loc[xcells]; ffuv = fuv.loc[vcells]

# # Silhouette scores and gene expression analysis

# ## read scanpy data 
# ## Filter scanpy objects to select cells for which we have merged_clusters (2 cells at least in each protocol)
# ## also remove multigenes
fvsc = {}; fxsc = {}
for bio in vsc: 
    fvsc[bio] = vsc[bio][[c in ffuv.index for c in vsc[bio].obs.index],:]
    fvsc[bio] = fvsc[bio][:, ['-' not in idx for idx in fvsc[bio].var.index] ]
    fvsc[bio].obs['merged_cluster'] = ffuv.loc[fvsc[bio].obs.index,'merged_cluster']
    print('vasa', bio, vsc[bio].shape, fvsc[bio].shape)
    
for bio in xsc:
    fxsc[bio] = xsc[bio][[c in xcells for c in xsc[bio].obs.index],:]
    fxsc[bio] = fxsc[bio][:, ['-' not in idx for idx in fxsc[bio].var.index] ]
    fxsc[bio].obs['merged_cluster'] = fux.loc[fxsc[bio].obs.index,'merged_cluster']
    print('10x', bio, xsc[bio].shape, fxsc[bio].shape)

# ## Silhouette score
print('average vasa silhouette score:', silhouette_score(fvsc['All-su'].obsm['X_pca'], fvsc['All-su'].obs['merged_cluster'], metric='manhattan'))
print('average 10x silhouette score:', silhouette_score(fxsc['All-su'].obsm['X_pca'], fxsc['All-su'].obs['merged_cluster'], metric='manhattan'))

pca_vasa = pd.DataFrame(fvsc['All-su'].obsm['X_pca'], index = fvsc['All-su'].obs.index)
pca_10x = pd.DataFrame(fxsc['All-su'].obsm['X_pca'], index = fxsc['All-su'].obs.index)
npca = fvsc['All-su'].uns['clean_manhattan']['params']['n_pca']
pca_vasa = pca_vasa.iloc[:, :npca]
npca = fvsc['All-su'].uns['clean_manhattan']['params']['n_pca']
pca_10x = pca_10x.iloc[:, :npca]

def my_silhouette(pcadf, labeldf):
    dist_array = distance.pdist(pcadf, metric = 'cityblock')
    dist_df = pd.DataFrame(distance.squareform(dist_array), columns = pcadf.index, index = pcadf.index)
    pandarallel.initialize(nb_workers = 8)
    adf = dist_df.parallel_apply(lambda x: dist_df.loc[x.name,labeldf[labeldf == labeldf.loc[x.name]].index].sum()/(len(labeldf[labeldf == labeldf.loc[x.name]].index)-1))
    pandarallel.initialize(nb_workers = 8)
    bdf = dist_df.parallel_apply(lambda x: {cl: dist_df.loc[x.name,labeldf[labeldf == cl].index].mean() for cl in set(labeldf) if cl != labeldf.loc[x.name]})
    bdf = pd.DataFrame(bdf.to_dict()).T
    sdf = pd.DataFrame({'b': bdf.min(axis=1), 'a': adf})
    sdf['s'] = sdf.apply(lambda x: (x['b']-x['a'])/np.max([x['b'],x['a']]), axis = 1)
    sdf = sdf.merge(bdf, how = 'inner', left_index = True, right_index = True)
    sdf['merged_cluster'] = labeldf.loc[sdf.index]
    return sdf

sdf_vasa = my_silhouette(pca_vasa, fvsc['All-su'].obs['merged_cluster'])
sdf_10x = my_silhouette(pca_10x , fxsc['All-su'].obs['merged_cluster'])

sdf_10x.to_csv(inputDirUMAPs + '/silhouette_mergedClusters_10x_'+timepoint+'.tsv', sep = '\t')
sdf_vasa.to_csv(inputDirUMAPs + '/silhouette_mergedClusters_VASA_'+timepoint+'.tsv', sep = '\t')

fig, ax = plaa.template_plot()
ax.hist(sdf_vasa['s'], bins = 100, alpha = 0.5, density= True, label = 'vasa')
ax.hist(sdf_10x['s'], bins = 100, alpha = 0.5, density = True, label = '10x')
ax.set_ylabel('density')
ax.set_xlabel('silhouette score')
ax.legend()
plt.savefig(inputDirUMAPs + '/histo_silhouette.pdf', bbox_inches = 'tight')
plt.close()

kk = {cl: {'vasa': sdf_vasa.loc[fvsc['All-su'].obs[fvsc['All-su'].obs['merged_cluster']==cl].index,'s'], 
      '10x': sdf_10x.loc[fxsc['All-su'].obs[fxsc['All-su'].obs['merged_cluster']==cl].index,'s']
      } for cl in sorted(set(fvsc['All-su'].obs['merged_cluster']))}

kk2 = []
for cl in sorted(set(fvsc['All-su'].obs['merged_cluster'])):
    kk2.append(kk[cl]['vasa'].values)
    kk2.append(kk[cl]['10x'].values)

msdf = pd.DataFrame({cl: {
      'vasa-mean': sdf_vasa.loc[fvsc['All-su'].obs[fvsc['All-su'].obs['merged_cluster']==cl].index,'s'].mean(), 
      '10x-mean': sdf_10x.loc[fxsc['All-su'].obs[fxsc['All-su'].obs['merged_cluster']==cl].index,'s'].mean(),
      'vasa-std': sdf_vasa.loc[fvsc['All-su'].obs[fvsc['All-su'].obs['merged_cluster']==cl].index,'s'].std(), 
      '10x-std': sdf_10x.loc[fxsc['All-su'].obs[fxsc['All-su'].obs['merged_cluster']==cl].index,'s'].std(),
      } for cl in sorted(set(fvsc['All-su'].obs['merged_cluster']))}).T

fig, ax = plaa.template_plot(figsize = (3.6*3*1.6,3))
for i, cl in enumerate(sorted(msdf.index)):
    if i == 0:
        ax.scatter([i+0.1], [msdf.loc[cl,'vasa-mean']], c = plaa.colors()[i], marker = 'o', label = 'vasa')
        ax.scatter([i-0.1], [msdf.loc[cl,'10x-mean']],  c = plaa.colors()[i], marker = 'x', label = '10x')
    else:
        ax.scatter([i+0.1], [msdf.loc[cl,'vasa-mean']], c = plaa.colors()[i], marker = 'o', label = '_nolabel_')
        ax.scatter([i-0.1], [msdf.loc[cl,'10x-mean']],  c = plaa.colors()[i], marker = 'x', label = '_nolabel_')
        
    ax.errorbar([i+0.1], [msdf.loc[cl,'vasa-mean']], yerr = [msdf.loc[cl,'vasa-std']], c = plaa.colors()[i])
    ax.errorbar([i-0.1], [msdf.loc[cl,'10x-mean']], yerr = [msdf.loc[cl,'10x-std']], c = plaa.colors()[i])
ax.set_xticks(range(len(msdf.index))); ax.set_xticklabels(msdf.index, rotation = 90)
ax.set_xlabel('clusters')
ax.set_ylabel('mean silhouette score')
ax.legend()
plt.savefig(inputDirUMAPs + '/scatter_silhouette.pdf', bbox_inches = 'tight')
plt.close()

### FIND A WAY TO REDUCE CLUSTER NUMBER IF NEEDED
### to do ###

# ## differential gene expression analysis
for bio in fvsc:
    print([(bio, fvsc[bio].raw.shape, fvsc[bio].shape) for bio in fvsc])

msc = {}
for bio in fxsc:
    rvsc = scaa.getRawData(fvsc[bio])
    rxsc = scaa.getRawData(fxsc[bio])
    msc[bio] = scaa.concatSCobj(rvsc, rxsc, ['vasa','10x'])
    print(bio, msc[bio].shape, rvsc.shape, rxsc.shape, rxsc.shape[0]+rvsc.shape[0]-msc[bio].shape[0])
    print('merged_clusters', len(set(msc[bio].obs['merged_cluster'])))

# ## Differential gene expression between clusters
dex_vasa = {}; dex_10x = {}
for bio in fvsc:
    print(bio,'vasa')
    dex_vasa[bio] = scaa.difGeneExpr(fvsc[bio], clusters = ['merged_cluster'])['merged_cluster']
for bio in fxsc:
    print(bio,'10x')
    dex_10x[bio] = scaa.difGeneExpr(fxsc[bio], clusters = ['merged_cluster'])['merged_cluster']

for bio in dex_vasa:
    os.system('mkdir ' + inputDirUMAPs + '/dexTables/')
    with pd.ExcelWriter(inputDirUMAPs + '/dexTables/dex_vasa_'+bio+'_'+timepoint+'.xlsx') as writer:  
        for cl in dex_vasa[bio]:
            dex_vasa[bio][cl].to_excel(writer, sheet_name=cl)

for bio in dex_10x:
    with pd.ExcelWriter(inputDirUMAPs + '/dexTables/dex_10x_'+bio+'_'+timepoint+'.xlsx') as writer:  
        for cl in dex_10x[bio]:
            dex_10x[bio][cl].to_excel(writer, sheet_name=cl)

# ## Differential gene expression analysis between protocols
dex_x2v = {}
dex_v2x = {}
for bio in dex_10x:
    dex_x2v[bio] = {cl: '' for cl in set(msc[bio].obs['merged_cluster'])}
    dex_v2x[bio] = {cl: '' for cl in set(msc[bio].obs['merged_cluster'])}
    keys = ['names', 'logfoldchanges', 'pvals', 'pvals_adj', 'scores']
    for cl in set(msc[bio].obs['merged_cluster']):
        tmp_df = msc[bio][msc[bio].obs['merged_cluster']==cl,:]
        sc.tl.rank_genes_groups(tmp_df, 'batch', groups=['vasa'], reference='10x', method='t-test', n_genes = 500)
        dex_v2x[bio][cl] = pd.DataFrame({k: [tmp_df.uns['rank_genes_groups'][k][i][0] for i in range(200)] for k in keys}).set_index('names')
        sc.tl.rank_genes_groups(tmp_df, 'batch', groups=['10x'], reference='vasa', method='t-test', n_genes = 500)
        dex_x2v[bio][cl] = pd.DataFrame({k: [tmp_df.uns['rank_genes_groups'][k][i][0] for i in range(200)] for k in keys}).set_index('names')
        with pd.ExcelWriter(inputDirUMAPs + '/dexTables/dex_'+bio+'_v2x_cl'+'_'.join([cl,timepoint])+'.xlsx') as writer:  
            dex_v2x[bio][cl].to_excel(writer, sheet_name='ref-10x')
            dex_x2v[bio][cl].to_excel(writer, sheet_name='ref-vasa')
        
fig, axs = plaa.template_plot(ncols = 2, figsize=(2*2*3*1.6,2*3))
for ax, df in zip(axs, [ux, uv]):
    ax.scatter(df['umap1'], df['umap2'], s = 1, color = 'silver')
    ax.grid(False); ax.set_xticks([]); ax.set_yticks([])
    
for ax, df in zip(axs, [fux, fuv]):
    for i, cl in enumerate(sorted(set(fux.loc[xcells]['merged_cluster']))):
        cells = df[df['merged_cluster']==cl].index
        ax.scatter(df.loc[cells,'umap1'], df.loc[cells,'umap2'], s = 5, color = plaa.colors()[i], label = cl)        
        ax.text(df.loc[cells,'umap1'].mean(),df.loc[cells,'umap2'].mean(), cl, fontsize = 8, va = 'center', ha = 'center') 
lgnd = axs[0].legend(loc = 2, bbox_to_anchor = (0,0), ncol = 10)
for handle in lgnd.legendHandles:
    handle.set_sizes([25.0])
axs[0].set_title('10x'); axs[1].set_title('VASA')
plt.savefig(inputDirUMAPs + '/umap_filtered-merged-clusters.pdf', bbox_inches = 'tight')
plt.close('all')

# ## Heatmap of differential gene expression analysis... what do we show here? 

lfc_dex_df = {}
for bio in dex_10x:
    lfc_dex_df[bio] = pd.DataFrame(columns = ['10x','vasa'], index = sorted(dex_10x[bio]))
    for cl in sorted(dex_10x[bio]):
        fdexvasa = dex_vasa[bio][cl][(dex_vasa[bio][cl]['logfoldchanges']>lfc_th)&(dex_vasa[bio][cl]['pvals_adj']<pv_th)]
        fdex10x = dex_10x[bio][cl][(dex_10x[bio][cl]['logfoldchanges']>lfc_th)&(dex_10x[bio][cl]['pvals_adj']<pv_th)]
        lfc_dex_df[bio].loc[cl] = [len(fdex10x), len(fdexvasa)]

fig, maxs = plaa.template_plot(ncols = 3, nrows = 3, figsize=(3*3*1.6,3*3))
axs = maxs.reshape(9)
for bio, ax in zip(lfc_dex_df, axs):
    for i, idx in enumerate(lfc_dex_df[bio].index):
        ax.text(lfc_dex_df[bio].loc[idx,'vasa'], lfc_dex_df[bio].loc[idx,'10x'], idx, va = 'center', ha = 'center')
        ax.scatter([lfc_dex_df[bio].loc[idx,'vasa']], [lfc_dex_df[bio].loc[idx,'10x']], c = plaa.colors()[i], s = 75)
    ax.plot(np.linspace(0,lfc_dex_df[bio].max().max(),200), np.linspace(0,lfc_dex_df[bio].max().max(),200), c = 'silver')
    ax.set_title(bio)
maxs[-1][1].set_xlabel('total diff expr genes (vasa)')
maxs[1][0].set_ylabel('total diff expr genes (10x)')
plt.savefig(inputDirUMAPs + '/scatter_detectedMarkers.pdf', bbox_inches = 'tight')
plt.close('all')

lfc_df = pd.DataFrame({bio: np.log2((lfc_dex_df[bio]['vasa']+1e-3)/(lfc_dex_df[bio]['10x']+1e-3)) for bio in lfc_dex_df})
fig, ax = plaa.template_plot(figsize = (3*1.6, 3*len(lfc_df)/20))
im = ax.imshow(lfc_df, aspect = 'auto', cmap = 'RdBu_r', vmin = -5, vmax = 5)
cbar = plt.colorbar(im, ax = ax, label = 'log2(vasa markers / 10x markers)', shrink = 0.8)
ax.set_xticks(range(len(lfc_df.columns))); ax.set_xticklabels(lfc_df.columns, rotation = 90)
ax.grid(False)
ax.set_ylabel('merged cluster'); ax.set_yticks(range(len(lfc_df.index))); ax.set_yticklabels(lfc_df.index)
cxra = cbar.get_ticks()
cbar.set_ticks(cxra)
cbar.set_ticklabels(['<'+str(cxra[0])] + [x for x in cxra[1:-1]] + ['>'+str(cxra[-1])])
plt.savefig(inputDirUMAPs + '/heatmap_detectedMarkers.pdf', bbox_inches = 'tight')



