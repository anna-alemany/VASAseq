#!/usr/bin/env python
# coding: utf-8
import sys, os
import numpy as np
import pandas as pd
from pandas.io.parsers import read_csv
from collections import Counter
import scanpy as sc
import scvelo as scv
import glob

try:
    timepoint = sys.argv[1]
except:
    sys.exit("Please, give timepoint (E65, E75, E85, E95)")

if timepoint == 'E65':
    inputfeatherdir = '/hpc/hub_oudenaarden/aalemany/vasaseq/NovaSeq/E65/all_covs/count_tables_filters/'
elif timepoint == 'E75':
    inputfeatherdir = '/hpc/hub_oudenaarden/aalemany/vasaseq/NovaSeq/E75/allCov/count_tables_with_filters/'
elif timepoint == 'E85':
    inputfeatherdir = '/hpc/hub_oudenaarden/aalemany/vasaseq/NovaSeq/E85/allCov/count_tables_with_filter/'
elif timepoint == 'E95':
    inputfeatherdir = '/hpc/hub_oudenaarden/aalemany/vasaseq/NovaSeq/E95/allCov/count_Tables_with_filters/'

inputScanpy = '../' + timepoint + '/res_scanpy_all_2021_v2'

fileU = glob.glob(inputfeatherdir + '/*shortGeneNames_uniaggGenes_unspliced.TranscriptCounts.feather')
fileS = glob.glob(inputfeatherdir + '/*shortGeneNames_uniaggGenes_spliced.TranscriptCounts.feather')

if len(fileU) != 1 or len(fileS) != 1:
    sys.exit("Input feather files not found" + str(len(fileU)) + ' ' + str(len(fileS)))

udf = pd.read_feather(fileU[0])
sdf = pd.read_feather(fileS[0])

udf = udf.set_index(udf.columns[0])
sdf = sdf.set_index(sdf.columns[0])

if os.path.isfile(inputScanpy + '/VASA_s_All_'+timepoint+'.h5ad.gz'):
    os.system('gunzip '+ inputScanpy + '/VASA_s_All_'+timepoint+'.h5ad.gz')
adata = sc.read_h5ad(inputScanpy + '/VASA_s_All_'+timepoint+'.h5ad')
os.system('gzip '+ inputScanpy + '/VASA_s_All_'+timepoint+'.h5ad')

if os.path.isfile(inputScanpy + '/VASA_su_All_'+timepoint+'.h5ad.gz'):
    os.system('gunzip '+ inputScanpy + '/VASA_su_All_'+timepoint+'.h5ad.gz')
all_sc = sc.read_h5ad(inputScanpy + '/VASA_su_All_'+timepoint+'.h5ad')
os.system('gzip '+ inputScanpy + '/VASA_su_All_'+timepoint+'.h5ad')

cells = Counter(list(adata.obs.index) + list(all_sc.obs.index))
keep = [idx for idx in cells if cells[idx] == 2]

adata = adata[[idx in keep for idx in adata.obs.index],:]
all_sc = all_sc[[idx in keep for idx in all_sc.obs.index],:]

sdf = 1e5*sdf/sdf.sum()
udf = 1e5*udf/udf.sum()

adata.layers['spliced'] = sdf.loc[adata.var.index, adata.obs.index].T
adata.layers['unspliced'] = udf.loc[adata.var.index, adata.obs.index].T

#scv.pl.proportions(adata, groupby = 'leiden_manhattan')
scv.pp.filter_genes(adata, min_shared_counts=15)
scv.pp.moments(adata, n_pcs = 40, n_neighbors = 20)
scv.tl.velocity(adata)
scv.tl.velocity_graph(adata)
adata.obsm['X_umap'] = all_sc.obsm['X_umap_clean_manhattan']
adata.uns['umap'] = all_sc.uns['umap']
adata.obs['leiden_clean_manhattan'] = all_sc.obs['leiden_clean_manhattan']
adata.uns['leiden'] = all_sc.uns['leiden']

adata.write(inputScanpy + '/VASA_RNAvelocity_' + timepoint + '.h5ad')

#scv.pl.velocity_embedding_stream(adata, basis='umap', color = 'leiden', size = 50, alpha = 1, figsize = (2*3*1.6,2*3))
#scv.pl.velocity_embedding(adata, arrow_length=3, arrow_size=2, dpi=120, size = 50, alpha = 1)
#scv.pl.velocity_embedding_grid(adata, size = 50, alpha = 1, arrow_size = 2, arrow_length = 5, figsize = (2*3*1.6,2*3))




