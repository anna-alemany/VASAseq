#!usr/bin/env python3
import sys, os
import numpy as np
import pandas as pd

def filterParams_allCov(timepoint):
    """Filtering parameters for VASA seq data, timepoints E65, E75, E85 and E95"""
    if timepoint == 'E65':
        sample_out = []
        reads_su_th = {'All': [10**4,10**6]}
        fracs_su_th = {'ProteinCoding': [.85,.95], 'lncRNA': [.01,.03], 'smallRNA': [.05,0.15]}
        n_pca = 20
        resolution_s = {'All': 0.7, 'ProteinCoding': 0.65, 'lncRNA': 0.75, 'smallRNA': 0.5, 'TF': 0.7, 'Cofactor': 0.7}
        resolution_su = {'All': 0.8, 'ProteinCoding': 0.65, 'lncRNA': 0.75, 'smallRNA': 0.5, 'TF': 0.7, 'Cofactor': 0.7}
        min_disp = 1
    elif timepoint == 'E75':
        sample_out = []
        reads_su_th = {'All': [10**3.5,10**6]}
        fracs_su_th = {'ProteinCoding': [.85,.95], 'lncRNA': [.01,.03], 'smallRNA': [.05,0.15]}
        n_pca = 40
        resolution_s = {'All': 1, 'ProteinCoding': 1, 'lncRNA': 1, 'smallRNA': 0.75, 'TF': 1, 'Cofactor': 0.75}
        resolution_su =  {'All': 1, 'ProteinCoding': 1, 'lncRNA': 1, 'smallRNA': 0.75, 'TF': 1, 'Cofactor': 0.75}
        min_disp = 0.5
    elif timepoint == 'E85':
        sample_out = ['E8.5-8_i21']
        reads_su_th = {'All': [10**3.5,10**6]}
        fracs_su_th = {'ProteinCoding': [.85,.95], 'lncRNA': [.01,.03], 'smallRNA': [.05,0.15]}
        n_pca = 40
        resolution_s = {'All': 1, 'ProteinCoding': 1, 'lncRNA': 1, 'smallRNA': 0.75, 'TF': 1, 'Cofactor': 0.75}
        resolution_su =  {'All': 1, 'ProteinCoding': 1, 'lncRNA': 1, 'smallRNA': 0.75, 'TF': 1, 'Cofactor': 0.75}
        min_disp = 0.5
    elif timepoint == 'E95':
        sample_out = []
        reads_su_th = {'All': [10**3,10**6]}
        fracs_su_th = {'ProteinCoding': [.85,.95], 'lncRNA': [.012,.03], 'smallRNA': [.05,0.15]}
        n_pca = 40
        resolution_s = {'All': 1, 'ProteinCoding': 1, 'lncRNA': 1, 'smallRNA': 0.75, 'TF': 1, 'Cofactor': 0.75}
        resolution_su = {'All': 1, 'ProteinCoding': 1, 'lncRNA': 1, 'smallRNA': 0.75, 'TF': 1, 'Cofactor': 0.75}
        min_disp = 0.5
    return sample_out, reads_su_th, fracs_su_th, n_pca, resolution_s, resolution_su, min_disp

def filterParams_highCov(timepoint):
    """Filtering parameters for VASA seq data, timepoints E65, E75, E85 and E95"""
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
    return sample_out, reads_su_th, fracs_su_th, n_pca, resolution_s, resolution_su, min_disp
