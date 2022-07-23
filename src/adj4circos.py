#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ============================================================================================= #
# adj4circos.py                                                                                 #
# Author: Juan Sebastian Diaz Boada                                                             #
# Creation Date: 29/10/2021                                                                     #
# ============================================================================================= #
""" Calculates the adjacency matrix that represents the appereances of combinations between
    different V and J usages of the TCRs. The adjacency matrix is exported as a .csv file to
    generate a circos plot in R.

    Parameters
    ----------
    infile : string.
        Path to the TCR meta-dataset.
    chain : String.
        Name of the chain allele from whom the sequence will be used for the plot. Can take the
        following values: (A,B,G or D).
    patient : int (optional).
        Number of the patient to select from the dataset. If zero, it uses the complete dataset.
        Defaults to zero.
    tissue : string (optional).
        Name of the tissue to select to make the plots. It can be 'MUSL', 'PB' or 'both'. Set as
        'both' by default.
    filter : flag.
        Include this flag to filter samples according to threshold.
    no--filter : flag.
        Include this flag to skip sample filtering. This is the default value of the flag.
    threshold : int.
        Cell threshold: If cells are to be filtered, discards the cells with counts below
        threshold. Set to 1 by default.

    Returns
    -------
    pd.Dataframe
        One-column dataframe with the isolated sequence or usage for each sample of the
        given chain allele.
"""
import argparse
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings("ignore")
#################################################################################################
##################################### PARAMETERS ################################################
#################################################################################################
# Argparser definition
parser = argparse.ArgumentParser(description='Filter parameters')
parser.add_argument('infile', type=str, help="Path of the file containing the TCR dataset")
parser.add_argument('chain', type=str, help="TCR chain to use. Can be 'A', 'B', 'G' or 'D'.")
parser.add_argument('-p','--patient', type=int, default=0, help="Number of the patient to plot. It \
has to be a number from 0-7. If 0, the plotting is done with all patients. Defaults to zero.")
parser.add_argument('-t','--tissue',type=str, default='both',help="Name of the tissue to be \
selected to do the analysis. If 'both', no filtering by tissue is done.")
parser.add_argument('--filter_cells', dest='filter', action='store_true',help="Filter out\ connections with a low number of samples.")
parser.add_argument('--no-filter', dest='filter', action='store_false',help="Do not filter out \
connections with a low number of samples. This is default behaviour.")
parser.set_defaults(filter=False)
parser.add_argument('-c','--threshold', type=int, default=1, help="Cell threshold: If cells are to be \
filtered, discards the cells with counts below threshold. Set to 1 by default.")
args = parser.parse_args()
# Chain
chain = args.chain
chain_options = ['A','B','G','D']
if not chain in chain_options:
    raise NameError("Invalid chain. Has to be either 'A','B','G' or 'D'.")
# Tissue
tissue = args.tissue
filter_tissue = True
tissue_options= ['MUSL','PB','both']
if not tissue in tissue_options:
    raise NameError("Invalid tissue. Has to be either 'MUSL','PB' or 'both'.")
elif tissue == 'both':
    filter_tissue = False
# Patient
pat = args.patient
filter_patient = True
if pat<0 or pat>7:
    raise ValueError("Invalid patient number. Has to be an int between 0-7.")
elif pat==0:
    filter_patient = False
    pat=''
else:
    pat = 'sc' + str(pat)
# Cell threshold
filter_cell_number = args.filter
thresh_cells = args.threshold
# Paths
infile = args.infile
output_path_prefix = '../data/output_data/circos_adj/'
outfile = output_path_prefix + 'circos_adjacency_' + chain + \
         filter_patient*('_' + pat) + filter_tissue*('_' + tissue) + \
         filter_cell_number*('_thresh_' + str(thresh_cells)) + '.csv'
#################################################################################################
##################################### PROCESSING ################################################
#################################################################################################
# LOAD DATA
DF = pd.read_csv(infile,sep='\t',index_col=0)
DF = DF.loc[DF['treatment.status']=='naive'] # Use only naive samples
# Filter by tissue
if filter_tissue:
    DF = DF.loc[DF['tissue']==tissue,:]
# Filter by patient
if filter_patient:
    DF = DF.loc[DF['patient']==pat]
# Select chains and productive samples
    # Allele 1
prod_col_1 = 'TR' + chain + '_1_productive'
DF_1 = DF.loc[DF[prod_col_1]==1,[chain + '_1_V',chain + '_1_J']].dropna(thresh=2)
    # Allele 2
prod_col_2 = 'TR' + chain + '_2_productive'
DF_2 = DF.loc[DF[prod_col_2]==1,[chain + '_2_V',chain + '_2_J']].dropna(thresh=2)
# Concatenate vertically both alleles
col_names = [chain + '_V', chain + '_J'] # V and J column names
DF_1.columns = col_names
DF_2.columns = col_names
DF_VJ = pd.concat([DF_1,DF_2],axis=0)
# Frequency calculation and filtering by cell threshold
DF_VJ['union'] = DF_VJ.iloc[:,0] + DF_VJ.iloc[:,1] # Combination of V and J
DF_VJ['freq'] = DF_VJ.loc[:,'union'].map(DF_VJ.loc[:,'union'].value_counts()) # Frequency VJ combination
DF_VJ = DF_VJ.sort_values(['freq'],ascending=False) # Sort in descending frequency order
if filter_cell_number:
    # Discard cells with fewer appearances than cell threshold
    DF_VJ = DF_VJ.loc[DF_VJ['freq']>thresh_cells]
# Overwrite frequency columns
DF_VJ.loc[:,'freq_v'] = DF_VJ.iloc[:,0].map(DF_VJ.iloc[:,0].value_counts())
DF_VJ.loc[:,'freq_j'] = DF_VJ.iloc[:,1].map(DF_VJ.iloc[:,1].value_counts())
DF_VJ.loc[:,'freq'] = DF_VJ.loc[:,'union'].map(DF_VJ.loc[:,'union'].value_counts())
# ADJACENCY MATRIX
# Unique alleles
V_als = DF_VJ.sort_values('freq_v',ascending=False).iloc[:,0].unique()
J_als = DF_VJ.sort_values('freq_j',ascending=False).iloc[:,1].unique()
# dictionaries of allele to index
v2idx = {v:i for i,v in enumerate(V_als)}
j2idx = {j:i for i,j in enumerate(J_als)}
# Creation of adjacency matrix
adj = np.zeros([len(V_als),len(J_als)],dtype=int)
for i in range(len(DF_VJ)):
    r = v2idx[DF_VJ.iloc[i,0]]
    c = j2idx[DF_VJ.iloc[i,1]]
    adj[r,c] = adj[r,c] + 1
df_adj = pd.DataFrame(adj,columns=j2idx.keys(),index=v2idx.keys(),dtype=int)
#################################################################################################
##################################### EXPORT DATA ###############################################
#################################################################################################
df_adj.to_csv(outfile,sep=',')
print("The dataset of {} chain".format(chain),filter_patient*"for patient {}".format(pat),\
"using {} tissue(s)".format(tissue),\
filter_cell_number*"filtered with a threshold of {}".format(thresh_cells),\
"has {} V alleles and {} J alleles.".format(len(V_als),len(J_als)))
