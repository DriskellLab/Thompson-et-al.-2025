#!/usr/bin/env python
# coding: utf-8

# In[1]:


#load packages
import stereo as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colormaps as cmaps
plt.figure(dpi=300)
import warnings
warnings.filterwarnings('ignore')


# In[2]:

data_path = '/scratch/user/sean.thompson/20241203_163132/12.2/PigP10_realign/Pig_P10_2/outs/visualization/visualization/B03701E2.tissue.gef'
st.io.read_gef_info(data_path)


# In[3]:


#load the data to generate Stereopy object from GEF file
#data = st.io.read_gef(file_path=data_path, bin_type='bins')
data = st.io.read_gef(file_path=data_path, bin_type='bins', bin_size=20, gene_name_index=True)
#preview data structure
data


# In[4]:


#extract the gene names from StereoExpData object for later analysis
df = data.to_df()#convert StereoExpData object to dataframe
#extract columns from the dataframe (gene expression matrix)
df.columns
#df.columns[0]#call individual gene name, starts from 0
df.columns.to_frame().T.to_csv('/scratch/user/sean.thompson/20241203_163132/12.3/PigP10/PigP10_bin20_genes.csv', index=False, header=False)


# In[5]:


#QC: looks at counts, n_genes, and % mitochondrial genes, similar to Seurat
data.tl.cal_qc()


# In[6]:


#visualize QC in aggregate
data.plt.violin()
plt.savefig(fname='/scratch/user/sean.thompson/20241203_163132/12.3/PigP10/1_QC_Violin.png', dpi=300)
data.plt.genes_count()
plt.savefig(fname='/scratch/user/sean.thompson/20241203_163132/12.3/PigP10/2_QC_Scatter.png', dpi=300)


# In[7]:


#Visualize QC on the spatial tissue map
data.plt.spatial_scatter(out_dpi=2500)
plt.savefig(fname='/scratch/user/sean.thompson/20241203_163132/12.3/PigP10/3_QC_SpatialScatter.png', dpi=1000)
#data.plt.spatial_scatter(base_image='/scratch/user/sean.thompson/20241203_163132/12.2/PigP10_realign/Pig_P10_2/outs/image/B03701E2_ssDNA_regist.tif', base_im_cmap='RGB', base_im_to_gray=True)


# In[8]:


#do QC filtration
#min_n_genes_by_counts: minimum number of total counts required for a cell to pass fitlering
#min_gene: minimum number of genes expressed required for a cell to pass filtering
#data.tl.filter_cells(
#        min_n_genes_by_counts=50,
#        min_gene=5,
#        pct_counts_mt=10,
#        inplace=True
#        )#in ~stereopy/stereo/core/st_pipeline.py:162, filter_cells is defined with parameters 'min_gene'/'max_gene'/'min_n_genes_by_counts' not 'min_genes'/'max_genes'/'min_counts' like in stereopy vignette
#data

#Visualize QC on the spatial tissue map
#data.plt.spatial_scatter(out_dpi=1000)
#data.plt.spatial_scatter(out_dpi=500)
#plt.savefig(fname='/scratch/user/sean.thompson/20241203_163132/12.3/PigP10/4_QC_Filtered_SpatialScatter.png', dpi=500)


# In[9]:


data.tl.raw_checkpoint()#save the stereopy data object before normalizing
data.tl.raw
#data.tl.reset_raw_data()#restore data object to the saved raw data object


#save the StereoExpData to file
st.io.write_h5ad(data, use_raw=True, use_result=True, key_record=None, output='/scratch/user/sean.thompson/20241203_163132/12.3/h5ad/PigP10_bin20.stereo.h5ad', split_batches=True)
#unfiltered StereoExpData saved to .stereo.h5ad from Kamiak HPC and analysis resumed in pt2 file on local hardware
