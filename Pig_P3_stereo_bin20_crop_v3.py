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

data_path = '/scratch/user/sean.thompson/20241203_163132/12.2/PigP3_realign/Pig_P3_3/outs/visualization/visualization/B03703F1.tissue.gef'
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
df.columns.to_frame().T.to_csv('/scratch/user/sean.thompson/20241203_163132/12.3/PigP3/PigP3_bin20_genes.csv', index=False, header=False)


# In[5]:


#QC: looks at counts, n_genes, and % mitochondrial genes, similar to Seurat
data.tl.cal_qc()


# In[6]:


#visualize QC in aggregate
data.plt.violin()
plt.savefig(fname='/scratch/user/sean.thompson/20241203_163132/12.3/PigP3/1_QC_Violin.png', dpi=300)
data.plt.genes_count()
plt.savefig(fname='/scratch/user/sean.thompson/20241203_163132/12.3/PigP3/2_QC_Scatter.png', dpi=300)


# In[7]:


#Visualize QC on the spatial tissue map
data.plt.spatial_scatter(out_dpi=2500)
plt.savefig(fname='/scratch/user/sean.thompson/20241203_163132/12.3/PigP3/3_QC_SpatialScatter.png', dpi=1000)
#data.plt.spatial_scatter(base_image='/scratch/user/sean.thompson/20241203_163132/12.2/PigP3_realign/Pig_P3_3/outs/image/B03703F1_ssDNA_regist.tif', base_im_cmap='RGB', base_im_to_gray=True)


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
#plt.savefig(fname='/scratch/user/sean.thompson/20241203_163132/12.3/PigP3/4_QC_Filtered_SpatialScatter.png', dpi=500)


# In[9]:


data.tl.raw_checkpoint()#save the stereopy data object before normalizing
data.tl.raw
#data.tl.reset_raw_data()#restore data object to the saved raw data object


#save the StereoExpData to file
st.io.write_h5ad(data, use_raw=True, use_result=True, key_record=None, output='/scratch/user/sean.thompson/20241203_163132/12.3/h5ad/PigP3_bin20.stereo.h5ad', split_batches=True)

#save the StereoExpData to anndata (Scanpy compatible)
# remember to set flavor as seurat
#adata = st.io.stereo_to_anndata(data, flavor='scanpy', sample_id='PigP3', output='/scratch/user/sean.thompson/20241203_163132/12.3/h5ad/PigP3_bin20.stereo.forScanpy.h5ad')
#st.io.write_h5ad(adata, use_raw=True, use_result=True, key_record=None, output='/scratch/user/sean.thompson/20241203_163132/12.3/h5ad/PigP3_bin20.stereo.forScanpy.h5ad', split_batches=True)
st.io.stereo_to_anndata(data, flavor='scanpy', sample_id='PigP3', output='/scratch/user/sean.thompson/20241203_163132/12.3/h5ad/PigP3_bin20.stereo.forScanpy.h5ad')

#del adata#clear the object to save memory

#save the StereoExpData to anndata (Seurat compatible)
# remember to set flavor as seurat
#adata = st.io.stereo_to_anndata(data, flavor='seurat', sample_id='PigP3', output='/scratch/user/sean.thompson/20241203_163132/12.3/h5ad/PigP3_bin20.stereo.forSeurat.h5ad')
#st.io.write_h5ad(adata, use_raw=True, use_result=True, key_record=None, output='/scratch/user/sean.thompson/20241203_163132/12.3/h5ad/PigP3_bin20.stereo.forSeurat.h5ad', split_batches=True)
st.io.stereo_to_anndata(data, flavor='seurat', sample_id='PigP3', output='/scratch/user/sean.thompson/20241203_163132/12.3/h5ad/PigP3_bin20.stereo.forSeurat.h5ad')

#del adata#clear the object to save memory




# In[10]:


# In[8]:


#do QC filtration
#min_n_genes_by_counts: minimum number of total counts required for a cell to pass fitlering
#min_gene: minimum number of genes expressed required for a cell to pass filtering
data.tl.filter_cells(
        min_n_genes_by_counts=50,
        min_gene=5,
        pct_counts_mt=10,
        inplace=True
        )#in ~stereopy/stereo/core/st_pipeline.py:162, filter_cells is defined with parameters 'min_gene'/'max_gene'/'min_n_genes_by_counts' not 'min_genes'/'max_genes'/'min_counts' like in stereopy vignette
data

#Visualize QC on the spatial tissue map
#data.plt.spatial_scatter(out_dpi=1000)
data.plt.spatial_scatter(out_dpi=500)
plt.savefig(fname='/scratch/user/sean.thompson/20241203_163132/12.3/PigP3/4_QC_Filtered_SpatialScatter.png', dpi=500)


# inplace is set to True by default
#normalize data
#data.tl.normalize_total()
#data.tl.log1p()


# In[11]:


#use scTransform to normalize data and identify variable genes
data.tl.sctransform(res_key='sctransform', inplace=True, filter_hvgs=False, n_cells=5000, n_genes=2000)
#https://stereopy.readthedocs.io/en/latest/Tutorials/scTransform.html


# In[12]:


# run PCA and identify number of PCs to use for clustering
data.tl.pca(use_highly_genes=False, hvg_res_key='highly_variable_genes', n_pcs=30, res_key='pca', svd_solver='arpack')
data.plt.elbow(pca_res_key='pca')
plt.savefig(fname='/scratch/user/sean.thompson/20241203_163132/12.3/PigP3/5_PC_ElbowPlot.png', dpi=300)


# In[13]:


data.tl.neighbors(pca_res_key='pca', n_pcs=20, res_key='neighbors', n_jobs=6)
data.tl.umap(pca_res_key='pca', neighbors_res_key='neighbors', res_key='umap', init_pos='spectral', spread=2.0)


# In[14]:


#visualize expression of marker genes across the UMAP
data.plt.umap(gene_names=['KRT15', 'KRT10'], res_key='umap', out_dpi=500)
plt.savefig(fname='/scratch/user/sean.thompson/20241203_163132/12.3/PigP3/6_KRT15_KRT10_UMAP.png', dpi=500)

data.plt.umap(gene_names=['COL1A1','COL3A1'], res_key='umap', out_dpi=500)
plt.savefig(fname='/scratch/user/sean.thompson/20241203_163132/12.3/PigP3/7_COL1A1_COL3A1_UMAP.png', dpi=500)

data.plt.umap(gene_names=['VEGFA','PDGFC'], res_key='umap', out_dpi=500)
plt.savefig(fname='/scratch/user/sean.thompson/20241203_163132/12.3/PigP3/8_VEGFA_PDGFC_UMAP.png', dpi=500)


# In[15]:


#Leiden clustering
data.tl.leiden(neighbors_res_key='neighbors',res_key='leiden', resolution=1)#resolution=1 is default


# In[16]:


#view Leiden clustering onto the mask of the tissue section
data.plt.cluster_scatter(res_key='leiden')
plt.savefig(fname='/scratch/user/sean.thompson/20241203_163132/12.3/PigP3/9_Spatial_Clustering.png', dpi=500)
data.plt.cluster_scatter(res_key='leiden', base_image='/scratch/user/sean.thompson/20241203_163132/12.2/PigP3_realign/Pig_P3_3/outs/image/B03703F1_ssDNA_regist.tif')
plt.savefig(fname='/scratch/user/sean.thompson/20241203_163132/12.3/PigP3/10_Spatial_Clustering_wMask.png', dpi=500)


# In[17]:


#view subset of Leiden clustering onto the mask of the tissue section
#data.plt.cluster_scatter(res_key='leiden', groups=['1', '2', '3', '4', '5', '8', '9', '10',
#                                                  '11', '12', '13', '14', '15', '16', '17', 
#                                                   '18', '19', '20', '21', '22'])
#plt.savefig(fname='/scratch/user/sean.thompson/20241203_163132/12.3/PigP3/11_filtered_Spatial_Clustering.png', dpi=500)
#data.plt.cluster_scatter(res_key='leiden', groups=['1', '2', '3', '4', '5', '8', '9', '10',
#                                                  '11', '12', '13', '14', '15', '16', '17', 
#                                                   '18', '19', '20', '21', '22'], base_image='/scratch/user/sean.thompson/20241203_163132/12.2/PigP3_realign/Pig_P3_3/outs/image/B03703F1_ssDNA_regist.tif')
#plt.savefig(fname='/scratch/user/sean.thompson/20241203_163132/12.3/PigP3/12_filtered_Spatial_Clustering_wMask.png', dpi=500)


# In[18]


data.plt.umap(res_key='umap', cluster_key='leiden')
plt.savefig(fname='/scratch/user/sean.thompson/20241203_163132/12.3/PigP3/13_Leiden_UMAP.png', dpi=500)


# In[19]


#view spatial gene expression across the tissue mask
data.plt.spatial_scatter_by_gene(gene_name='KRT14', palette='linear_grey_10_95_c0',out_dpi=500, height=1000, width=1000)
plt.savefig(fname='/scratch/user/sean.thompson/20241203_163132/12.3/PigP3/KRT14_Spatial_Expression.png', dpi=500)
data.plt.spatial_scatter_by_gene(gene_name='KRT15', palette='linear_grey_10_95_c0',out_dpi=500, height=1000, width=1000)
plt.savefig(fname='/scratch/user/sean.thompson/20241203_163132/12.3/PigP3/KRT15_Spatial_Expression.png', dpi=500)
data.plt.spatial_scatter_by_gene(gene_name='KRT10', palette='linear_grey_10_95_c0',out_dpi=500, height=1000, width=1000)
plt.savefig(fname='/scratch/user/sean.thompson/20241203_163132/12.3/PigP3/KRT10_Spatial_Expression.png', dpi=500)
data.plt.spatial_scatter_by_gene(gene_name='COL1A1', palette='linear_grey_10_95_c0',out_dpi=500, height=1000, width=1000)
plt.savefig(fname='/scratch/user/sean.thompson/20241203_163132/12.3/PigP3/COL1A1_Spatial_Expression.png', dpi=500)
data.plt.spatial_scatter_by_gene(gene_name='COL3A1', palette='linear_grey_10_95_c0',out_dpi=500, height=1000, width=1000)
plt.savefig(fname='/scratch/user/sean.thompson/20241203_163132/12.3/PigP3/COL3A1_Spatial_Expression.png', dpi=500)
data.plt.spatial_scatter_by_gene(gene_name='VEGFA', palette='linear_grey_10_95_c0',out_dpi=500, height=1000, width=1000)
plt.savefig(fname='/scratch/user/sean.thompson/20241203_163132/12.3/PigP3/VEGFA_Spatial_Expression.png', dpi=500)
data.plt.spatial_scatter_by_gene(gene_name='PDGFC', palette='linear_grey_10_95_c0',out_dpi=500, height=1000, width=1000)
plt.savefig(fname='/scratch/user/sean.thompson/20241203_163132/12.3/PigP3/PDGFC_Spatial_Expression.png', dpi=500)
data.plt.spatial_scatter_by_gene(gene_name='PECAM1', palette='linear_grey_10_95_c0',out_dpi=500, height=1000, width=1000)
plt.savefig(fname='/scratch/user/sean.thompson/20241203_163132/12.3/PigP3/PECAM1_Spatial_Expression.png', dpi=500)
data.plt.spatial_scatter_by_gene(gene_name='CDH5', palette='linear_grey_10_95_c0',out_dpi=500, height=1000, width=1000)
plt.savefig(fname='/scratch/user/sean.thompson/20241203_163132/12.3/PigP3/CDH5_Spatial_Expression.png', dpi=500)
data.plt.spatial_scatter_by_gene(gene_name='ACTA2', palette='linear_grey_10_95_c0',out_dpi=500, height=1000, width=1000)
plt.savefig(fname='/scratch/user/sean.thompson/20241203_163132/12.3/PigP3/ACTA2_Spatial_Expression.png', dpi=500)
data.plt.spatial_scatter_by_gene(gene_name='RGS5', palette='linear_grey_10_95_c0',out_dpi=500, height=1000, width=1000)
plt.savefig(fname='/scratch/user/sean.thompson/20241203_163132/12.3/PigP3/RGS5_Spatial_Expression.png', dpi=500)


# In[20]


#Leiden clustering optimization
#0.8 (1.0 default)
data.tl.leiden(neighbors_res_key='neighbors',res_key='leiden', resolution=0.8)#resolution=1 is default
#view Leiden clustering onto the mask of the tissue section
data.plt.cluster_scatter(res_key='leiden')
plt.savefig(fname='/scratch/user/sean.thompson/20241203_163132/12.3/PigP3/14_Spatial_Clustering_res0.8.png', dpi=500)
data.plt.cluster_scatter(res_key='leiden', base_image='/scratch/user/sean.thompson/20241203_163132/12.2/PigP3_realign/Pig_P3_3/outs/image/B03703F1_ssDNA_regist.tif')
plt.savefig(fname='/scratch/user/sean.thompson/20241203_163132/12.3/PigP3/15_Spatial_Clustering_wMask_res0.8.png', dpi=500)
data.plt.umap(res_key='umap', cluster_key='leiden')
plt.savefig(fname='/scratch/user/sean.thompson/20241203_163132/12.3/PigP3/16_Leiden_UMAP_res0.8.png', dpi=500)

#0.6 (1.0 default)
data.tl.leiden(neighbors_res_key='neighbors',res_key='leiden', resolution=0.6)#resolution=1 is default
#view Leiden clustering onto the mask of the tissue section
data.plt.cluster_scatter(res_key='leiden')
plt.savefig(fname='/scratch/user/sean.thompson/20241203_163132/12.3/PigP3/17_Spatial_Clustering_res0.6.png', dpi=500)
data.plt.cluster_scatter(res_key='leiden', base_image='/scratch/user/sean.thompson/20241203_163132/12.2/PigP3_realign/Pig_P3_3/outs/image/B03703F1_ssDNA_regist.tif')
plt.savefig(fname='/scratch/user/sean.thompson/20241203_163132/12.3/PigP3/18_Spatial_Clustering_wMask_res0.6.png', dpi=500)
data.plt.umap(res_key='umap', cluster_key='leiden')
plt.savefig(fname='/scratch/user/sean.thompson/20241203_163132/12.3/PigP3/19_Leiden_UMAP_res0.8.png', dpi=500)

#0.4 (1.0 default)
data.tl.leiden(neighbors_res_key='neighbors',res_key='leiden', resolution=0.4)#resolution=1 is default
#view Leiden clustering onto the mask of the tissue section
data.plt.cluster_scatter(res_key='leiden')
plt.savefig(fname='/scratch/user/sean.thompson/20241203_163132/12.3/PigP3/20_Spatial_Clustering_res0.4.png', dpi=500)
data.plt.cluster_scatter(res_key='leiden', base_image='/scratch/user/sean.thompson/20241203_163132/12.2/PigP3_realign/Pig_P3_3/outs/image/B03703F1_ssDNA_regist.tif')
plt.savefig(fname='/scratch/user/sean.thompson/20241203_163132/12.3/PigP3/21_Spatial_Clustering_wMask_res0.4.png', dpi=500)
data.plt.umap(res_key='umap', cluster_key='leiden')
plt.savefig(fname='/scratch/user/sean.thompson/20241203_163132/12.3/PigP3/22_Leiden_UMAP_res0.4.png', dpi=500)

#1.2 (1.0 default)
data.tl.leiden(neighbors_res_key='neighbors',res_key='leiden', resolution=1.2)#resolution=1 is default
#view Leiden clustering onto the mask of the tissue section
data.plt.cluster_scatter(res_key='leiden')
plt.savefig(fname='/scratch/user/sean.thompson/20241203_163132/12.3/PigP3/23_Spatial_Clustering_res1.2.png', dpi=500)
data.plt.cluster_scatter(res_key='leiden', base_image='/scratch/user/sean.thompson/20241203_163132/12.2/PigP3_realign/Pig_P3_3/outs/image/B03703F1_ssDNA_regist.tif')
plt.savefig(fname='/scratch/user/sean.thompson/20241203_163132/12.3/PigP3/24_Spatial_Clustering_wMask_res1.2.png', dpi=500)
data.plt.umap(res_key='umap', cluster_key='leiden')
plt.savefig(fname='/scratch/user/sean.thompson/20241203_163132/12.3/PigP3/25_Leiden_UMAP_res1.2.png', dpi=500)


# In[21]


#save the StereoExpData to file
st.io.write_h5ad(data, use_raw=True, use_result=True, key_record=None, output='/scratch/user/sean.thompson/20241203_163132/12.3/h5ad/PigP3_bin20_normalized.stereo.h5ad', split_batches=True)

#save the StereoExpData to anndata (Seurat compatible)
# remember to set flavor as seurat
#adata = st.io.stereo_to_anndata(data, flavor='seurat', sample_id='PigP3', output='/scratch/user/sean.thompson/20241203_163132/12.3/h5ad/PigP3_bin20_normalized.stereo.forSeurat.h5ad')
#st.io.write_h5ad(adata, use_raw=True, use_result=True, key_record=None, output='/scratch/user/sean.thompson/20241203_163132/12.3/h5ad/PigP3_bin20_normalized.stereo.forSeurat.h5ad', split_batches=True)
st.io.stereo_to_anndata(data, flavor='seurat', sample_id='PigP3', output='/scratch/user/sean.thompson/20241203_163132/12.3/h5ad/PigP3_bin20_normalized.stereo.forSeurat.h5ad')
#del adata#clear the object to save memory

#save the StereoExpData to anndata (Scanpy compatible)
# remember to set flavor as seurat
#adata = st.io.stereo_to_anndata(data, flavor='scanpy', sample_id='PigP3', output='/scratch/user/sean.thompson/20241203_163132/12.3/h5ad/PigP3_bin20_normalized.stereo.forScanpy.h5ad')
#st.io.write_h5ad(adata, use_raw=True, use_result=True, key_record=None, output='/scratch/user/sean.thompson/20241203_163132/12.3/h5ad/PigP3_bin20_normalized.stereo.forScanpy.h5ad', split_batches=True)
st.io.stereo_to_anndata(data, flavor='scanpy', sample_id='PigP3', output='/scratch/user/sean.thompson/20241203_163132/12.3/h5ad/PigP3_bin20_normalized.stereo.forScanpy.h5ad')
#del adata#clear the object to save memory


# In[22]


#load the h5ad
#data = st.io.read_h5ad(file_path='/scratch/user/sean.thompson/20241203_163132/12.3/PigP3_bin50.stereo.h5ad', flavor='stereopy')

#validate proper read-in
#data.plt.spatial_scatter_by_gene(gene_name='COL1A1', palette='linear_grey_10_95_c0',out_dpi=500, height=1000, width=1000)
#data.plt.cluster_scatter(res_key='leiden', groups=['1', '2', '3', '4', '5', '8', '9', '10',
#                                                  '11', '12', '13', '14', '15', '16', '17', 
#                                                   '18', '19', '20', '21', '22'])
#data.plt.cluster_scatter(res_key='leiden', groups=['1', '2', '3', '4', '5', '8', '9', '10',
#                                                  '11', '12', '13', '14', '15', '16', '17', 
#                                                   '18', '19', '20', '21', '22'], base_image='/scratch/user/sean.thompson/20241203_163132/12.2/PigP3_realign/Pig_P3_3/outs/image/B03703F1_ssDNA_regist.tif')


