import umap
import anndata as ad
from sklearn.preprocessing import StandardScaler
import numpy as np 
import time
import matplotlib.pyplot as pyplot
import seaborn as sns


seed =385
np.random.seed(seed)

# Load dataset
adata = ad.read_h5ad('../data/tb_adata_v1.h5ad')

# Save raw data for later
adata.raw = adata

print('Loaded dataset. ')

# Scale all columns to a standard gaussian, mu = 0, sigma = 1
scaler = StandardScaler()

# Replace dataset with scaled version 
adata.X = scaler.fit_transform(adata.X)

print('Finished normalization ... \n')

# Extract datasets
adata_expression = adata[:, adata.var.data_type == 'expression']
adata_tn_bin = adata[:, adata.var.data_type == 'tnseq_bin']
adata_tn_fc = adata[:, adata.var.data_type == 'tnseq_fc']

# Also get data combinations 
adata_exp_bin = adata[:, adata.var.data_type.isin(['expression', 'tnseq_bin'])]
adata_exp_fc = adata[:, adata.var.data_type.isin(['expression', 'tnseq_fc'])]

t_start = time.time()

print('Starting UMAPs....')
# Get umap coordinates for each one 
u_exp = umap.UMAP(random_state = seed, min_dist = 0, n_components = 3).fit_transform(adata_expression.X)
u_bin= umap.UMAP(random_state = seed, min_dist = 0, n_components = 3).fit_transform(adata_tn_bin.X)
u_fc = umap.UMAP(random_state = seed, min_dist = 0, n_components = 3).fit_transform(adata_tn_fc.X)
u_expbin = umap.UMAP(random_state = seed, min_dist = 0, n_components = 3).fit_transform(adata_exp_bin.X)
u_expfc = umap.UMAP(random_state = seed, min_dist = 0, n_components = 3).fit_transform(adata_exp_fc.X)

t_end = time.time()

print(f'Computation of all UMAP coordinates took {(t_end - t_start)/60 } minutes.\n')


# Save and export ... 
adata.obsm['umap_exp'] = u_exp
adata.obsm['umap_tnseq_bin'] = u_bin
adata.obsm['umap_tnseq_fc'] = u_fc
adata.obsm['umap_exp_tnseq_bin'] = u_expbin
adata.obsm['umap_exp_tnseq_fc'] = u_expfc

print('Exported umap coordinates as .obsm objects....')

fname = 'tb_adata_v2.h5ad'
print('Saving to: ', fname)
adata.write('../data/' + fname)


