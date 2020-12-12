import matplotlib.pyplot as plt
import anndata as ad
from tbredox import viz
import colorcet as cc
import numpy as np
import pandas as pd
import warnings 
warnings.filterwarnings('ignore')
#%matplotlib inline
#%config inlinebackend.figure_format = 'retina'

viz.set_plotting_style_plt()


a = ad.read_h5ad('../data/tb_adata_v2.h5ad')
print('Data loaded. Making the cluster colors.')

# Make an array of hex values using colorcet

n_clus = a.obs.cluster_labels.max()
cmap_dict = dict(zip(np.arange(n_clus + 1), cc.glasbey[:n_clus + 1]))
cluster_colors = list(map(lambda x: cmap_dict[x], a.obs.cluster_labels.values))

plot_combinations = list(a.obsm.keys())

height = len(plot_combinations)
width = 1

print('Initializing plot...')
# Initialize figure and subplots
fig, ax = plt.subplots(height,width, figsize = (4, height * 4))

fig.subplots_adjust(hspace=0.80, wspace=.2)

#Make 2d scatter plots
for ix, combination in enumerate(plot_combinations): 
    ax[ix].scatter(
        a.obsm[combination][:, 0],
        a.obsm[combination][:, 1],
        alpha = 0.4,
        c = cluster_colors,
        s = 3
    )

    ax[ix].set_title(combination)
    ax[ix].axis('off')

plt.show()