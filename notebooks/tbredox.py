import numpy as np

# Interactive plotting
import bokeh_catplot
import bokeh.io
import holoviews as hv
import colorcet as cc
from bokeh.layouts import gridplot
import hvplot.pandas
import holoviews.operation.datashader
import datashader as ds

import warnings
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rcParams



rcParams['axes.titlepad'] = 20 


def set_plotting_style_plt():
      
    tw = 1.5
    rc = {'lines.linewidth': 2,
        'axes.labelsize': 18,
        'axes.titlesize': 21,
        'xtick.major' : 16,
        'ytick.major' : 16,
        'xtick.major.width': tw,
        'xtick.minor.width': tw,
        'ytick.major.width': tw,
        'ytick.minor.width': tw,
        'xtick.labelsize': 'large',
        'ytick.labelsize': 'large',
        'font.family': 'sans',
        'weight':'bold',
        'grid.linestyle': ':',
        'grid.linewidth': 1.5,
        'grid.color': '#ffffff',
        'mathtext.fontset': 'stixsans',
        'mathtext.sf': 'fantasy',
        'legend.frameon': True,
        'legend.fontsize': 12, 
       "xtick.direction": "in","ytick.direction": "in"}



    plt.rc('text.latex', preamble=r'\usepackage{sfmath}')
    plt.rc('mathtext', fontset='stixsans', sf='sans')
    sns.set_style('ticks', rc=rc)

    #sns.set_palette("colorblind", color_codes=True)
    sns.set_context('notebook', rc=rc)




def plot_genes_colombos_expression(adata, condition_str, list_of_genes, 
								   condition_col= 'annotation', gene_col = 'Rv_ID',
								   annot_col = 'Test description', condition_list = None): 

	"""
	Returns a KDE plot for each expression experiment matching a given pattern. 
	
	Params
	------
	adata (ad.AnnData)

	condition_str (str)
		Pattern for which to extract expression conditions. Examples : 'CHOLESTEROL',
		'pH', 'mut', etc. We recommend checking the expression conditions manually a priori
		to avoid plotting unrelevant datasets for your purposes. 

	list_of_genes (list)
		List of genes to extract data from. Rv_ID format is recommended.

	condition_col (str, default = 'annotation')
		Column in the adata.var dataframe to make the text matching with. 

	gene_col (str, default = 'Rv_ID')
		Column name in the adata.obs dataframe to extract genes from. 

	
	annot_col (str, default = 'Test description')
		Annotation column from the adata.obs dataframe to make the 

	condition_list(list, default = None) 
		List of expression conditions to extract data from. 
		This is an optional argument to use instead of the condition_str.

	"""

	n_genes = len(list_of_genes)

	# Get sample indices 

	sample_ixs = adata.var[adata.var[condition_col].str.contains(condition_str)].index.to_list()
	n_samples = len(sample_ixs)

	# Get fold change values for genes of interest 

	log_fc_genes = adata[adata.obs[gene_col].isin(list_of_genes), sample_ixs].X.A

	# Get fold-change distributions for all genes in those conditions
	log_fc_all = adata[:, sample_ixs].X.A

	# Initialize color palette
	pal = sns.color_palette('colorblind', n_colors = n_genes)

	fig = plt.figure(figsize = (5, n_samples * 4))

	fig.subplots_adjust(hspace = 1.5, wspace = 0.3)

	# Stack KDE plots per sample vertically
	for i in range(n_samples):

	    ax = fig.add_subplot(n_samples, 1, i+1)
	    
	    ax = sns.distplot(log_fc_all[:, i], color = 'grey')
	    
	    for j in range(n_genes):
	        ax.axvline(log_fc_genes[j, i], label = list_of_genes[j], linewidth = 3, color = pal[j])
	        ax.set_title(adata.var.loc[sample_ixs[i]][annot_col])
	    
	    ax.set_xlabel(r'$log_{10}(FC)$')
	    ax.set_xlim(-5, 5)
	    plt.legend()