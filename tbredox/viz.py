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
    sns.set_context('notebook', rc=rc)


def set_style_bokeh():

    '''
    Formats bokeh plotting environment similar
    to that used in Physical Biology of the Cell, 2nd edition.
    Based on @gchure and @mrazomej's work.

    '''

    theme_json = {'attrs':
            {'Figure': {
                'background_fill_color': '#ffffff',
                'outline_line_color': '#000000',
            },
            'Axis': {
            'axis_line_color': "slategray",
            'major_tick_in': 7,
            'major_tick_line_width': 2,
            'major_tick_line_color': "white",
            'minor_tick_line_color': "grey",
            'axis_label_text_font': 'Helvetica Neue',
            'axis_label_text_font_style': 'normal'
            },
            'Grid': {
                'grid_line_color': 'white',
            },
            'Legend': {
                'background_fill_color': '#E3DCD0',
                'border_line_color': 'slategray',
                'border_line_width': 1.5,
                'background_fill_alpha': 0.5
            },
            'Text': {
                'text_font_style': 'normal',
               'text_font': 'Helvetica'
            },
            'Title': {
                'background_fill_color': '#FFEDC0',
                'text_font_style': 'normal',
                'align': 'center',
                'text_font': 'Helvetica Neue',
                'offset': 2,
            }}}

    theme = bokeh.themes.Theme(json = theme_json)
    bokeh.io.curdoc().theme = theme 


    #return theme_json


def plot_sample_datashade(df, sample_name, vars_, sample_col = 'sample_id', **kwargs): 
    """
    Returns an hvplot object with the sample colored by an indicator 
    variable colored using datashading. This function is devised for 
    visualizing datasets with millions of datapoints. 

    Params
    ------
    df (pd.DataFrame)
        Annotated pandas dataframe. 

    sample_name (str)
        Name of the sample to be colored. 

    vars (list)
        Name of the xy variables for the scatter plot.

    sample_col(str, default = 'sample_id')
        Name of the column for which the sample_name will be selected from.

    kwargs 
        All kwargs go directly to format the hvplot object.  


    Returns 
    -------

    shader_plot ()
        Scatter plot colored by sample name using datashader. 

    """

    # Assert sample in sample_col
    assert sample_name in df[sample_col].values 

    # Assert there are more than two vars to plot with 

    assert len(vars_) >= 2 

    # Make binary indicator var
    indicator_variable = [1 if smpl== sample_name else 0 for smpl in df[sample_col]]

    # Add variable to dataframe 
    df[sample_name] = indicator_variable

    # Initialize plot for two variables  
    if len(vars_) == 2: 
        var_1, var_2  = vars_ 

        shader_plot = df.hvplot.scatter(
            x = var_1, 
            y = var_2, 
            c = sample_name, 
            width = 600,
            datashade = True, 
            **kwargs 
        )

    # Initialize plot for multiple variables 
    else : 

        shader_plot = df.hvplot.scatter(
            x = vars_[0], 
            y = vars_[1:], 
            c = sample_name, 
            datashade = True, 
            **kwargs 
        )
    return shader_plot 



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