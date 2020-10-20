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
