import numpy as np
import scipy.stats as st
import pandas as pd 
from sklearn.neighbors import NearestNeighbors
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import tqdm


def get_neighbors(data, n_neighbors = 5, ix = None, **kwargs): 
    """
    Given an input dataframe with a given coordinate system, 
    the function returns the indices of the nearest neighbors
    of a dataset or a specific data point. 


    Params 
    ------
    data (pd.DataFrame)
        Dataframe whose columns represent a given coordinate system where distances
        can be calculated.

    n_neighbors(int, default = 5)
        Number of nearest neighbors. 
        
    ix (int, default = None)
        Index of data point for which to find the nearest neighbors.
        Defaults to None, to return the indices of all data points. 

    Returns 
    -------
    neighbor_ixs (np.array)
        Array of the indices of the nearest neighbors of shape (n_datapoints, n_neighbors).
    """

    # Initalize Nearest Neighbors object
    knn = NearestNeighbors(n_neighbors = n_neighbors).fit(data)

    neighbors = knn.kneighbors()

    neighbor_ixs = neighbors[1]

    # Loop to return ixs of a specific point.
    if ix is not None: 
        return neighbor_ixs[ix]

    return neighbor_ixs


def choose_k(data, max_clus):
    "Helper function to choose the number of clusters for KMeans."
    
    k_tests = np.arange(2, max_clus)
    wcss = []
    silhouettes_ = []
    for k in tqdm.tqdm(k_tests): 
        clustering = KMeans(k, random_state= 42).fit(data)
        wcss.append(clustering.inertia_)
        silhouettes_.append(silhouette_score(data, clustering.predict(data)))
        
    return wcss, silhouettes_


def fisher_enrichment_test(df_annot, annotation, cluster, clus_col_name = 'cluster_labels'): 
    """
    Returns a repor dataframe with the top 5 enriched functions
    for a given subset of data. This function is especially suited
    for statistical enrichment tests after clustering. 
    
    Params 
    ------
    df_annot (pd.DataFrame)
        Annotated dataframe containing the 'annotation' column 
        and a 'clus_col_name' column. 
    
    annotation (str)
        Annotation to make the enrichment test on. In the case 
        of gene set enrichment this could be a Gene Ontology 
        or COG annotation. 
    
    cluster (int or str)
        Cluster (or in general group of data points) to test.
    
    col_clus_name (str)
        Name of the cluster column in the df_annot dataframe.
    
    Returns 
    -------
    df_report (pd.DataFrame)
        Report dataframe with pvalues and annotation names. 
    
    """
    # Get subset of completely annotated genes 
    df_test = df_annot[pd.notnull(df_annot[annotation])]

    # Number of genes with valid annotation 
    M = df_test.shape[0]

    # Extract data for given cluster
    df_clus = df_test[df_test[clus_col_name] == cluster]

    # Get top 5 categories to test
    cats = df_clus[annotation].value_counts().head().index.to_list()
    
    # Number of genes in the cluster (sample size)
    N = df_clus.shape[0]
    
    # Initialize pvalue array 
    pvals = np.empty(len(cats))
    
    # Loop through the top categories
    for i, cat in enumerate(cats): 
        
        df_cat = df_test[df_test[annotation] == cat]
        
        # Total number of genes that map to given category (total number of white balls)
        n = df_cat.shape[0]
        
        # Number of genes inside cluster that map to given category (number of white balls in sample)
        x = df_clus[df_clus[annotation] == cat].shape[0]
        
        # Sweep through the probabilities from x to n 
        pmfs = st.hypergeom.pmf(k = np.arange(x, n + 1), N = N, n = n, M = M)
        
        # Compute pval
        pvals[i] = pmfs.sum()
    
    # Save results
    df_report = pd.DataFrame(
        {'categories': cats, 'pval': pvals}
    )

    df_report['cluster'] = cluster
    df_report['annot'] = annotation

    return df_report 


def get_missing_columns(data):    
    '''
    Returns a dataframe of the missing values in each column 
    with its corresponding dtype.
    
    Params 
    ------
    data (pd.DataFrame)
        Input dataframe to be tested. 
        
    Returns 
    --------
    missing_cols_df(pd.DataFrame)
        Report dataframe with % of missing values. 
    '''
    
    # Generate a DataFrame with the % of missing values for each column
    df_missing_values = (data.isnull().sum(axis = 0) / len(data) * 100)\
                        .sort_values(ascending = False)\
                        .to_frame('% missing_values').reset_index()
    
    # Generate a DataFrame that indicated the data type for each column
    df_feature_type = data.dtypes.to_frame('feature_type').reset_index()
    
    # Merge frames
    missing_cols_df = pd.merge(df_feature_type, df_missing_values, on = 'index',
                         how = 'inner')

    missing_cols_df.sort_values(['% missing_values', 'feature_type'], inplace = True)
    
    return missing_cols_df