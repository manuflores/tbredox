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
    Returns a report dataframe with the top 5 enriched functions
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


def confusion_matrix(pred_labels, true_labels): 
    """
    Returns a confusion matrix from a multiclass classification 
    set of labels. 
    
    Params 
    ------
    pred_labels (array-like):
        List of labels as predicted by a classification algorithm. 
    
    true_labels (array-like): 
        List of ground truth labels. 
    
    Returns 
    -------
    conf_mat (array-like): 
        Confusion matrix. 
    """
    
    n_labels = max(np.max(pred_labels), np.max(true_labels)) + 1
    conf_mat = np.zeros(shape = (n_labels, n_labels))
    
    for (i, j) in zip(pred_labels, true_labels): 
        conf_mat[i,j] +=1
        
    return conf_mat



def quantile_norm(X):
    """Normalize the columns of X to each have the same distribution.

    Given an expression matrix (microarray data, read counts, etc) of M genes
    by N samples, quantile normalization ensures all samples have the same
    spread of data (by construction).

    The data across each row are averaged to obtain an average column. Each
    column quantile is replaced with the corresponding quantile of the average
    column.

    Parameters
    ----------
    X : 2D array of float, shape (M, N)
        The input data, with M rows (genes/features) and N columns (samples).

    Returns
    -------
    Xn : 2D array of float, shape (M, N)
        The normalized data.
    """
    # compute the quantiles
    quantiles = np.mean(np.sort(X, axis=0), axis=1)

    # compute the column-wise ranks. Each observation is replaced with its
    # rank in that column: the smallest observation is replaced by 1, the
    # second-smallest by 2, ..., and the largest by M, the number of rows.
    ranks = np.apply_along_axis(st.rankdata, 0, X)

    # convert ranks to integer indices from 0 to M-1
    rank_indices = ranks.astype(int) - 1

    # index the quantiles for each rank with the ranks matrix
    Xn = quantiles[rank_indices]

    return(Xn)


def element_wise_entropy(px):
    """
    Returns a numpy array with element wise entropy calculated as -pi*log_2(p_i).
    
    Params
    ------
    px (np.array)
        Array of individual probabilities, i.e. a probability vector or distribution.
    
    Returns
    -------
    entropy (np.array)
        Array of element-wise entropies.
    """
    if isinstance(px, list):
        px = np.array(px)
        
    # Make a copy of input array
    entropy = px.copy()
    
    # Get indices of nonzero probability values
    nz = np.nonzero(entropy)
    
    # Compute -pi*log_2(p_i) element-wise
    entropy[nz] *= - np.log2(entropy[nz])
    
    return entropy


def entropy(ps): 
    "Returns the entropy of a probability distribution `ps`."
    # Get nonzero indices
    nz = np.nonzero(ps)

    # Compute entropy for nonzero indices
    entropy = np.sum(-ps[nz]*np.log2(ps[nz]))
    
    return entropy