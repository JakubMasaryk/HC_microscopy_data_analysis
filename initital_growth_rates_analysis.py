# ## __Analysis of the Average OD600 Change in the Early Exponential Phase__

# #### __general description__
# * __growth assays__ on selected mutant strain: __temperature-sensitive mutants__ (TS collection) in presence of __arsenic__
# * focus on the __early exponential phase__: app. __5-15 hours of growth__
# * __linear regression__ line fitted to the data points (characterozed by a timepoint (h) and OD600) between 5-15 hours and __slope calculated__ (__slope__ eqauls and __average OD600 increase per hour__)
# * 4D __data__ (slope control, slope exposed, p-value control and p-value exposed) also __clustered__ by 2 algorithms: __K-Means and DBSCAN__

# #### __input data__
# * __OD measurements datasets__: well-timepoint-OD600
# * __plate layout dataset:__ well-strain-repeat-As concentration
# * __load from a file__: __specify__ file pathways in __'file pathways'__ section and __define__ in __'load_process_merge' function__ (also choose __'files'__ for source)
# * __load from a cloud__: __specify__ file labels in __'file labels'__ (pre-defined) section and __define__ in __'load_process_merge' function__ (also choose __'cloud'__ for source)

# #### __inputs__

# * __define the temporal range__ in the __'analyse slopes' section__ as an argument of 'analyse_slopes' function (pre-set to 5-15 hours)

# ## __Libraries__
import numpy as np
import pandas as pd
import seaborn as sns
from b2sdk.v2 import InMemoryAccountInfo, B2Api
from io import BytesIO
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import ttest_ind
from sklearn.cluster import KMeans
from sklearn.cluster import DBSCAN
from sklearn.metrics import silhouette_score as ss


# ## __Cloud Access__

# #### __Backlblaze B2 authentication__
data_bucket_name= 'bioscreen-data'
bucket_key_id= '003b5f880f95dd40000000009';
bucket_key= 'K003+OGX7qB+jAEv5GZrF93mHKj1NcU'


# ## __Functions__

# #### __statistics__
#margin of error: t-distribution, CL 95%, confidence intervals: mean +/- margin of error
#input: list of OD600 values from individual biological replicates (for a single timepoint)
def t_margin_of_error_cl95(data):
    cl=0.95
    std= np.std(data, ddof=1)
    n=len(data)
    
    std_err= std/np.sqrt(n)
    t_score = stats.t.ppf((1 + cl) / 2, df=n - 1)  
    t_margin_of_error = t_score * std_err
    return t_margin_of_error

#t-test: independent (two sample), assumes equal variance (by default), two-tailed
def single_t_test(column1, column2):
    try:
        # Run t-test
        t_stat, p_value = ttest_ind(column1, column2)
        return p_value
    except Exception:
        return np.nan
    
#calculate slope from 2d data: x- time, y- OD600
def slope(x, y):
    try:
        slope, _ = np.polyfit(x, y, 1)
        return slope
    except Exception:
        return np.NaN
    
#0-1 normalization of the numerical columns (prep for clustering)
def normalize_numerical_columns(dataset):  
    #extract numerical columns
    numeric_cols = dataset.select_dtypes(include="number").columns
    #normalize each column
    for col in numeric_cols:
        #extract min and max values
        col_min = dataset[col].min()
        col_max = dataset[col].max()
        #constant columns (min=max)
        if col_max == col_min:
            dataset[col] = 0
        #0-1 norm
        else:
            dataset[col] = (dataset[col] - col_min) / (col_max - col_min)

    return dataset


#k-means clustering
#input NORMALIZED (use 'normalize_numerical_columns') dataset for clustering
def k_means(data, k= 3):
    #filter down to numerical columns
    data= data.select_dtypes(include='number')
    
    try:
        kmeans_clustering= KMeans(n_clusters= k)
        kmeans_clustering.fit(data)
        clustered_data= data.assign(clusters= kmeans_clustering.labels_)
        print('clustering succesful')
        return clustered_data
    except Exception as ex:
        print(f'clustering failed: {ex}')
        return data.assign(k_means_clusters= 0)
    
    
#dbscan
#input NORMALIZED (use 'normalize_numerical_columns') dataset for clustering
def dbscan(data, epsilon= 0.5, no_of_points= 5):
    #filter down to numerical columns
    data= data.select_dtypes(include='number')
    
    try:
        dbscan_clustering= DBSCAN(eps= epsilon, min_samples= no_of_points)
        dbscan_clustering.fit(data)
        clustered_data= data.assign(clusters= dbscan_clustering.labels_)
        print('clustering succesful')
        return clustered_data
    except Exception as ex:
        print(f'clustering failed: {ex}')
        return data.assign(dbscan_clusters= 0)
    
    
#performs elbow plot by performing a k-means clustering using multiple k values
#calculates and visualizes WCSS for each k value
#input NORMALIZED (use 'normalize_numerical_columns') dataset for clustering
def elbow_plot(data, min_k= 1, max_k= 10):
    
    #filter down to numerical columns
    data= data.select_dtypes(include='number')
    #remove the WT row
    data= data.loc[data.index != 'WT']
    #listo to store wcss values
    wcss = []
    ks= []
    #define range of teste k
    k_range = range(min_k, max_k+1)

    for k in k_range:
        try:
            kmeans_clustering= KMeans(n_clusters= k)
            kmeans_clustering.fit(data)
            wcss.append(kmeans_clustering.inertia_)
            ks.append(k)
        
        except Exception as ex:
            print(f'clustering failed: {ex}\n k value {k} skipped')
        
    fig, ax= plt.subplots()   
    ax.plot(ks,
            wcss,
            marker='o',
            mfc= 'white',
            mew= 2.5,
            lw= 2.5,
            color= 'red')
    ax.set_xticks(ks)
    ax.set_xlabel('k', weight= 'bold', fontsize= 13)
    ax.set_ylabel('WCSS', weight= 'bold', fontsize= 13)
    

#evaluates DBSCAN parameters based on the silhouette score
#input NORMALIZED (use 'normalize_numerical_columns') dataset for clustering
def DBSCAN_param_opt(data, no_of_combinations= 20, min_eps=0.1, max_eps=1, min_min_samples=3, max_min_samples=10):
    
    #filter down to numerical columns
    data= data.select_dtypes(include='number')
    #remove the WT row
    data= data.loc[data.index != 'WT']
    
    #parameter combinations
    try:
        epsilons= np.linspace(min_eps, max_eps, no_of_combinations)
        min_samples= np.linspace(min_min_samples, max_min_samples, no_of_combinations)
        combos= [[round(a, 2), round(b)] for a in epsilons for b in min_samples]
    except Exception as ex:
        raise RuntimeError(f"unable to combine DBSCAN parameters: {ex}")
    
    #extracting the silhouette scores
    sil_scores= []
    for combo in combos:
        try:
            #dbscan
            eps, min_sam= combo[0], combo[1]
            dbscan= DBSCAN(eps=eps, min_samples=min_sam)
            dbscan.fit(data)
            #no. of clusters
            clusters= set(dbscan.labels_)
            clusters.discard(-1)
            no_of_clusters= len(clusters)
            #conditional append (at least 2 clusters must be present (excluding outliers) to calculate the ss)
            if no_of_clusters >1:
                sil_scores.append(ss(data, dbscan.labels_))
            else:
                sil_scores.append(np.NaN)
        except Exception as ex:
            sil_scores.append(np.NaN)
            print(f'combination epsilon {eps} and min. samples {min_sam} skipped: {ex}')
    
    #storing ss in DF        
    combos_scores= pd.DataFrame({'Combo': combos, 
                                 'SilhouetteScore':sil_scores}).dropna()
    
    return combos_scores.sort_values('SilhouetteScore', ascending= False)


# #### __data processing__
### load both the bioscreen data and corresponding plate layout, process and merge
# load from files on local hardrive (source= 'files'), file pathways need to be defined
#load from Backblaze B2 cloud storage (source= 'cloud'), file names (as listed on the cloud) need to be defined
def load_process_merge(bioscreen_data_pathway= None,
                       plate_layout_pathway= None,
                       bioscreen_data_cloud_label= None,
                       plate_layout_cloud_label= None,
                       source= 'files',
                       bucket_name= data_bucket_name, 
                       key_id= bucket_key_id, 
                       key= bucket_key):
    
    if source== 'files':
        #bioscreen-data file load
        try:
            data= pd.read_csv(bioscreen_data_pathway,
                              header= 2)
        except FileNotFoundError:
            raise
        except Exception as ex:
            raise RuntimeError(f"Failed to load bioscreen-data file: {ex}")

        #plate-layout file load
        try:
            layout= pd.read_excel(plate_layout_pathway)
        except FileNotFoundError:
            raise
        except Exception as ex:
            raise RuntimeError(f"Failed to load plate-layout file: {ex}")
    
    elif source== 'cloud':
        #access the cloud
        try:
            #init. B2 SDK
            info = InMemoryAccountInfo()
            b2_api = B2Api(info)
            # authentication
            application_key_id = key_id
            application_key = key
            b2_api.authorize_account("production", application_key_id, application_key)
            # get the bucket
            bucket = b2_api.get_bucket_by_name(bucket_name)
        except Exception as ex:
            raise RuntimeError(f'Failed to access the cloud storage: {ex}')
        #load the data
        try:
            # all_file_names= [file_info.file_name for file_info, _ in bucket.ls(recursive= True)]
            # return all_file_names
            ###data load from Backblaze B2###
            #store data into in-memory object ('in_memory_data')
            in_memory_data = BytesIO() #inmemory storage object 1
            in_memory_plate_layout = BytesIO() #inmemory storage object 2
            
            bucket.download_file_by_name('temperature-sensitive-alleles/' + bioscreen_data_cloud_label).save(in_memory_data) #download data from specified bucket-subbucket ('temperature-sensitive-alleles/')-file into 'in_memory_storage'
            bucket.download_file_by_name('temperature-sensitive-alleles/' + plate_layout_cloud_label).save(in_memory_plate_layout) #download plate layout from specified bucket-subbucket ('temperature-sensitive-alleles/')-file into 'in_memory_storage'
            in_memory_data.seek(0) #rewind to the beginning
            in_memory_plate_layout.seek(0) #rewind to the beginning
            # Load into pandas df
            data = pd.read_csv(in_memory_data,
                               header= 2)
            layout= pd.read_excel(in_memory_plate_layout)

        except Exception as ex:
            raise RuntimeError(f'Failed to access the data: {ex}')
            
    else:
        raise ValueError(f"Invalid source argument: '{source}'. Expected: 'files' or 'cloud'.")
        
    #reorganize the bioscreen dataset (drop blank, melt, convert time to hours and unify the dtype of well with layout)
    try:
        data= data.drop(columns= 'Blank', 
                        errors='ignore')
        data= data.melt(id_vars= 'Time',
                        var_name= 'Well',
                        value_name= 'OD600')  
        data= data.assign(Hours= data.Time.apply(lambda x: x.split(':')[0]),
                          Minutes= data.Time.apply(lambda x: x.split(':')[1])) #split the time column to h and m
        data= data.astype({'Well':'Int32', 'OD600':'Float32', 'Hours':'Int32', 'Minutes':'Int32'}) #adjust dtypes
        data= data.assign(Hours= data.Hours + data.Minutes/60) #sum the hour timepoint
        data= data.drop(columns= ['Time', 'Minutes'])
    except Exception as ex:
        raise RuntimeError(f'Failed to transform the bioscreen dataset: {ex}')

    #merge to the layout dataset
    try:
        data= data.merge(layout,
                         how= 'inner',
                         on= 'Well')
    except Exception as ex:
        raise RuntimeError(f'Failed to merge the bioscreen and the layout datasets: {ex}')

    #reorder columns
    try:
        data= data.reindex(columns= ['Well', 'Hours', 'Strain', 'BiologicalRepeat', 'TechnicalRepeat', 'AsConcentration', 'OD600'])
    except Exception as ex:
        raise RuntimeError(f'Failed to reorder columns of the final dataset: {ex}')

    return data


###averaging repeats (both technical and biological) + basic stats
def tech_rep_average(dataset):
    
    #average the technical repeats
    try:
        data= dataset.groupby(['Strain', 'AsConcentration', 'Hours', 'BiologicalRepeat'])[['OD600']].mean().reset_index()
        data= data.sort_values(['AsConcentration', 'Strain', 'BiologicalRepeat', 'Hours'])
    except Exception as ex:
        raise RuntimeError(f'Failed to average the technical repeats: {ex}')
        
    return data


###drop technical replicates with low-quality data
def drop_wells(dataset, wells_to_drop):
    if not wells_to_drop:
        return dataset
    else:
        return dataset.loc[~dataset.Well.isin(wells_to_drop)]
    

#drop strain (index) from data fro cliustering
#input: data for clustering + list of alleles (mutants) to drop
#mostly for apparent outliers
def drop_strain(dataset_for_clustering, mutants_to_drop):
    if not mutants_to_drop:
        return dataset_for_clustering
    else:
        return dataset_for_clustering.loc[~dataset_for_clustering.index.isin(mutants_to_drop)]


###calculating slopes for specific time range (each strain-condition-repeat), group by strain-condition, store individual slope valuesin list and: mean, std, moe + ttest WT vs mutants for each condition
#for controls (each condition-repeat) slope calculated from 4 sets of datapoints (one for each experiment)
#input: complete dataset with averaged technical replicates
def analyse_slopes(dataset, start, end):
    #define timerange
    data= dataset.loc[(dataset.Hours>=start)&
                      (dataset.Hours<=end)]  

    try:
        #extract slopes for each strain-condition-repeat using the 'slope' function
        data= data.groupby(['Strain', 'AsConcentration', 'BiologicalRepeat']).apply(lambda x: slope(x['Hours'], x['OD600']), include_groups=False).reset_index()
        #group the repeats for each strain-condition into a list  (0- temporary label for the slope column)
        data= data.groupby(['Strain', 'AsConcentration']).agg({0: list}).reset_index()
        data.columns= ['Strain', 'AsConcentration', 'Slopes']
        #mean, std and moe of slopes
        data= data.assign(SlopesMean= data.Slopes.apply(lambda x: np.array(x).mean()),
                          SlopesSTD= data.Slopes.apply(lambda x: np.array(x).std(ddof= 1)),
                          SlopesMOE95= data.Slopes.apply(lambda x: t_margin_of_error_cl95(x)))
        #control data (+ relabel the slopes column to 'ControlSlopes')
        control= data.loc[data.Strain== 'WT', ['AsConcentration', 'Slopes']]
        control.columns= ['AsConcentration', 'ControlSlopes']
        #merge control dataset to the original dataset based on condition ('AsConcentration')
        data= data.merge(control, how= 'left', on= 'AsConcentration')
        #ttest slopes to the correspondng controlslopes
        data= data.assign(SlopesPValue= data.apply(lambda x: single_t_test(x['Slopes'], x['ControlSlopes']), axis= 1))
        #drop control slopes list
        data= data.drop(columns=['Slopes', 'ControlSlopes'])
        #evaluate significance 
        data= data.assign(Significance= np.where(data.SlopesPValue >= 0.05, '',
                                                 np.where(data.SlopesPValue >= 0.01, '*',
                                                          np.where(data.SlopesPValue >= 0.001, '**', '***'))))
        
        return data
    except Exception as ex:
        raise RuntimeError(f'Failed to analyse the slopes: {ex}')
        

#adds a new column 'Cluster'
#assign a cluster name to each mutant that corresponds to the image-based screening labels
def assign_cluster(dataset):
    
    #cluster
    actin_cluster= ['act1-105', 'act1-108', 'act1-119', 'act1-121', 'act1-124', 'act1-133', 'act1-155', 'act1-3', 
                'arc15-10', 'arc35-6', 'arp2-14', 'las17-14', 'last17-1', 'pfy1-14', 'pfy1-4', 'srv2-ts',]
    tubulin_cluster= ['spc110-221', 'spc34 41-1', 'tub4-Y445D', 'tub4-ΔDSY', 'ask1-2', 'dam1-19']
    vesicular_transport_cluster= ['bos1-1', 'sec12-1', 'sec13-1', 'sec17-1', 'sec23-1', 'sed5-1', 'yip1-4', 'yip1-40', 'ypt1-3']
    proteasome_cluster= ['pre2-127', 'rpn1-821', 'rpn6-1', 'rpn7-3', 'rpt3-1', 'rpt6-20']
    nucleus_import_cluster= ['gsp1-P162L',  'nup57-E17', 'srm1-ts', ]
    
    try:
        dataset = dataset.assign(
            Cluster=np.where(
                dataset.Strain.isin(actin_cluster), 'actin',
                np.where(
                    dataset.Strain.isin(tubulin_cluster), 'tubulin',
                    np.where(
                        dataset.Strain.isin(vesicular_transport_cluster), 'vesicular transport',
                        np.where(
                            dataset.Strain.isin(proteasome_cluster), 'proteasome 3',
                            np.where(
                                dataset.Strain.isin(nucleus_import_cluster),
                                'nuclear import',
                                'CTRL'))))))
        return dataset.sort_values(['Cluster', 'Strain', 'AsConcentration'])
         
    except Exception as ex:
        raise RuntimeError(f'Failed to assign cluster: {ex}')


# #### __plotting__
def plot(data, x_offset= 0.15, export= False, path_for_export= r"C:\Users\Jakub\Desktop\fig.png"):
    
    try:
        fig, ax= plt.subplots(figsize= (20, 9.6))

        #plot params
        width= .3
        x= np.arange(0, len(data.Strain.unique()))
        
        #split the data by condition
        control= data.loc[data.AsConcentration== 0].reset_index()
        _0_5_mM= data.loc[data.AsConcentration== 0.5]
        _1_mM= data.loc[data.AsConcentration== 1]
        
        #cluster-based coloring
        colors= ['#FEF247' if cluster== 'actin' else 
                 '#85CC00' if cluster== 'tubulin' else
                 '#ED872D' if cluster== 'vesicular transport' else
                 '#FF4646' if cluster== 'proteasome 3' else
                 '#DBF796' if cluster== 'nuclear import' else
                 'Gray' for cluster in control.Cluster]
        
        #plot
        ax.bar(x-width,
               control.SlopesMean,
               yerr= control.SlopesMOE95,
               color= colors,
               width= width,            
               edgecolor= '#00527C',
               capsize= 2,
               error_kw= {'elinewidth':.75},
               lw= 1)
        ax.bar(x,
               _0_5_mM.SlopesMean,
               yerr= _0_5_mM.SlopesMOE95,
               color= colors,
               alpha= 0.5,
               width= width,            
               edgecolor= '#00527C',
               capsize= 2,
               error_kw= {'elinewidth':.75},
               lw= 1)
        ax.bar(x+width,
               _1_mM.SlopesMean,
               yerr= _1_mM.SlopesMOE95,
               color= colors,
               alpha= 0.2,
               width= width,            
               edgecolor= '#00527C',
               capsize= 2,
               error_kw= {'elinewidth':.75},
               lw= 1)

        #significance symbols
        #control
        x_coor_ctrl= x - width + x_offset
        y_coor_ctrl= control.SlopesMean + control.SlopesMOE95 + 0.0025
        coordinates_ctrl= [[x, y] for x, y in zip(x_coor_ctrl, y_coor_ctrl)]
        for i, c in enumerate(coordinates_ctrl):
            cx,cy = c[0], c[1]
            ax.text(cx, cy, f'{control.Significance.iloc[i]}', fontsize= 13, ha= 'center', rotation = 90, weight= 'bold')
            
        #0.5mM As
        x_coor_05= x + x_offset
        y_coor_05= _0_5_mM.SlopesMean + _0_5_mM.SlopesMOE95 + 0.0025
        coordinates_05= [[x, y] for x, y in zip(x_coor_05, y_coor_05)]
        for i, c in enumerate(coordinates_05):
            cx,cy = c[0], c[1]
            ax.text(cx, cy, f'{_0_5_mM.Significance.iloc[i]}', fontsize= 13, ha= 'center', rotation = 90, weight= 'bold')
            
        #1mM As
        x_coor_1= x + width + x_offset
        y_coor_1= _1_mM.SlopesMean + _1_mM.SlopesMOE95 + 0.0025
        coordinates_1= [[x, y] for x, y in zip(x_coor_1, y_coor_1)]
        for i, c in enumerate(coordinates_1):
            cx,cy = c[0], c[1]
            ax.text(cx, cy, f'{_1_mM.Significance.iloc[i]}', fontsize= 13, ha= 'center', rotation = 90, weight= 'bold')
            
        ax.set_xticks(x)
        ax.set_xticklabels(data.Strain.unique(), rotation= 45, fontsize= 12)
        ax.set_xlabel('Strain', weight= 'bold', fontsize= 12)
        ax.set_ylabel('OD$_{600}$·h$^{-1}$', weight= 'bold', fontsize= 15)
        ax.set_ylim(0, 0.14)
        
        #export
        if export== True:
            plt.savefig(path_for_export, dpi= 1000)
        elif export== False:
            pass;
        else:
            raise ValueError(f"Invalid export argument: '{export}'. Expected: boolean ('True' or 'False').")
            
    except Exception as ex:
        print(f'failed to visualise {ex}')
        
    
#visualize clusters in a 2D scatter plot
#input: clutered data, p-values exposed and p-vaues non-exposed (from data_for_clustering- non normalized)
#datapoint size according to significance category ('-',' *', '**', '***')
#visuals: datapoint representing alleles too different from control under control conditions ('**', '***') lowered alpha
#         datapoint size according to significance category ('-',' *', '**', '***') exposed conditions, more significant difference-larger datapoint
def visualize_clusters(clustered_data, nonexposed_pvalue_series, exposed_p_value_series, export= False, path_for_export= r"C:\Users\Jakub\Desktop\fig.png"):
    try:
        #plot
        fig, ax= plt.subplots(figsize= (10, 10))
        
        #mutant (clustered) data
        ax.scatter(clustered_data.iloc[:, 0],
                   clustered_data.iloc[:, 1],
                   c= clustered_data.clusters,
                   s= [1000 if p_value < 0.001 else 500 if p_value < 0.01 else 250 if p_value < 0.05 else 80 for p_value in exposed_p_value_series],
                   alpha= [0.1 if p_value < 0.01 else 0.75 for p_value in nonexposed_pvalue_series])
        ax.set_ylabel('exposed slope', weight= 'bold', fontsize= 13)
        ax.set_xlabel('non-exposed slope', weight= 'bold', fontsize= 13)
        #labels without WT
        clustered_data_mut= clustered_data.loc[clustered_data.index != 'WT']
        for xi, yi, label in zip(clustered_data_mut.iloc[:, 0], clustered_data_mut.iloc[:, 1], clustered_data_mut.index):
            ax.text(xi, yi-0.025, label, fontsize= 9, ha='center', va='top', weight= 'bold', alpha= 0.33)
        #WT label
        clustered_data_wt= clustered_data.loc[clustered_data.index == 'WT']
        x= clustered_data_wt.iloc[0, 0]
        y= clustered_data_wt.iloc[0, 1]
        ax.text(x-0.01, y, 'WT', fontsize= 15 , ha='right', va='top', weight= 'bold')
        #export
        if export== True:
            plt.savefig(path_for_export, dpi= 1000)
        elif export== False:
            pass;
        else:
            raise ValueError(f"Invalid export argument: '{export}'. Expected: boolean ('True' or 'False').")
    except Exception as ex:
        print(f'cluster-visualisation failed: {ex}')


# ## __Data Preparation and Statistical Analysis__

# __file pathways__
#bioscreen 1
data_path_bsc1_05_mM= r""
layout_path_bsc1_05_mM= r""
data_path_bsc1_1_mM= r""
layout_path_bsc1_1_mM= r""

#bioscreen 2
data_path_bsc2_05_mM= r""
layout_path_bsc2_05_mM= r""
data_path_bsc2_1_mM= r""
layout_path_bsc2_1_mM= r""

#bioscreen 3
data_path_bsc3_05_mM= r""
layout_path_bsc3_05_mM= r""
data_path_bsc3_1_mM= r""
layout_path_bsc3_1_mM= r""

#bioscreen 4
data_path_bsc4_05_mM= r""
layout_path_bsc4_05_mM= r""
data_path_bsc4_1_mM= r""
layout_path_bsc4_1_mM= r""


# __file labels__ 
# * for cloud
#bioscreen 1
data_label_bsc1_05_mM= '20251210_bioscreen1_0.5mM.csv'
layout_label_bsc1_05_mM= '20251210_bioscreen1_0.5mM_plate_layout.xlsx'
data_label_bsc1_1_mM= '20251210_bioscreen1_1mM.csv'
layout_label_bsc1_1_mM= '20251210_bioscreen1_1mM_plate_layout.xlsx'

#bioscreen 2
data_label_bsc2_05_mM= '20251212_bioscreen2_0.5mM.csv'
layout_label_bsc2_05_mM= '20251212_bioscreen2_0.5mM_plate_layout.xlsx'
data_label_bsc2_1_mM= '20251212_bioscreen2_1mM.csv'
layout_label_bsc2_1_mM= '20251212_bioscreen2_1mM_plate_layout.xlsx'

#bioscreen 3
data_label_bsc3_05_mM= '20251215_bioscreen3_0.5mM.csv'
layout_label_bsc3_05_mM= '20251215_bioscreen3_0.5mM_plate_layout.xlsx'
data_label_bsc3_1_mM= '20251215_bioscreen3_1mM.csv'
layout_label_bsc3_1_mM= '20251215_bioscreen3_1mM_plate_layout.xlsx'

#bioscreen 4
data_label_bsc4_05_mM= '20251217_bioscreen4_0.5mM.csv'
layout_label_bsc4_05_mM= '20251217_bioscreen4_0.5mM_plate_layout.xlsx'
data_label_bsc4_1_mM= '20251217_bioscreen4_1mM.csv'
layout_label_bsc4_1_mM= '20251217_bioscreen4_1mM_plate_layout.xlsx'


# __load and merge the datasets__
#bioscreen 1
data_05_1= load_process_merge(bioscreen_data_cloud_label= data_label_bsc1_05_mM,
                              plate_layout_cloud_label= layout_label_bsc1_05_mM,
                              source= 'cloud')
data_1_1= load_process_merge(bioscreen_data_cloud_label= data_label_bsc1_1_mM,
                             plate_layout_cloud_label= layout_label_bsc1_1_mM,
                             source= 'cloud')

#bioscreen 2
data_05_2= load_process_merge(bioscreen_data_cloud_label= data_label_bsc2_05_mM,
                              plate_layout_cloud_label= layout_label_bsc2_05_mM,
                              source= 'cloud')
data_1_2= load_process_merge(bioscreen_data_cloud_label= data_label_bsc2_1_mM,
                             plate_layout_cloud_label= layout_label_bsc2_1_mM,
                             source= 'cloud')

#bioscreen 3
data_05_3= load_process_merge(bioscreen_data_cloud_label= data_label_bsc3_05_mM,
                              plate_layout_cloud_label= layout_label_bsc3_05_mM,
                              source= 'cloud')
data_1_3= load_process_merge(bioscreen_data_cloud_label= data_label_bsc3_1_mM,
                             plate_layout_cloud_label= layout_label_bsc3_1_mM,
                             source= 'cloud')

#bioscreen 4
data_05_4= load_process_merge(bioscreen_data_cloud_label= data_label_bsc4_05_mM,
                              plate_layout_cloud_label= layout_label_bsc4_05_mM,
                              source= 'cloud')
data_1_4= load_process_merge(bioscreen_data_cloud_label= data_label_bsc4_1_mM,
                             plate_layout_cloud_label= layout_label_bsc4_1_mM,
                             source= 'cloud')


# __drop low-quality technical replicates__
#bioscreen 1
data_1_1= drop_wells(data_1_1, [1, 11])

#bioscreen 2
data_1_2= drop_wells(data_1_2, [62])

#bioscreen 3
data_05_3= drop_wells(data_05_3, [19])
data_1_3= drop_wells(data_1_3, [6, 11, 12, 16, 23, 25, 27, 39, 56])

#bioscreen 4
data_1_4= drop_wells(data_1_4, [64])


# __average the technical replicates__
#bioscreen 1
data_05_1= tech_rep_average(data_05_1)
data_1_1= tech_rep_average(data_1_1)

#bioscreen 2
data_05_2= tech_rep_average(data_05_2)
data_1_2= tech_rep_average(data_1_2)

#bioscreen 3
data_05_3= tech_rep_average(data_05_3)
data_1_3= tech_rep_average(data_1_3)

#bioscreen 4
data_05_4= tech_rep_average(data_05_4)
data_1_4= tech_rep_average(data_1_4)


# __concatenate into a single dataset__
data= pd.concat([data_05_1, data_1_1,
                 data_05_2, data_1_2,
                 data_05_3, data_1_3,
                 data_05_4, data_1_4],
                 axis= 0)


# __analyse slopes__
data= analyse_slopes(data, 5, 15)


# __assign cluster from the image-based screening__
data= assign_cluster(data)


# __remove unnecessary clusters__
data= data.loc[data.Cluster != 'nuclear import']


# ## __Visualisation__

# __full dataset__
plot(data)


# __filtered dataset: based on p-value category in control (non-exposed) cells ('no significance' and '*')__
_filter= data.loc[(data.AsConcentration== 0) &
                  (data.Significance.isin(['', '*'])),
                  'Strain']
filtered_data= data.loc[data.Strain.isin(_filter)]

plot(filtered_data,
     x_offset= 0.05)


# ## __clustering__

# * 4D __data__ (control slope, exposed slope, p-value control conditions and p-value exposed cells) __clustered__ by __K-Means__ and __DBSCAN__
# * data __normalized__ (0-1)
# * clustering __parameters optimized__ by __Elbow Plot__ (K-Means) or __Silhouette Score__ (DBSCAN)
# * __input__ dataset: 'data', __final dataset__ from __previous section__

# #### __data prep__
#pivot around the As concentration
data_for_clustering= data.pivot_table(index= 'Strain',
                                      columns= 'AsConcentration',
                                      values= ['SlopesMean', 'SlopesPValue']).reset_index()
#relabel the columns
data_for_clustering.columns= ['strain', 'slope_0mM', 'slope_05mM', 'slope_1mM', 'slope_pvalue_0mM', 'slope_pvalue_05mM', 'slope_pvalue_1mM']
#set strain as index
data_for_clustering= data_for_clustering.set_index('strain')


# __drop strains__
#drop strain (index) from data fro cliustering
#input: data for clustering + list of alleles (mutants) to drop
#mostly for apparent outliers
def drop_strain(dataset_for_clustering, mutants_to_drop):
    if not mutants_to_drop:
        return dataset_for_clustering
    else:
        return dataset_for_clustering.loc[~dataset_for_clustering.index.isin(mutants_to_drop)]

data_for_clustering= drop_strain(data_for_clustering, ['act1-108'])
# data_for_clustering


# __split by a condition ([As])__
_05mM_data_for_clustering= data_for_clustering.loc[:, ['slope_0mM', 'slope_05mM', 'slope_pvalue_0mM', 'slope_pvalue_05mM']]
_1mM_data_for_clustering= data_for_clustering.loc[:, ['slope_0mM', 'slope_1mM', 'slope_pvalue_0mM', 'slope_pvalue_1mM']]


# #### __k-means__

# * __0.5 mM arsenic__
elbow_plot(normalize_numerical_columns(_05mM_data_for_clustering))

k_means_clustered= k_means(normalize_numerical_columns(_05mM_data_for_clustering), k= 6)
# k_means_clustered

visualize_clusters(k_means_clustered, _05mM_data_for_clustering.iloc[:, -2], _05mM_data_for_clustering.iloc[:, -1])


# * __1 mM arsenic__
elbow_plot(normalize_numerical_columns(_1mM_data_for_clustering))

k_means_clustered2= k_means(normalize_numerical_columns(_1mM_data_for_clustering), k= 6)
# k_means_clustered2

visualize_clusters(k_means_clustered2, _1mM_data_for_clustering.iloc[:, -2], _1mM_data_for_clustering.iloc[:, -1])


# #### __DBSCAN__

# * __0.5 mM arsenic__
DBSCAN_param_opt(normalize_numerical_columns(_05mM_data_for_clustering), no_of_combinations= 20, min_eps=0.1, max_eps=1, min_min_samples=3, max_min_samples=10).head(3)

dbscan_clustered= dbscan(normalize_numerical_columns(_05mM_data_for_clustering), epsilon= 0.19, no_of_points= 6)
# dbscan_clustered

visualize_clusters(dbscan_clustered, _05mM_data_for_clustering.iloc[:, -2], _05mM_data_for_clustering.iloc[:, -1])


# * __1mM arsenic__
DBSCAN_param_opt(normalize_numerical_columns(_1mM_data_for_clustering), no_of_combinations= 20, min_eps=0.1, max_eps=1, min_min_samples=3, max_min_samples=10).head(3)

dbscan_clustered2= dbscan(normalize_numerical_columns(_1mM_data_for_clustering), epsilon= 0.34, no_of_points= 3)
# dbscan_clustered2

visualize_clusters(dbscan_clustered2, _1mM_data_for_clustering.iloc[:, -2], _1mM_data_for_clustering.iloc[:, -1])
