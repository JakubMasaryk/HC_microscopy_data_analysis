#!/usr/bin/env python
# coding: utf-8

# * __libraries__

# In[4]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from scipy import stats
from scipy.stats import ttest_ind
import pwlf
from sqlalchemy import create_engine


# ---------------------------------------------------------

# * __params__

# In[7]:


plt.rcParams["legend.frameon"] = False
plt.rcParams['legend.fontsize'] = 15

plt.rcParams['axes.labelsize'] = 20
plt.rcParams['axes.labelweight'] = 'bold'

plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18  

plt.rcParams['font.size'] = 16

plt.rcParams['figure.dpi'] = 1000


# ---------------------------------------------------------------------------------

# * __inputs__

# In[10]:


#microscopy parameters
initital_timepoints_skippped= 1
microscopy_interval= 3.5
microscopy_initital_delay= 7

#mysql server connection parameters
username= 'root'
password= 'poef.qve5353'
hostname= '127.0.0.1'
port= '3306'

#paths to files
path_to_raw_file= r"C:\Users\Jakub\Desktop\figures\Figure_S1\data\Fig1_S1_data.csv"
path_to_plate_file= r"C:\Users\Jakub\Desktop\figures\Figure_S1\data\Fig1_S1_plate_layout.xlsx"


# ------------------------------------------------------------------

# * __data processiong functions__

# In[13]:


#load from a file, data processing
def raw_data_Load_and_processing_file(path, initial_delay, frequency, initial_timepoints_skipped):
    dataset= pd.read_csv(path,
                         usecols= ['WELL LABEL', 'T', 'Cells Count wv1', 'Granules Cells with Org wv2', 'Granules Org per Cell wv2', 'Granules Area wv2'],
                         converters= {'WELL LABEL':lambda x: x.replace(' - ', '0') if len(x) == 5 else x.replace(' - ', '')})
    dataset= dataset.assign(Timepoint_minutes=dataset['T'] * frequency - (frequency - initial_delay),
                            Timepoint_hours= lambda x: x['Timepoint_minutes']/60,
                            Percentage= (dataset['Granules Cells with Org wv2']/dataset['Cells Count wv1'])*100)
    dataset.columns= ['Well', 'Timepoint', 'NumberOfCells', 'NumberOfCellsContainingAggregates', 'AverageNumberOfAggregatesPerCell', 'AverageSizeOfSingleAggregates', 'TimepointMinutes', 'TimepointHours', 'PercentageOfCellsContainingAggregates']    
    dataset= dataset.reindex(columns= ['Well', 'Timepoint', 'TimepointHours', 'TimepointMinutes', 'NumberOfCells', 'NumberOfCellsContainingAggregates','PercentageOfCellsContainingAggregates', 'AverageNumberOfAggregatesPerCell', 'AverageSizeOfSingleAggregates'])    
    dataset=dataset.loc[dataset.Timepoint > initial_timepoints_skipped]
    dataset= dataset.loc[dataset.Well.isin(['N03', 'N04', 'O03', 'O04', 'P03', 'P04'])]
    return dataset

#load from db
def raw_data_Load_and_processing_db(initial_timepoints_skipped):
    
    #mysql server connection
    connection_string = f"mysql+pymysql://{username}:{password}@{hostname}:{port}/hc_microscopy_data_v2"
    engine = create_engine(connection_string) 
    
    #query to obtain the desired data
    query = "call p_wt_characterisation_data (%s, %s)"
    param1= initial_timepoints_skipped
    param2= 'basic'
    data= pd.read_sql(query, engine, params= (param1,param2,))
    
    #unifying the column names (with 'file load')
    data.columns= ['Well', 'Timepoint', 'TimepointHours', 'TimepointMinutes', 'NumberOfCells', 'NumberOfCellsContainingAggregates', 'PercentageOfCellsContainingAggregates', 'AverageNumberOfAggregatesPerCell', 'AverageSizeOfSingleAggregates']
    
    #unifying NaNs (with 'file load')
    data= data.assign(AverageSizeOfSingleAggregates= np.where((data.NumberOfCellsContainingAggregates==0)&(data.PercentageOfCellsContainingAggregates==0)&(data.AverageNumberOfAggregatesPerCell==0), np.NaN, data.AverageSizeOfSingleAggregates))
    
    return data

# 'db' to load from mysql database, 'raw file' to load from a file
def data_load(source):
    if source=='db':
        data= raw_data_Load_and_processing_db(initital_timepoints_skippped)
        return data
    elif source=='raw file':
        data= raw_data_Load_and_processing_file(path_to_raw_file, microscopy_initital_delay, microscopy_interval, initital_timepoints_skippped)
        return data
    else:
        raise ValueError(f"Invalid source input: '{source}'. Expected: 'db' or 'raw file'.")

#filling in a missing values ba linear interpolation
def missing_values(data):
    if data.isna().sum().sum() > 0:
        return data.interpolate()
    else:
        return data
    
#allows filtering data to fit into defined time-range (hours)
def time_range_hours(data, start= 0, end= 8):
    data= data.loc[(data.TimepointHours>=start)&(data.TimepointHours<=end)]
    return data

#margin of error: t-distribution, CL 95%, confidence intervals: mean +/- margin of error
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
    if  column1== np.NaN or column2== np.NaN:
        return np.NaN
    else:
        t_stat, p_value = ttest_ind(column1, column2)
    return p_value

#grouping repeats for each mutant-condition-timepoint (into list), calculating mean, std and margins of error (for CL 95%) for each group
def repeats_group_mean_std_moe95(data):
    data= data.groupby(['Strain', 'Conditions', 'Timepoint', 'TimepointHours', 'TimepointMinutes'])[['PercentageOfCellsContainingAggregates', 'AverageNumberOfAggregatesPerCell', 'AverageSizeOfSingleAggregates']].agg({'PercentageOfCellsContainingAggregates':list,
                                                                                                                                                                                                                        'AverageNumberOfAggregatesPerCell':list,
                                                                                                                                                                                                                        'AverageSizeOfSingleAggregates':list}).reset_index()
    data= data.assign(PercentageOfCellsContainingAggregatesMean= data.PercentageOfCellsContainingAggregates.apply(lambda x: np.array(x).mean()),
                      PercentageOfCellsContainingAggregatesSTD= data.PercentageOfCellsContainingAggregates.apply(lambda x: np.array(x).std(ddof= 1)),
                      PercentageOfCellsContainingAggregatesMOE95= data.PercentageOfCellsContainingAggregates.apply(lambda x: t_margin_of_error_cl95(x)),
                      AverageNumberOfAggregatesPerCellMean= data.AverageNumberOfAggregatesPerCell.apply(lambda x: np.array(x).mean()),
                      AverageNumberOfAggregatesPerCellSTD= data.AverageNumberOfAggregatesPerCell.apply(lambda x: np.array(x).std(ddof= 1)),
                      AverageNumberOfAggregatesPerCellMOE95= data.AverageNumberOfAggregatesPerCell.apply(lambda x: t_margin_of_error_cl95(x)),
                      AverageSizeOfSingleAggregatesMean= data.AverageSizeOfSingleAggregates.apply(lambda x: np.array(x).mean()),
                      AverageSizeOfSingleAggregatesSTD= data.AverageSizeOfSingleAggregates.apply(lambda x: np.array(x).std(ddof= 1)),
                      AverageSizeOfSingleAggregatesMOE95= data.AverageSizeOfSingleAggregates.apply(lambda x: t_margin_of_error_cl95(x)))
    return data

#stage bins calculation (time ranges for individual stages)
#calculates correlation coefficient (cc) for all three stages for different combinations of formation-end and relocation & fusion ends
#selects a combination with the lowest deviation of cc from 1 as perfect positive correlation (for formation and clearance) and from -1 as a perfect negative correlation (for relocation &fusion)
def stage_bins(data, step= 2, minimal_formation_length= 30, maximal_formation_end= 120, minimal_rf_length= 60, maximal_rf_end= 300):
    
    #potential time ranges 
    potential_formation_end= np.arange(data.TimepointMinutes.min() + minimal_formation_length, maximal_formation_end + step, step)
    f_rf_end_combinations= [[x, i] for x in potential_formation_end for i in np.arange(x + minimal_rf_length, maximal_rf_end + step, step)]   
    
    #empt df, to append the cc results
    break_points_df= pd.DataFrame(columns= ['f_end', 'rf_end', 'f_cc', 'rf_cc', 'c_cc']) #empty df
    
    #cc calculation for all potential time ranges
    for breakpoints in f_rf_end_combinations:
        f_end= breakpoints[0]
        rf_end= breakpoints[1]

        formation_data= data.loc[data.TimepointMinutes <= f_end]
        relocation_and_fusion_data= data.loc[(data.TimepointMinutes > f_end) & ((data.TimepointMinutes <= rf_end))]
        clearance_data= data.loc[data.TimepointMinutes > rf_end]

        formation_cc= round(formation_data.AverageNumberOfAggregatesPerCellMean.corr(formation_data.AverageSizeOfSingleAggregatesMean, method= 'pearson', min_periods= 2), 2)
        relocation_and_fusion_cc= round(relocation_and_fusion_data.AverageNumberOfAggregatesPerCellMean.corr(relocation_and_fusion_data.AverageSizeOfSingleAggregatesMean, method= 'pearson', min_periods= 2), 2)    
        clearance_cc= round(clearance_data.AverageNumberOfAggregatesPerCellMean.corr(clearance_data.AverageSizeOfSingleAggregatesMean, method= 'pearson', min_periods= 2), 2)

        single_row= pd.DataFrame([[f_end, rf_end, formation_cc, relocation_and_fusion_cc, clearance_cc]], columns= ['f_end', 'rf_end', 'f_cc', 'rf_cc', 'c_cc'])
        break_points_df= pd.concat([break_points_df, single_row], axis= 0)
    

    #calculation of total deviation from perfect positive/negative cc (total_dev=∣x−1∣+∣x+1∣+∣x−1∣)
    break_points_df= break_points_df.assign(total_dev= np.abs(break_points_df.f_cc - 1) + np.abs(break_points_df.rf_cc + 1) + np.abs(break_points_df.c_cc - 1) )
    
    #selection of time ranges with the minimal deviation
    break_points_df= break_points_df.loc[break_points_df.total_dev==break_points_df.total_dev.min()]
    
    #average the timepoints (in case of multiple entries for the same min total_dev)
    break_points_df= break_points_df.loc[:, ['f_end', 'rf_end']].mean()
    
    #define stage bins: min set to 0 and max to 100 (used in pd.cut to define stages)
    stage_bins= [0, break_points_df.loc['f_end'], break_points_df.loc['rf_end'], 1000]
    
    return stage_bins


# * __visualisation functions__

# In[15]:


def Figure_S1(data, stage_bins, export=False):
    fig, ax= plt.subplots(1, 2, figsize= (19.6, 7.2), constrained_layout= True)

    #data split
    control_data= data.loc[data.Conditions=='control']
    exposed_data= data.loc[data.Conditions!='control']

    #assign stages
    exposed_data= exposed_data.assign(Stage= pd.cut(exposed_data.TimepointMinutes, 
                                                    bins= stage_bins,
                                                    labels= ['Formation', 'Relocation & Fusion', 'Clearance']))
    
    #control vsualisation
    ax[0].scatter(control_data.AverageNumberOfAggregatesPerCellMean,
                  control_data.AverageSizeOfSingleAggregatesMean,
                  edgecolor= 'black',
                  lw=.15,
                  color= '#ECECEC')
    ax[0].set_xlabel('avg. no. of agg. per cell', weight= 'bold')
    ax[0].set_ylabel('avg. size of a single agg.', weight= 'bold')
    
    #exposed + control visualisation
    ax[1].scatter(control_data.AverageNumberOfAggregatesPerCellMean,
                  control_data.AverageSizeOfSingleAggregatesMean,
                  edgecolor= 'black',
                  lw=.15,
                  color= '#ECECEC',
                  alpha= 0.4,
                  label= 'control cells')  
    # cc= round(exposed_data.loc[exposed_data.Stage=='Formation', 'AverageNumberOfAggregatesPerCellMean'].corr(exposed_data.loc[exposed_data.Stage=='Formation', 'AverageSizeOfSingleAggregatesMean'], method= 'pearson', min_periods= 2), 2)
    ax[1].scatter(exposed_data.loc[exposed_data.Stage=='Formation', 'AverageNumberOfAggregatesPerCellMean'],
                exposed_data.loc[exposed_data.Stage=='Formation', 'AverageSizeOfSingleAggregatesMean'],
                edgecolor= 'black',
                lw=.15,
                color= '#DFE9F5',
                label= f'As-exposed cells, formation')
    x= exposed_data.loc[exposed_data.Stage=='Formation', 'AverageNumberOfAggregatesPerCellMean']
    y= exposed_data.loc[exposed_data.Stage=='Formation', 'AverageSizeOfSingleAggregatesMean']
    m, b= np.polyfit(x, y, 1)
    ax[1].plot(x, m*x+b,
             color= '#00527C',
             linewidth= 2.5)

    # cc= round(exposed_data.loc[exposed_data.Stage=='Relocation & Fusion', 'AverageNumberOfAggregatesPerCellMean'].corr(exposed_data.loc[exposed_data.Stage=='Relocation & Fusion', 'AverageSizeOfSingleAggregatesMean'], method= 'pearson', min_periods= 2), 2)
    ax[1].scatter(exposed_data.loc[exposed_data.Stage=='Relocation & Fusion', 'AverageNumberOfAggregatesPerCellMean'],
                exposed_data.loc[exposed_data.Stage=='Relocation & Fusion', 'AverageSizeOfSingleAggregatesMean'],
                edgecolor= 'black',
                lw=.15,
                color= '#FEE7CC',
                label= f'As-exposed cells, relocation & fusion')
    x= exposed_data.loc[exposed_data.Stage=='Relocation & Fusion', 'AverageNumberOfAggregatesPerCellMean']
    y= exposed_data.loc[exposed_data.Stage=='Relocation & Fusion', 'AverageSizeOfSingleAggregatesMean']
    m, b= np.polyfit(x, y, 1)
    ax[1].plot(x, m*x+b,
             color= '#FF781F',
             linewidth= 2.5)

    # cc= round(exposed_data.loc[exposed_data.Stage=='Clearance', 'AverageNumberOfAggregatesPerCellMean'].corr(exposed_data.loc[exposed_data.Stage=='Clearance', 'AverageSizeOfSingleAggregatesMean'], method= 'pearson', min_periods= 2), 2)
    ax[1].scatter(exposed_data.loc[exposed_data.Stage=='Clearance', 'AverageNumberOfAggregatesPerCellMean'],
                exposed_data.loc[exposed_data.Stage=='Clearance', 'AverageSizeOfSingleAggregatesMean'],
                edgecolor= 'black',
                lw=.15,
                color= '#D6FED2',
                label= f'As-exposed cells, clearance')
    x= exposed_data.loc[exposed_data.Stage=='Clearance', 'AverageNumberOfAggregatesPerCellMean']
    y= exposed_data.loc[exposed_data.Stage=='Clearance', 'AverageSizeOfSingleAggregatesMean']
    m, b= np.polyfit(x, y, 1)
    ax[1].plot(x, m*x+b,
             color= '#0C8001',
             linewidth= 2.5)
    ax[1].set_xlabel('avg. no. of agg. per cell', weight= 'bold')
    ax[1].set_ylabel('avg. size of a single agg.', weight= 'bold')
    ax[1].legend(frameon= False)

    
    #export
    if export== True:
        plt.savefig(r"C:\Users\Jakub\Desktop\FigureS1.png", dpi= 1000)
    elif export== False:
        pass;
    else:
        raise ValueError(f"Invalid export argument: '{export}'. Expected: boolean ('True' or 'False').")


# ------------------------------------------------------------

# * __WT analysis__

# In[18]:


_20250106= data_load('raw file')


# In[19]:


_20250106_plate= pd.read_excel(path_to_plate_file)
_20250106= _20250106.merge(_20250106_plate, how= 'left', on='Well')


# In[20]:


_20250106= missing_values(_20250106)


# In[21]:


_20250106= repeats_group_mean_std_moe95(_20250106)


# In[22]:


calculated_stage_bins= stage_bins(_20250106)


# In[23]:


_20250106.head(5)


# ----------------------------------------------------------------------------------------

# * __Figure S1__

# In[26]:


Figure_S1(_20250106, calculated_stage_bins, export= False)


# In[ ]:




