#!/usr/bin/env python
# coding: utf-8

# __read me__

# * __database as a data source__
#     - select 'db' as argument in 'hc_data_load' (variable 'hc_data')
#     - fill in MySQL server connection parameters ('username', 'password', 'hostname', 'port')
#     - fill in pathway to bioscreen data ('path_bioscreen_data')
#     - fill in pathway for figure export ('figure_export_path') and bioscreen-supplementary table export ('supp_table_export_path')
#         - to export set the 'export' argument to 'True' (in 'follow_up_figure_basic' and 'bioscreen_supp_fig')
#     - fill in selected mutated genes ('formated_gene_list')
#     - fill in 'selected_mutants' in 'follow_up_figure_basic' function- mutants visualised as a scatterplot (WT included by default, do not include)

# * __file as a data source__
#     - select 'file' as argument in 'hc_data_load' (variable 'hc_data')
#     - fill in pathway to hc microscopy data ('path_microscopy_data')
#     - fill in pathway to bioscreen data ('path_bioscreen_data')
#     - fill in pathway for figure export ('figure_export_path') and bioscreen-supplementary table export ('supp_table_export_path')
#         - to export set the 'export' argument to 'True' (in 'follow_up_figure_basic' and 'bioscreen_supp_fig')
#     - fill in selected mutated genes ('formated_gene_list')
#     - fill in 'selected_mutants' in 'follow_up_figure_basic' function- mutants visualised as a scatterplot (WT included by default, do not include)

# __libraries__

# In[30]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.cm as cm
import seaborn as sns
from scipy import stats
from scipy.stats import ttest_ind
import pwlf
from sqlalchemy import create_engine
import matplotlib.gridspec as gridspec
import string
import sys


# __params__

# In[32]:


plt.rcParams["legend.frameon"] = False
plt.rcParams['legend.fontsize'] = 15
plt.rcParams['axes.labelsize'] = 20
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18  
plt.rcParams['font.size'] = 16
plt.rcParams['figure.dpi'] = 1000


# __inputs__

# In[12]:


#initial timepoints skipped (generally low quality data)
initital_timepoints_skipped= 0

#p value threshold
p_value_thr= 0.05

#selected timepoints
formation= 10 #60 min
rf= 25 #150 min
clearance= 68 # 408

#mysql server connection parameters
username= 'root'
password= 'poef.qve5353'
hostname= '127.0.0.1'
port= '3306' 

#list of genes of interest, replcace the "GENE1-∞" with selected genes (e.g., "ACT1", "TPM1" etc...)
#formated for mysql stored procedure (just fill in the names)
formated_gene_list= '["HRD1", "UBR1", "SLX8", "RAD6"]'
#standard list, based on formated, used to filter bioscreen data (to selected genes)
gene_list= formated_gene_list[1:-1].split(',') #split the formatted gene list
gene_list= [x.rstrip(' ').lstrip(' ').rstrip('"').lstrip('"').lower() for x in gene_list] #format the elements/genes


# In[35]:


####paths to imported files####

#pathway to microscopy dataset (needed if data loaded from raw, processed file)
path_microscopy_data= r"C:\Users\Jakub\Desktop\figures\Figure_X_ubiquitin_ligases\data\ubiquitin_ligases_microscopy_processed_data.csv"
#pathway to bioscreen dataset (always needed, bioscreen data not part of the database)
path_bioscreen_data= r"C:\Users\Jakub\Desktop\figures\Figure_X_ubiquitin_ligases\data\ub_lig_mutants_bioscreen_processed_data.xlsx"


# In[36]:


####paths for export####

#
figure_export_path= r"C:\Users\Jakub\Desktop\fig_ub_ligases.png"
#
supp_table_export_path= r"C:\Users\Jakub\Desktop\supp_tab_ub_ligases.xlsx"


# ----------------------------------------------------------------------------------------------------------------------

# __functions__

# * __data load and processing__

# In[40]:


#microscopy data load
#source 'db' to load from database or 'file' to load from a file
def hc_data_load(source, path=path_microscopy_data, skipped_tmpts= initital_timepoints_skipped, gene_list= formated_gene_list, mysql_username= username, mysql_password= password, mysql_hostname= hostname, mysql_port= port):
    if source== 'db':
        try:
            #mysql server connection
            connection_string = f"mysql+pymysql://{mysql_username}:{mysql_password}@{mysql_hostname}:{mysql_port}/hc_microscopy_data_v2"
            engine = create_engine(connection_string)
            #query to obtain the desired sbw data
            query = "call p_follow_up_mutants_sbw (%s, %s)"
            data= pd.read_sql(query, engine, params= (skipped_tmpts, gene_list,))
            return data
        except Exception as ex:
            raise RuntimeError(f'data NOT loaded from the database, error: {ex}') #stop the script if data fails to load
    elif source== 'file':
        try:
            data= pd.read_csv(path)
            data= data.loc[data.Timepoint > skipped_tmpts]
            return data
        except Exception as ex:
            raise RuntimeError(f'data file NOT loaded, error: {ex}') #stop the script if data fails to load
    else:
        raise ValueError(f"Invalid 'source' argument: '{source}'. Expected: 'db' for database or 'file' for raw file.")

#bioscreen data load (always froma file)
def bsc_data_load(path= path_bioscreen_data):
    try:
        data= pd.read_excel(path)
        data = data.loc[:, ~data.columns.str.contains('^Unnamed')] #remove the index column ('Unnamed: 0')
        data= data.loc[((data.Strain== 'WT') | (data.Strain.isin(gene_list))) &  #filter to WT and selected genes
                       (data.AsConcentration.isin([0, 0.5]))] #filter only to control and exposed to 0.5 mM As
        return data
    except Exception as ex:
        raise RuntimeError(f'data file NOT loaded, error: {ex}') #stop the script if data fails to load

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


# * __statistics__

# In[42]:


#margin of error: t-distribution, CL 95%, confidence intervals: mean +/- margin of error
def t_margin_of_error_cl95(data):
    cl=0.95
    std= np.std(data, ddof=1)
    n=len(data)
    
    std_err= std/np.sqrt(n)
    t_score = stats.t.ppf((1 + cl) / 2, df=n - 1)  
    t_margin_of_error = t_score * std_err
    return t_margin_of_error
    # return (np.array(data).mean()+t_margin_of_error, np.array(data).mean()-t_margin_of_error)

    
#grouping repeats for each mutant-condition-timepoint (into list), calculating mean, std and margins of error (for CL 95%) for each group
def repeats_group_mean_std_moe95_two_conditions(data):
    data= data.groupby(['Mutant', 'Condition', 'Timepoint', 'TimepointHours', 'TimepointMinutes'])[['PercentageOfCellsContainingAggregates', 'AverageNumberOfAggregatesPerCell', 'AverageSizeOfSingleAggregates']].agg({'PercentageOfCellsContainingAggregates':list,
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


#t-test: independent (two sample), assumes equal variance (by default), two-tailed
def single_t_test(column, ctrl_column):
    t_stat, p_value = ttest_ind(column, ctrl_column)
    return p_value


#row-by-row (data column vs. control column), for 2 conditions
#significance: n-no sig. difference, * p_value < threshold
def t_test_two_conditions(data, ctrl_label='WT', ctrl_condition_label= 'control', exposed_condition_label= 'As-exposed', p_values_thr= p_value_thr):
    control_control_data= data.loc[(data.Mutant==ctrl_label)&(data.Condition==ctrl_condition_label), ['Condition', 'TimepointHours', 'TimepointMinutes', 'PercentageOfCellsContainingAggregates', 'AverageNumberOfAggregatesPerCell', 'AverageSizeOfSingleAggregates']]
    control_control_data.columns= ['Condition','TimepointHours', 'TimepointMinutes', 'PercentageOfCellsContainingAggregatesCtrlData', 'AverageNumberOfAggregatesPerCellCtrlData', 'AverageSizeOfSingleAggregatesCtrlData']
    data= data.merge(control_control_data, how= 'left', on= ['Condition', 'TimepointHours', 'TimepointMinutes'])
    control_exposed_data= data.loc[(data.Mutant==ctrl_label)&(data.Condition==exposed_condition_label), ['Condition', 'TimepointHours', 'TimepointMinutes', 'PercentageOfCellsContainingAggregates', 'AverageNumberOfAggregatesPerCell', 'AverageSizeOfSingleAggregates']]
    control_exposed_data.columns= ['Condition','TimepointHours', 'TimepointMinutes', 'PercentageOfCellsContainingAggregatesCtrlData', 'AverageNumberOfAggregatesPerCellCtrlData', 'AverageSizeOfSingleAggregatesCtrlData']
    data= data.merge(control_exposed_data, how= 'left', on= ['Condition', 'TimepointHours', 'TimepointMinutes'])
    data=data.assign(PercentageOfCellsContainingAggregatesCtrlData= data.PercentageOfCellsContainingAggregatesCtrlData_x.combine_first(data.PercentageOfCellsContainingAggregatesCtrlData_y),
                     AverageNumberOfAggregatesPerCellCtrlData= data.AverageNumberOfAggregatesPerCellCtrlData_x.combine_first(data.AverageNumberOfAggregatesPerCellCtrlData_y),
                     AverageSizeOfSingleAggregatesCtrlData= data.AverageSizeOfSingleAggregatesCtrlData_x.combine_first(data.AverageSizeOfSingleAggregatesCtrlData_y))
    data=data.drop(['PercentageOfCellsContainingAggregatesCtrlData_x', 'PercentageOfCellsContainingAggregatesCtrlData_y', 'AverageNumberOfAggregatesPerCellCtrlData_x', 'AverageNumberOfAggregatesPerCellCtrlData_y', 'AverageSizeOfSingleAggregatesCtrlData_x', 'AverageSizeOfSingleAggregatesCtrlData_y'], axis= 1)
    data= data.assign(PercentageOfCellsContainingAggregatesPValues= data.apply(lambda x: single_t_test(x['PercentageOfCellsContainingAggregates'], x['PercentageOfCellsContainingAggregatesCtrlData']), axis= 1),
                      AverageNumberOfAggregatesPerCellPValues= data.apply(lambda x: single_t_test(x['AverageNumberOfAggregatesPerCell'], x['AverageNumberOfAggregatesPerCellCtrlData']), axis= 1),
                      AverageSizeOfSingleAggregatesPValues= data.apply(lambda x: single_t_test(x['AverageSizeOfSingleAggregates'], x['AverageSizeOfSingleAggregatesCtrlData']), axis= 1))
    data= data.assign(PercentageOfCellsContainingAggregatesSignificance= np.where(data.PercentageOfCellsContainingAggregatesPValues==1, '', np.where(data.PercentageOfCellsContainingAggregatesPValues<p_values_thr, '*', 'n')),
                      AverageNumberOfAggregatesPerCellSignificance= np.where(data.AverageNumberOfAggregatesPerCellPValues==1, '', np.where(data.AverageNumberOfAggregatesPerCellPValues<p_values_thr, '*', 'n')),
                      AverageSizeOfSingleAggregatesSignificance= np.where(data.AverageSizeOfSingleAggregatesPValues==1, '', np.where(data.AverageSizeOfSingleAggregatesPValues<p_values_thr, '*', 'n')))
    data= data.drop(columns= ['PercentageOfCellsContainingAggregates', 'AverageNumberOfAggregatesPerCell', 'AverageSizeOfSingleAggregates', 'PercentageOfCellsContainingAggregatesCtrlData', 'AverageNumberOfAggregatesPerCellCtrlData', 'AverageSizeOfSingleAggregatesCtrlData'])
    data= data.reindex(columns= ['Mutant', 'Condition', 'Timepoint', 'TimepointHours', 'TimepointMinutes', 'PercentageOfCellsContainingAggregatesMean', 'PercentageOfCellsContainingAggregatesSTD', 'PercentageOfCellsContainingAggregatesMOE95', 'PercentageOfCellsContainingAggregatesPValues', 'PercentageOfCellsContainingAggregatesSignificance',
                                 'AverageNumberOfAggregatesPerCellMean', 'AverageNumberOfAggregatesPerCellSTD', 'AverageNumberOfAggregatesPerCellMOE95', 'AverageNumberOfAggregatesPerCellPValues', 'AverageNumberOfAggregatesPerCellSignificance', 'AverageSizeOfSingleAggregatesMean', 'AverageSizeOfSingleAggregatesSTD',
                                 'AverageSizeOfSingleAggregatesMOE95', 'AverageSizeOfSingleAggregatesPValues', 'AverageSizeOfSingleAggregatesSignificance'])
    return data


#order data by mutant names (and time) putting control (WT) on top (2 conditions)
def order_mutants_two_conditions(data, control_label= 'WT'):
    control= data.loc[data.Mutant==control_label]
    mutants= data.loc[data.Mutant!=control_label]
    
    control= control.sort_values(by=['Mutant', 'Condition', 'Timepoint'], ascending= [True, False, True]).reset_index(drop=True)
    mutants= mutants.sort_values(by=['Mutant', 'Condition', 'Timepoint'], ascending= [True, False, True]).reset_index(drop=True)
    
    data= pd.concat([control, mutants], ignore_index= True)
    return data


#calculates slopes for each growth curves within defined timerange, each biological repeat and ttest against wt control
#applied on the exponential phase, app. 10-20 h of growth
def exp_phase_slopes(data, as_concentration, timepoint1, timepoint2, p_value_threshold):  
    #filter by timpoints, select needed columns
    data= data.loc[(data.Hours >= timepoint1) & 
                   (data.Hours <= timepoint2) &
                   (data.AsConcentration == as_concentration),
                   ['Strain', 'Hours', 'OD600']]
    
    #deal with the list formating issue after reimport from excel: remove '[]', convert to list, remove empty elements ('')
    data=data.assign(OD600= data.OD600.str.rstrip(' ]').str.lstrip('[ ').apply(lambda x: x.split(' ')).apply(lambda x: [i for i in x if i != '']))
    
    #split the repeats to a separate columns
    data= data.assign(rep1= data.OD600.apply(lambda x: x[0]),
                      rep2= data.OD600.apply(lambda x: x[1]),
                      rep3= data.OD600.apply(lambda x: x[2]))
    
    #converting individual values from the list back to float
    data= data.astype({'rep1':'float64',
                       'rep2':'float64',
                       'rep3':'float64'})
    
    #fill potential NaNs by the average of the other two columns
    data= data.assign(rep1= data.rep1.fillna(data.loc[:, ['rep2', 'rep3']].mean(axis= 1)),
                      rep2= data.rep2.fillna(data.loc[:, ['rep1', 'rep3']].mean(axis= 1)),
                      rep3= data.rep3.fillna(data.loc[:, ['rep1', 'rep2']].mean(axis= 1)))
  
    #empty dataframe for slopes
    slopes= pd.DataFrame(columns= ['strain', 'slope_repeat1', 'slope_repeat2', 'slope_repeat3'])
    
    #calculate slope for each strain-repeat, append to the 'slope' df
    for strain in data.Strain.unique():
        slope_data = data[data.Strain== strain]
        
        x=  np.array((slope_data.Hours).astype('float32'))
        y1= np.array((slope_data.rep1).astype('float32'))
        y2= np.array((slope_data.rep2).astype('float32'))
        y3= np.array((slope_data.rep3).astype('float32'))
        
        slope1, intercept1 = np.polyfit(x, y1, 1)
        slope2, intercept2 = np.polyfit(x, y2, 1)
        slope3, intercept3 = np.polyfit(x, y3, 1)
        
        new_row= pd.DataFrame([{'strain':strain, 'slope_repeat1':slope1, 'slope_repeat2':slope2, 'slope_repeat3':slope3}])   
        slopes= pd.concat([slopes, new_row], axis= 0)
    
    ##statistics
    #concat slopes into lists (input into ttest function)
    slopes= slopes.assign(slope= slopes.loc[:, ['slope_repeat1', 'slope_repeat2', 'slope_repeat3']].mean(axis= 1),
                          slopes= slopes.loc[:, ['slope_repeat1', 'slope_repeat2', 'slope_repeat3']].values.tolist())
    
    #separate control and mutant slope-data
    ctrl= slopes.loc[slopes.strain== 'WT']  
    ctrl.columns= ['strain', 'slope_repeat1', 'slope_repeat2', 'slope_repeat3', 'slope', 'control_slopes']
    mut= slopes.loc[slopes.strain!= 'WT']
    
    #merge the control data to each mutant
    final_dataset= mut.merge(ctrl.loc[:, 'control_slopes'], how= 'left', left_index= True, right_index= True)
    
    #t test
    final_dataset= final_dataset.assign(p_value= round(final_dataset.apply(lambda x: single_t_test(x['slopes'], x['control_slopes']), axis= 1), 6))
    
    #significance
    final_dataset= final_dataset.assign(significance= np.where(final_dataset.p_value < p_value_threshold, '*', ''))
    
    #filtering down to needed columns
    final_dataset= final_dataset.loc[:, ['strain', 'slope', 'p_value', 'significance']]
    
    #concat back the control
    ctrl= ctrl.loc[:, ['strain', 'slope']]
    ctrl= ctrl.assign(p_value= 0,
                      significance= '-')
    
    final_dataset= pd.concat([ctrl, final_dataset], axis= 0)
    
    return final_dataset.reset_index(drop= True)


#compares OD wt control vs. mutant in a selected timepoint (ttest, 3 repeats)
#applied on a timepoint in the stationary phase, app. 48h
def stat_phase_OD(data, as_concentration, timepoint, p_value_threshold): 
    #filter by timpoint, select needed columns
    data= data.loc[(data.Hours == timepoint) &
                   (data.AsConcentration == as_concentration),
                   ['Strain', 'Hours', 'OD600']]
    
    #deal with the list formating issue after reimport from excel: remove '[]', convert to list, remove empty elements ('')
    data=data.assign(OD600= data.OD600.str.rstrip(' ]').str.lstrip('[ ').apply(lambda x: x.split(' ')).apply(lambda x: [i for i in x if i != '']))
    
    #converting to np.array and adjusting the dtype to float
    data=data.assign(OD600= data.OD600.apply(lambda x: np.array(x)))
    data=data.assign(OD600= data.OD600.apply(lambda x: x.astype('float64')))
    
    #filling potential with the mean of the others
    data=data.assign(OD600= data.OD600.apply(lambda x: np.where(np.isnan(np.array(x)), np.mean(np.array(x)[~np.isnan(np.array(x))]),np.array(x))))
    
    #split control and mutant data
    ctrl= data.loc[data.Strain=='WT']
    ctrl.columns= ['Strain', 'Hours', 'OD600_control']
    mut= data.loc[data.Strain!='WT']
    
    #average the OD600 from the repeats
    #assign control data to each mutant
    mut= mut.assign(OD= mut.loc[:, 'OD600'].apply(lambda x: np.array(x).mean()))
    mut= mut.merge(ctrl.loc[:, ['Hours', 'OD600_control']], how= 'left', on= 'Hours')
    
    #ttest
    mut=mut.assign(p_value= mut.apply(lambda x: single_t_test(x['OD600'], x['OD600_control']), axis= 1))
    mut= mut.assign(significance= np.where(mut.p_value < p_value_threshold, '*', ''))
    
    #filter down to needed columns
    mut= mut.loc[:, ['Strain', 'OD', 'p_value', 'significance']]
    
    #formatting control data
    ctrl= ctrl.assign(OD= ctrl.OD600_control.apply(lambda x: np.array(x).mean()),
                      p_value= 0,
                      significance= '-')
    ctrl= ctrl.loc[:, ['Strain', 'OD', 'p_value', 'significance']]
    
    #concat control and mutant data- final dataset
    final_dataset= pd.concat([ctrl, mut], axis= 0)
    
    #unify column labels (with slopes)
    final_dataset.columns= ['strain', 'OD', 'p_value', 'significance']
    
    return final_dataset.reset_index(drop= True)


#stage bins calculation (time ranges for individual stages)
#calculates correlation coefficient (cc) for all three stages for different combinations of formation-end and relocation & fusion ends
#selects a combination with the lowest deviation of cc from 1 as perfect positive correlation (for formation and clearance) and from -1 as a perfect negative correlation (for relocation &fusion)
def stage_bins(data, step= 2, minimal_formation_length= 30, maximal_formation_end= 120, minimal_rf_length= 60, maximal_rf_end= 300):
    
    #only exposed, control data
    data= data.loc[(data.Condition=='As-exposed') & (data.Mutant== 'WT')]
    
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

def stages(data, stage_bins):
    data= data.assign(Stage= pd.cut(data.TimepointMinutes,
                                    bins= stage_bins,
                                    labels= ['Formation', 'Relocation & Fusion', 'Clearance']))
    return data


# * __visualisation__

# In[44]:


def follow_up_figure_basic(hc_data, bsc_data, selected_mutants, subfigures_label=False, start_label= 'A', timepoint_formation= formation, timepoint_rf= rf, timepoint_cl= clearance, error_bars= 'moe', control_label= 'control', exposed_label= 'As-exposed', export= False, export_path= figure_export_path):
    fig= plt.figure(figsize= (22, 15), constrained_layout= True)
    gs= gridspec.GridSpec(9, 8, figure=fig)
    colors= sns.color_palette("tab10", n_colors= len(bsc_data.Strain.unique()))
    ctrl_patch = mpatches.Patch(color='#DFE9F5', label=control_label)
    exposed_patch = mpatches.Patch(color='#FEE7CC', label=exposed_label)
    
    ax1= fig.add_subplot(gs[0:2, 0:4])
    ax2= fig.add_subplot(gs[2:4, 0:4])    
    ax3= fig.add_subplot(gs[4:6, 0:4])    
    
    ax4= fig.add_subplot(gs[0:3, 4:6])
    ax5= fig.add_subplot(gs[0:3, 6:8])
    ax6= fig.add_subplot(gs[3:6, 4:6])
    ax7= fig.add_subplot(gs[3:6, 6:8])
    
    ax8= fig.add_subplot(gs[6:9, 0:4])
    ax9= fig.add_subplot(gs[6:9, 4:8]) 
    
    #ax1: percentage of cells containing aggregates - formation
    x= list()
    width= 0.40
    
    #significance symbol coordinates: control
    stat_significant_diff_strains_ctrl= list(hc_data.loc[(hc_data.Timepoint==timepoint_formation)&(hc_data.Condition==control_label)&(hc_data.PercentageOfCellsContainingAggregatesSignificance=='*'), 'Mutant'])
    asterisk_indexes_ctrl= [i for i, x in enumerate(list(hc_data.Mutant.unique())) if x in stat_significant_diff_strains_ctrl] #indexes for strain with significant dofference
    asterisk_x_positions_ctrl= [i+1-width/2 for i in asterisk_indexes_ctrl] #positions (indexes +1)
    
    if error_bars== 'std':
        asterisk_y_positions_ctrl= list(hc_data.loc[(hc_data.Timepoint==timepoint_formation)&(hc_data.Condition==control_label)&(hc_data.PercentageOfCellsContainingAggregatesSignificance=='*'), 'PercentageOfCellsContainingAggregatesMean'] + 
                                        hc_data.loc[(hc_data.Timepoint==timepoint_formation)&(hc_data.Condition==control_label)&(hc_data.PercentageOfCellsContainingAggregatesSignificance=='*'), 'PercentageOfCellsContainingAggregatesSTD'] + 
                                        0.75)
    elif error_bars== 'moe':
        asterisk_y_positions_ctrl= list(hc_data.loc[(hc_data.Timepoint==timepoint_formation)&(hc_data.Condition==control_label)&(hc_data.PercentageOfCellsContainingAggregatesSignificance=='*'), 'PercentageOfCellsContainingAggregatesMean'] + 
                                        hc_data.loc[(hc_data.Timepoint==timepoint_formation)&(hc_data.Condition==control_label)&(hc_data.PercentageOfCellsContainingAggregatesSignificance=='*'), 'PercentageOfCellsContainingAggregatesMOE95'] + 
                                        0.75) 
    else:
        raise ValueError(f"Invalid error_bars argument: '{error_bars}'. Expected: 'std' or 'moe'.")
        
    coordinates_ctrl= [[x, y] for x, y in zip(asterisk_x_positions_ctrl, asterisk_y_positions_ctrl)]
        
    #significance symbol coordinates: As-exposed
    stat_significant_diff_strains_exposed= list(hc_data.loc[(hc_data.Timepoint==timepoint_formation)&(hc_data.Condition==exposed_label)&(hc_data.PercentageOfCellsContainingAggregatesSignificance=='*'), 'Mutant'])
    asterisk_indexes_exposed= [i for i, x in enumerate(list(hc_data.Mutant.unique())) if x in stat_significant_diff_strains_exposed] #indexes for strain with significant dofference
    asterisk_x_positions_exposed= [i+1+width/2 for i in asterisk_indexes_exposed] #positions (indexes +1)
    
    if error_bars== 'std':
        asterisk_y_positions_exposed= list(hc_data.loc[(hc_data.Timepoint==timepoint_formation)&(hc_data.Condition==exposed_label)&(hc_data.PercentageOfCellsContainingAggregatesSignificance=='*'), 'PercentageOfCellsContainingAggregatesMean'] + 
                                           hc_data.loc[(hc_data.Timepoint==timepoint_formation)&(hc_data.Condition==exposed_label)&(hc_data.PercentageOfCellsContainingAggregatesSignificance=='*'), 'PercentageOfCellsContainingAggregatesSTD'] + 
                                           0.75)
    elif error_bars== 'moe':
        asterisk_y_positions_exposed= list(hc_data.loc[(hc_data.Timepoint==timepoint_formation)&(hc_data.Condition==exposed_label)&(hc_data.PercentageOfCellsContainingAggregatesSignificance=='*'), 'PercentageOfCellsContainingAggregatesMean'] + 
                                           hc_data.loc[(hc_data.Timepoint==timepoint_formation)&(hc_data.Condition==exposed_label)&(hc_data.PercentageOfCellsContainingAggregatesSignificance=='*'), 'PercentageOfCellsContainingAggregatesMOE95'] + 
                                           0.75) 
    else:
        raise ValueError(f"Invalid error_bars argument: '{error_bars}'. Expected: 'std' or 'moe'.")
        
    coordinates_exposed= [[x, y] for x, y in zip(asterisk_x_positions_exposed, asterisk_y_positions_exposed)]
    
    for i, strain in enumerate(hc_data.Mutant.unique()):
        a= i+1
        x.append(a)
        
        control_value= hc_data.loc[(hc_data.Mutant==strain)&(hc_data.Timepoint==timepoint_formation)&(hc_data.Condition==control_label), 'PercentageOfCellsContainingAggregatesMean']
        exposed_value= hc_data.loc[(hc_data.Mutant==strain)&(hc_data.Timepoint==timepoint_formation)&(hc_data.Condition==exposed_label), 'PercentageOfCellsContainingAggregatesMean']
        
        if error_bars== 'std':
            error_bar_ctrl= hc_data.loc[(hc_data.Mutant==strain)&(hc_data.Timepoint==timepoint_formation)&(hc_data.Condition==control_label), 'PercentageOfCellsContainingAggregatesSTD']
            error_bar_exposed= hc_data.loc[(hc_data.Mutant==strain)&(hc_data.Timepoint==timepoint_formation)&(hc_data.Condition==exposed_label), 'PercentageOfCellsContainingAggregatesSTD']
        elif error_bars== 'moe':
            error_bar_ctrl= hc_data.loc[(hc_data.Mutant==strain)&(hc_data.Timepoint==timepoint_formation)&(hc_data.Condition==control_label), 'PercentageOfCellsContainingAggregatesMOE95']
            error_bar_exposed= hc_data.loc[(hc_data.Mutant==strain)&(hc_data.Timepoint==timepoint_formation)&(hc_data.Condition==exposed_label), 'PercentageOfCellsContainingAggregatesMOE95']
        else:
            raise ValueError(f"Invalid error_bars argument: '{error_bars}'. Expected: 'std' or 'moe'.")
    
        ax1.bar(a-width/2,
                  control_value,
                  yerr= error_bar_ctrl,
                  color= '#DFE9F5',
                  edgecolor= '#00527C',
                  width=width,
                  linewidth= .66,
                  error_kw={'ecolor': 'black', 'elinewidth': .66, 'capsize': 2, 'capthick': .66, 'alpha': 0.8})
        ax1.bar(a+width/2,
                  exposed_value,
                  yerr= error_bar_exposed,
                  width=width,
                  color= '#FEE7CC',
                  edgecolor= '#FF781F',
                  linewidth= .66,
                  error_kw={'ecolor': 'black', 'elinewidth': .66, 'capsize': 2, 'capthick': .66, 'alpha': 0.8})
    ax1.set_xticks(x)
    ax1.set_xticklabels(hc_data.Mutant.unique())
    styles= ['normal' if strain=='WT' else 'italic' for strain in list(hc_data.Mutant.unique())]
    for label, style in zip(ax1.get_xticklabels(), styles):
        label.set_fontstyle(style)
    ax1.set_ylim(-10, 145)
    ax1.set_xlabel('', weight= 'bold')
    ax1.set_ylabel('percentage', weight= 'bold')
    ax1.legend(handles=[ctrl_patch, exposed_patch], loc= 'upper right', ncol= 2, frameon= False)
    for c in coordinates_ctrl:
        x, y= c[0], c[1]
        ax1.text(x, y, '*', weight= 'bold', ha= 'center')
    for c in coordinates_exposed:
        x, y= c[0], c[1]
        ax1.text(x, y, '*', weight= 'bold', ha= 'center')
    ax1.text(0.1, 0.875, 
             f'formation',
             transform=ax1.transAxes,
             fontweight= 'bold', fontsize= 14,
             verticalalignment='center',
             horizontalalignment='left')
    
    #ax2: percentage of cells containing aggregates - formation
    x= list()
    width= 0.40
    
    #significance symbol coordinates: control
    stat_significant_diff_strains_ctrl= list(hc_data.loc[(hc_data.Timepoint==timepoint_rf)&(hc_data.Condition==control_label)&(hc_data.PercentageOfCellsContainingAggregatesSignificance=='*'), 'Mutant'])
    asterisk_indexes_ctrl= [i for i, x in enumerate(list(hc_data.Mutant.unique())) if x in stat_significant_diff_strains_ctrl] #indexes for strain with significant dofference
    asterisk_x_positions_ctrl= [i+1-width/2 for i in asterisk_indexes_ctrl] #positions (indexes +1)
    
    if error_bars== 'std':
        asterisk_y_positions_ctrl= list(hc_data.loc[(hc_data.Timepoint==timepoint_rf)&(hc_data.Condition==control_label)&(hc_data.PercentageOfCellsContainingAggregatesSignificance=='*'), 'PercentageOfCellsContainingAggregatesMean'] + 
                                        hc_data.loc[(hc_data.Timepoint==timepoint_rf)&(hc_data.Condition==control_label)&(hc_data.PercentageOfCellsContainingAggregatesSignificance=='*'), 'PercentageOfCellsContainingAggregatesSTD'] + 
                                        0.75)
    elif error_bars== 'moe':
        asterisk_y_positions_ctrl= list(hc_data.loc[(hc_data.Timepoint==timepoint_rf)&(hc_data.Condition==control_label)&(hc_data.PercentageOfCellsContainingAggregatesSignificance=='*'), 'PercentageOfCellsContainingAggregatesMean'] + 
                                        hc_data.loc[(hc_data.Timepoint==timepoint_rf)&(hc_data.Condition==control_label)&(hc_data.PercentageOfCellsContainingAggregatesSignificance=='*'), 'PercentageOfCellsContainingAggregatesMOE95'] + 
                                        0.75) 
    else:
        raise ValueError(f"Invalid error_bars argument: '{error_bars}'. Expected: 'std' or 'moe'.")
        
    coordinates_ctrl= [[x, y] for x, y in zip(asterisk_x_positions_ctrl, asterisk_y_positions_ctrl)]
        
    #significance symbol coordinates: As-exposed
    stat_significant_diff_strains_exposed= list(hc_data.loc[(hc_data.Timepoint==timepoint_rf)&(hc_data.Condition==exposed_label)&(hc_data.PercentageOfCellsContainingAggregatesSignificance=='*'), 'Mutant'])
    asterisk_indexes_exposed= [i for i, x in enumerate(list(hc_data.Mutant.unique())) if x in stat_significant_diff_strains_exposed] #indexes for strain with significant dofference
    asterisk_x_positions_exposed= [i+1+width/2 for i in asterisk_indexes_exposed] #positions (indexes +1)
    
    if error_bars== 'std':
        asterisk_y_positions_exposed= list(hc_data.loc[(hc_data.Timepoint==timepoint_rf)&(hc_data.Condition==exposed_label)&(hc_data.PercentageOfCellsContainingAggregatesSignificance=='*'), 'PercentageOfCellsContainingAggregatesMean'] + 
                                           hc_data.loc[(hc_data.Timepoint==timepoint_rf)&(hc_data.Condition==exposed_label)&(hc_data.PercentageOfCellsContainingAggregatesSignificance=='*'), 'PercentageOfCellsContainingAggregatesSTD'] + 
                                           0.75)
    elif error_bars== 'moe':
        asterisk_y_positions_exposed= list(hc_data.loc[(hc_data.Timepoint==timepoint_rf)&(hc_data.Condition==exposed_label)&(hc_data.PercentageOfCellsContainingAggregatesSignificance=='*'), 'PercentageOfCellsContainingAggregatesMean'] + 
                                           hc_data.loc[(hc_data.Timepoint==timepoint_rf)&(hc_data.Condition==exposed_label)&(hc_data.PercentageOfCellsContainingAggregatesSignificance=='*'), 'PercentageOfCellsContainingAggregatesMOE95'] + 
                                           0.75) 
    else:
        raise ValueError(f"Invalid error_bars argument: '{error_bars}'. Expected: 'std' or 'moe'.")
        
    coordinates_exposed= [[x, y] for x, y in zip(asterisk_x_positions_exposed, asterisk_y_positions_exposed)]
    
    for i, strain in enumerate(hc_data.Mutant.unique()):
        a= i+1
        x.append(a)
        
        control_value= hc_data.loc[(hc_data.Mutant==strain)&(hc_data.Timepoint==timepoint_rf)&(hc_data.Condition==control_label), 'PercentageOfCellsContainingAggregatesMean']
        exposed_value= hc_data.loc[(hc_data.Mutant==strain)&(hc_data.Timepoint==timepoint_rf)&(hc_data.Condition==exposed_label), 'PercentageOfCellsContainingAggregatesMean']
        
        if error_bars== 'std':
            error_bar_ctrl= hc_data.loc[(hc_data.Mutant==strain)&(hc_data.Timepoint==timepoint_rf)&(hc_data.Condition==control_label), 'PercentageOfCellsContainingAggregatesSTD']
            error_bar_exposed= hc_data.loc[(hc_data.Mutant==strain)&(hc_data.Timepoint==timepoint_rf)&(hc_data.Condition==exposed_label), 'PercentageOfCellsContainingAggregatesSTD']
        elif error_bars== 'moe':
            error_bar_ctrl= hc_data.loc[(hc_data.Mutant==strain)&(hc_data.Timepoint==timepoint_rf)&(hc_data.Condition==control_label), 'PercentageOfCellsContainingAggregatesMOE95']
            error_bar_exposed= hc_data.loc[(hc_data.Mutant==strain)&(hc_data.Timepoint==timepoint_rf)&(hc_data.Condition==exposed_label), 'PercentageOfCellsContainingAggregatesMOE95']
        else:
            raise ValueError(f"Invalid error_bars argument: '{error_bars}'. Expected: 'std' or 'moe'.")
    
        ax2.bar(a-width/2,
                  control_value,
                  yerr= error_bar_ctrl,
                  color= '#DFE9F5',
                  edgecolor= '#00527C',
                  width=width,
                  linewidth= .66,
                  error_kw={'ecolor': 'black', 'elinewidth': .66, 'capsize': 2, 'capthick': .66, 'alpha': 0.8})
        ax2.bar(a+width/2,
                  exposed_value,
                  yerr= error_bar_exposed,
                  width=width,
                  color= '#FEE7CC',
                  edgecolor= '#FF781F',
                  linewidth= .66,
                  error_kw={'ecolor': 'black', 'elinewidth': .66, 'capsize': 2, 'capthick': .66, 'alpha': 0.8})
    ax2.set_xticks(x)
    ax2.set_xticklabels(hc_data.Mutant.unique())
    styles= ['normal' if strain=='WT' else 'italic' for strain in list(hc_data.Mutant.unique())]
    for label, style in zip(ax2.get_xticklabels(), styles):
        label.set_fontstyle(style)
    ax2.set_ylim(-10, 145)
    ax2.set_xlabel('', weight= 'bold')
    ax2.set_ylabel('percentage', weight= 'bold')
    ax2.legend(handles=[ctrl_patch, exposed_patch], loc= 'upper right', ncol= 2, frameon= False)
    for c in coordinates_ctrl:
        x, y= c[0], c[1]
        ax2.text(x, y, '*', weight= 'bold', ha= 'center')
    for c in coordinates_exposed:
        x, y= c[0], c[1]
        ax2.text(x, y, '*', weight= 'bold', ha= 'center')
    ax2.text(0.1, 0.875, 
             f'relocation & fusion',
             transform=ax2.transAxes,
             fontweight= 'bold', fontsize= 14,
             verticalalignment='center',
             horizontalalignment='left')
    
    #ax3: percentage of cells containing aggregates - clearance
    x= list()
    width= 0.40
    
    #significance symbol coordinates: control
    stat_significant_diff_strains_ctrl= list(hc_data.loc[(hc_data.Timepoint==timepoint_cl)&(hc_data.Condition==control_label)&(hc_data.PercentageOfCellsContainingAggregatesSignificance=='*'), 'Mutant'])
    asterisk_indexes_ctrl= [i for i, x in enumerate(list(hc_data.Mutant.unique())) if x in stat_significant_diff_strains_ctrl] #indexes for strain with significant dofference
    asterisk_x_positions_ctrl= [i+1-width/2 for i in asterisk_indexes_ctrl] #positions (indexes +1)
    
    if error_bars== 'std':
        asterisk_y_positions_ctrl= list(hc_data.loc[(hc_data.Timepoint==timepoint_cl)&(hc_data.Condition==control_label)&(hc_data.PercentageOfCellsContainingAggregatesSignificance=='*'), 'PercentageOfCellsContainingAggregatesMean'] + 
                                        hc_data.loc[(hc_data.Timepoint==timepoint_cl)&(hc_data.Condition==control_label)&(hc_data.PercentageOfCellsContainingAggregatesSignificance=='*'), 'PercentageOfCellsContainingAggregatesSTD'] + 
                                        0.75)
    elif error_bars== 'moe':
        asterisk_y_positions_ctrl= list(hc_data.loc[(hc_data.Timepoint==timepoint_cl)&(hc_data.Condition==control_label)&(hc_data.PercentageOfCellsContainingAggregatesSignificance=='*'), 'PercentageOfCellsContainingAggregatesMean'] + 
                                        hc_data.loc[(hc_data.Timepoint==timepoint_cl)&(hc_data.Condition==control_label)&(hc_data.PercentageOfCellsContainingAggregatesSignificance=='*'), 'PercentageOfCellsContainingAggregatesMOE95'] + 
                                        0.75) 
    else:
        raise ValueError(f"Invalid error_bars argument: '{error_bars}'. Expected: 'std' or 'moe'.")
        
    coordinates_ctrl= [[x, y] for x, y in zip(asterisk_x_positions_ctrl, asterisk_y_positions_ctrl)]
        
    #significance symbol coordinates: As-exposed
    stat_significant_diff_strains_exposed= list(hc_data.loc[(hc_data.Timepoint==timepoint_cl)&(hc_data.Condition==exposed_label)&(hc_data.PercentageOfCellsContainingAggregatesSignificance=='*'), 'Mutant'])
    asterisk_indexes_exposed= [i for i, x in enumerate(list(hc_data.Mutant.unique())) if x in stat_significant_diff_strains_exposed] #indexes for strain with significant dofference
    asterisk_x_positions_exposed= [i+1+width/2 for i in asterisk_indexes_exposed] #positions (indexes +1)
    
    if error_bars== 'std':
        asterisk_y_positions_exposed= list(hc_data.loc[(hc_data.Timepoint==timepoint_cl)&(hc_data.Condition==exposed_label)&(hc_data.PercentageOfCellsContainingAggregatesSignificance=='*'), 'PercentageOfCellsContainingAggregatesMean'] + 
                                           hc_data.loc[(hc_data.Timepoint==timepoint_cl)&(hc_data.Condition==exposed_label)&(hc_data.PercentageOfCellsContainingAggregatesSignificance=='*'), 'PercentageOfCellsContainingAggregatesSTD'] + 
                                           0.75)
    elif error_bars== 'moe':
        asterisk_y_positions_exposed= list(hc_data.loc[(hc_data.Timepoint==timepoint_cl)&(hc_data.Condition==exposed_label)&(hc_data.PercentageOfCellsContainingAggregatesSignificance=='*'), 'PercentageOfCellsContainingAggregatesMean'] + 
                                           hc_data.loc[(hc_data.Timepoint==timepoint_cl)&(hc_data.Condition==exposed_label)&(hc_data.PercentageOfCellsContainingAggregatesSignificance=='*'), 'PercentageOfCellsContainingAggregatesMOE95'] + 
                                           0.75) 
    else:
        raise ValueError(f"Invalid error_bars argument: '{error_bars}'. Expected: 'std' or 'moe'.")
        
    coordinates_exposed= [[x, y] for x, y in zip(asterisk_x_positions_exposed, asterisk_y_positions_exposed)]
    
    for i, strain in enumerate(hc_data.Mutant.unique()):
        a= i+1
        x.append(a)
        
        control_value= hc_data.loc[(hc_data.Mutant==strain)&(hc_data.Timepoint==timepoint_cl)&(hc_data.Condition==control_label), 'PercentageOfCellsContainingAggregatesMean']
        exposed_value= hc_data.loc[(hc_data.Mutant==strain)&(hc_data.Timepoint==timepoint_cl)&(hc_data.Condition==exposed_label), 'PercentageOfCellsContainingAggregatesMean']
        
        if error_bars== 'std':
            error_bar_ctrl= hc_data.loc[(hc_data.Mutant==strain)&(hc_data.Timepoint==timepoint_cl)&(hc_data.Condition==control_label), 'PercentageOfCellsContainingAggregatesSTD']
            error_bar_exposed= hc_data.loc[(hc_data.Mutant==strain)&(hc_data.Timepoint==timepoint_cl)&(hc_data.Condition==exposed_label), 'PercentageOfCellsContainingAggregatesSTD']
        elif error_bars== 'moe':
            error_bar_ctrl= hc_data.loc[(hc_data.Mutant==strain)&(hc_data.Timepoint==timepoint_cl)&(hc_data.Condition==control_label), 'PercentageOfCellsContainingAggregatesMOE95']
            error_bar_exposed= hc_data.loc[(hc_data.Mutant==strain)&(hc_data.Timepoint==timepoint_cl)&(hc_data.Condition==exposed_label), 'PercentageOfCellsContainingAggregatesMOE95']
        else:
            raise ValueError(f"Invalid error_bars argument: '{error_bars}'. Expected: 'std' or 'moe'.")
    
        ax3.bar(a-width/2,
                  control_value,
                  yerr= error_bar_ctrl,
                  color= '#DFE9F5',
                  edgecolor= '#00527C',
                  width=width,
                  linewidth= .66,
                  error_kw={'ecolor': 'black', 'elinewidth': .66, 'capsize': 2, 'capthick': .66, 'alpha': 0.8})
        ax3.bar(a+width/2,
                  exposed_value,
                  yerr= error_bar_exposed,
                  width=width,
                  color= '#FEE7CC',
                  edgecolor= '#FF781F',
                  linewidth= .66,
                  error_kw={'ecolor': 'black', 'elinewidth': .66, 'capsize': 2, 'capthick': .66, 'alpha': 0.8})
    ax3.set_xticks(x)
    ax3.set_xticklabels(hc_data.Mutant.unique())
    styles= ['normal' if strain=='WT' else 'italic' for strain in list(hc_data.Mutant.unique())]
    for label, style in zip(ax3.get_xticklabels(), styles):
        label.set_fontstyle(style)
    ax3.set_ylim(-10, 145)
    ax3.set_xlabel('', weight= 'bold')
    ax3.set_ylabel('percentage', weight= 'bold')
    ax3.legend(handles=[ctrl_patch, exposed_patch], loc= 'upper right', ncol= 2, frameon= False)
    for c in coordinates_ctrl:
        x, y= c[0], c[1]
        ax3.text(x, y, '*', weight= 'bold', ha= 'center')
    for c in coordinates_exposed:
        x, y= c[0], c[1]
        ax3.text(x, y, '*', weight= 'bold', ha= 'center')
    ax3.text(0.1, 0.875, 
             f'clearance',
             transform=ax3.transAxes,
             fontweight= 'bold', fontsize= 14,
             verticalalignment='center',
             horizontalalignment='left')
    
    #ax4: size-number covariance: WT
    wt_data= hc_data.loc[(hc_data.Mutant=='WT') & (hc_data.Condition==exposed_label)]
    cc= round(wt_data.loc[wt_data.Stage=='Formation', 'AverageNumberOfAggregatesPerCellMean'].corr(wt_data.loc[wt_data.Stage=='Formation', 'AverageSizeOfSingleAggregatesMean'], method= 'pearson', min_periods= 2), 2)
    ax4.scatter(wt_data.loc[wt_data.Stage=='Formation', 'AverageNumberOfAggregatesPerCellMean'],
                wt_data.loc[wt_data.Stage=='Formation', 'AverageSizeOfSingleAggregatesMean'],
                edgecolor= 'black',
                lw=.15,
                color= '#DFE9F5',
                label= f'FO, {cc}')
    x= wt_data.loc[wt_data.Stage=='Formation', 'AverageNumberOfAggregatesPerCellMean']
    y= wt_data.loc[wt_data.Stage=='Formation', 'AverageSizeOfSingleAggregatesMean']
    m, b= np.polyfit(x, y, 1)
    ax4.plot(x, 
             m*x+b,
             color= '#00527C',
             linewidth= 2.5)
    
    cc= round(wt_data.loc[wt_data.Stage=='Relocation & Fusion', 'AverageNumberOfAggregatesPerCellMean'].corr(wt_data.loc[wt_data.Stage=='Relocation & Fusion', 'AverageSizeOfSingleAggregatesMean'], method= 'pearson', min_periods= 2), 2)
    ax4.scatter(wt_data.loc[wt_data.Stage=='Relocation & Fusion', 'AverageNumberOfAggregatesPerCellMean'],
                wt_data.loc[wt_data.Stage=='Relocation & Fusion', 'AverageSizeOfSingleAggregatesMean'],
                edgecolor= 'black',
                lw=.15,
                color= '#FEE7CC',
                label= f'RF, {cc}')
    x= wt_data.loc[wt_data.Stage=='Relocation & Fusion', 'AverageNumberOfAggregatesPerCellMean']
    y= wt_data.loc[wt_data.Stage=='Relocation & Fusion', 'AverageSizeOfSingleAggregatesMean']
    m, b= np.polyfit(x, y, 1)
    ax4.plot(x, 
             m*x+b,
             color= '#FF781F',
             linewidth= 2.5)
    
    cc= round(wt_data.loc[wt_data.Stage=='Clearance', 'AverageNumberOfAggregatesPerCellMean'].corr(wt_data.loc[wt_data.Stage=='Clearance', 'AverageSizeOfSingleAggregatesMean'], method= 'pearson', min_periods= 2), 2)
    ax4.scatter(wt_data.loc[wt_data.Stage=='Clearance', 'AverageNumberOfAggregatesPerCellMean'],
                wt_data.loc[wt_data.Stage=='Clearance', 'AverageSizeOfSingleAggregatesMean'],
                edgecolor= 'black',
                lw=.15,
                color= '#D6FED2',
                label= f'CL, {cc}')
    x= wt_data.loc[wt_data.Stage=='Clearance', 'AverageNumberOfAggregatesPerCellMean']
    y= wt_data.loc[wt_data.Stage=='Clearance', 'AverageSizeOfSingleAggregatesMean']
    m, b= np.polyfit(x, y, 1)
    ax4.plot(x, 
             m*x+b,
             color= '#0C8001',
             linewidth= 2.5)
    wt_data_ctrl= hc_data.loc[(hc_data.Mutant=='WT') & (hc_data.Condition==control_label)]
    ax4.scatter(wt_data_ctrl.loc[:, 'AverageNumberOfAggregatesPerCellMean'],
                wt_data_ctrl.loc[:, 'AverageSizeOfSingleAggregatesMean'],
                edgecolor= 'black',
                lw=.15,
                color= '#ECECEC',
                label= f'control')
    ax4.set_xlabel('agg. no.', weight= 'bold')
    ax4.set_ylabel('agg. size', weight= 'bold')
    ax4.legend(frameon= True, fontsize= 13, loc= 'upper right')
    ax4.set_xlim(0, 3.2)
    ax4.set_xticks([0, 0.8, 1.6, 2.4, 3.2])
    ax4.set_ylim(0.12, 0.5)
    ax4.set_yticks([0.1, 0.2, 0.3, 0.4, 0.5])
    ax4.text(00.85, 0.1, 
         f'WT',
         transform=ax4.transAxes,
         fontweight= 'bold', fontsize= 20,
         verticalalignment='center',
         horizontalalignment='center')
    
    #ax5: size-number covariance: selected mutant no. 1
    mutant1_data= hc_data.loc[(hc_data.Mutant==selected_mutants[0]) & (hc_data.Condition==exposed_label)]
    cc= round(mutant1_data.loc[mutant1_data.Stage=='Formation', 'AverageNumberOfAggregatesPerCellMean'].corr(mutant1_data.loc[mutant1_data.Stage=='Formation', 'AverageSizeOfSingleAggregatesMean'], method= 'pearson', min_periods= 2), 2)
    ax5.scatter(mutant1_data.loc[mutant1_data.Stage=='Formation', 'AverageNumberOfAggregatesPerCellMean'],
                mutant1_data.loc[mutant1_data.Stage=='Formation', 'AverageSizeOfSingleAggregatesMean'],
                edgecolor= 'black',
                lw=.15,
                color= '#DFE9F5',
                label= f'FO, {cc}')
    x= mutant1_data.loc[mutant1_data.Stage=='Formation', 'AverageNumberOfAggregatesPerCellMean']
    y= mutant1_data.loc[mutant1_data.Stage=='Formation', 'AverageSizeOfSingleAggregatesMean']
    m, b= np.polyfit(x, y, 1)
    ax5.plot(x, 
             m*x+b,
             color= '#00527C',
             linewidth= 2.5)
    
    cc= round(mutant1_data.loc[mutant1_data.Stage=='Relocation & Fusion', 'AverageNumberOfAggregatesPerCellMean'].corr(mutant1_data.loc[mutant1_data.Stage=='Relocation & Fusion', 'AverageSizeOfSingleAggregatesMean'], method= 'pearson', min_periods= 2), 2)
    ax5.scatter(mutant1_data.loc[mutant1_data.Stage=='Relocation & Fusion', 'AverageNumberOfAggregatesPerCellMean'],
                mutant1_data.loc[mutant1_data.Stage=='Relocation & Fusion', 'AverageSizeOfSingleAggregatesMean'],
                edgecolor= 'black',
                lw=.15,
                color= '#FEE7CC',
                label= f'RF, {cc}')
    x= mutant1_data.loc[mutant1_data.Stage=='Relocation & Fusion', 'AverageNumberOfAggregatesPerCellMean']
    y= mutant1_data.loc[mutant1_data.Stage=='Relocation & Fusion', 'AverageSizeOfSingleAggregatesMean']
    m, b= np.polyfit(x, y, 1)
    ax5.plot(x, 
             m*x+b,
             color= '#FF781F',
             linewidth= 2.5)
    
    cc= round(mutant1_data.loc[mutant1_data.Stage=='Clearance', 'AverageNumberOfAggregatesPerCellMean'].corr(mutant1_data.loc[mutant1_data.Stage=='Clearance', 'AverageSizeOfSingleAggregatesMean'], method= 'pearson', min_periods= 2), 2)
    ax5.scatter(mutant1_data.loc[mutant1_data.Stage=='Clearance', 'AverageNumberOfAggregatesPerCellMean'],
                mutant1_data.loc[mutant1_data.Stage=='Clearance', 'AverageSizeOfSingleAggregatesMean'],
                edgecolor= 'black',
                lw=.15,
                color= '#D6FED2',
                label= f'CL, {cc}')
    x= mutant1_data.loc[mutant1_data.Stage=='Clearance', 'AverageNumberOfAggregatesPerCellMean']
    y= mutant1_data.loc[mutant1_data.Stage=='Clearance', 'AverageSizeOfSingleAggregatesMean']
    m, b= np.polyfit(x, y, 1)
    ax5.plot(x, 
             m*x+b,
             color= '#0C8001',
             linewidth= 2.5)
    mutant1_data_ctrl= hc_data.loc[(hc_data.Mutant==selected_mutants[0]) & (hc_data.Condition==control_label)]
    ax5.scatter(mutant1_data_ctrl.loc[:, 'AverageNumberOfAggregatesPerCellMean'],
                mutant1_data_ctrl.loc[:, 'AverageSizeOfSingleAggregatesMean'],
                edgecolor= 'black',
                lw=.15,
                color= '#ECECEC',
                label= f'control')
    ax5.set_xlabel('agg. no.', weight= 'bold')
    ax5.set_ylabel('agg. size', weight= 'bold')
    ax5.legend(frameon= True, fontsize= 13, loc= 'upper right')
    ax5.set_xlim(0, 3.2)
    ax5.set_xticks([0, 0.8, 1.6, 2.4, 3.2])
    ax5.set_ylim(0.12, 0.5)
    ax5.set_yticks([0.1, 0.2, 0.3, 0.4, 0.5])
    ax5.text(0.85, 0.1, 
         f'{selected_mutants[0]}',
         transform=ax5.transAxes,
         fontweight='bold', style='italic', fontsize= 20,
         verticalalignment='center',
         horizontalalignment='center')
    
    #ax6: size-number covariance: selected mutant no. 2
    mutant1_data= hc_data.loc[(hc_data.Mutant==selected_mutants[1]) & (hc_data.Condition==exposed_label)]
    cc= round(mutant1_data.loc[mutant1_data.Stage=='Formation', 'AverageNumberOfAggregatesPerCellMean'].corr(mutant1_data.loc[mutant1_data.Stage=='Formation', 'AverageSizeOfSingleAggregatesMean'], method= 'pearson', min_periods= 2), 2)
    ax6.scatter(mutant1_data.loc[mutant1_data.Stage=='Formation', 'AverageNumberOfAggregatesPerCellMean'],
                mutant1_data.loc[mutant1_data.Stage=='Formation', 'AverageSizeOfSingleAggregatesMean'],
                edgecolor= 'black',
                lw=.15,
                color= '#DFE9F5',
                label= f'FO, {cc}')
    x= mutant1_data.loc[mutant1_data.Stage=='Formation', 'AverageNumberOfAggregatesPerCellMean']
    y= mutant1_data.loc[mutant1_data.Stage=='Formation', 'AverageSizeOfSingleAggregatesMean']
    m, b= np.polyfit(x, y, 1)
    ax6.plot(x, 
             m*x+b,
             color= '#00527C',
             linewidth= 2.5)
    
    cc= round(mutant1_data.loc[mutant1_data.Stage=='Relocation & Fusion', 'AverageNumberOfAggregatesPerCellMean'].corr(mutant1_data.loc[mutant1_data.Stage=='Relocation & Fusion', 'AverageSizeOfSingleAggregatesMean'], method= 'pearson', min_periods= 2), 2)
    ax6.scatter(mutant1_data.loc[mutant1_data.Stage=='Relocation & Fusion', 'AverageNumberOfAggregatesPerCellMean'],
                mutant1_data.loc[mutant1_data.Stage=='Relocation & Fusion', 'AverageSizeOfSingleAggregatesMean'],
                edgecolor= 'black',
                lw=.15,
                color= '#FEE7CC',
                label= f'RF, {cc}')
    x= mutant1_data.loc[mutant1_data.Stage=='Relocation & Fusion', 'AverageNumberOfAggregatesPerCellMean']
    y= mutant1_data.loc[mutant1_data.Stage=='Relocation & Fusion', 'AverageSizeOfSingleAggregatesMean']
    m, b= np.polyfit(x, y, 1)
    ax6.plot(x, 
             m*x+b,
             color= '#FF781F',
             linewidth= 2.5)
    
    cc= round(mutant1_data.loc[mutant1_data.Stage=='Clearance', 'AverageNumberOfAggregatesPerCellMean'].corr(mutant1_data.loc[mutant1_data.Stage=='Clearance', 'AverageSizeOfSingleAggregatesMean'], method= 'pearson', min_periods= 2), 2)
    ax6.scatter(mutant1_data.loc[mutant1_data.Stage=='Clearance', 'AverageNumberOfAggregatesPerCellMean'],
                mutant1_data.loc[mutant1_data.Stage=='Clearance', 'AverageSizeOfSingleAggregatesMean'],
                edgecolor= 'black',
                lw=.15,
                color= '#D6FED2',
                label= f'CL, {cc}')
    x= mutant1_data.loc[mutant1_data.Stage=='Clearance', 'AverageNumberOfAggregatesPerCellMean']
    y= mutant1_data.loc[mutant1_data.Stage=='Clearance', 'AverageSizeOfSingleAggregatesMean']
    m, b= np.polyfit(x, y, 1)
    ax6.plot(x, 
             m*x+b,
             color= '#0C8001',
             linewidth= 2.5)
    mutant1_data_ctrl= hc_data.loc[(hc_data.Mutant==selected_mutants[1]) & (hc_data.Condition==control_label)]
    ax6.scatter(mutant1_data_ctrl.loc[:, 'AverageNumberOfAggregatesPerCellMean'],
                mutant1_data_ctrl.loc[:, 'AverageSizeOfSingleAggregatesMean'],
                edgecolor= 'black',
                lw=.15,
                color= '#ECECEC',
                label= f'control')
    ax6.set_xlabel('agg. no.', weight= 'bold')
    ax6.set_ylabel('agg. size', weight= 'bold')
    ax6.legend(frameon= True, fontsize= 13, loc= 'upper right')
    ax6.set_xlim(0, 3.2)
    ax6.set_xticks([0, 0.8, 1.6, 2.4, 3.2])
    ax6.set_ylim(0.12, 0.5)
    ax6.set_yticks([0.1, 0.2, 0.3, 0.4, 0.5])
    ax6.text(0.85, 0.1, 
         f'{selected_mutants[1]}',
         transform=ax6.transAxes,
         fontweight='bold', style='italic', fontsize= 20,
         verticalalignment='center',
         horizontalalignment='center')
    
    #ax7: size-number covariance: selected mutant no. 3
    if len(selected_mutants) < 3:
        ax7.axis('off')
    elif len(selected_mutants)==3:
        mutant1_data= hc_data.loc[(hc_data.Mutant==selected_mutants[2]) & (hc_data.Condition==exposed_label)]
        cc= round(mutant1_data.loc[mutant1_data.Stage=='Formation', 'AverageNumberOfAggregatesPerCellMean'].corr(mutant1_data.loc[mutant1_data.Stage=='Formation', 'AverageSizeOfSingleAggregatesMean'], method= 'pearson', min_periods= 2), 2)
        ax7.scatter(mutant1_data.loc[mutant1_data.Stage=='Formation', 'AverageNumberOfAggregatesPerCellMean'],
                    mutant1_data.loc[mutant1_data.Stage=='Formation', 'AverageSizeOfSingleAggregatesMean'],
                    edgecolor= 'black',
                    lw=.15,
                    color= '#DFE9F5',
                    label= f'FO, {cc}')
        x= mutant1_data.loc[mutant1_data.Stage=='Formation', 'AverageNumberOfAggregatesPerCellMean']
        y= mutant1_data.loc[mutant1_data.Stage=='Formation', 'AverageSizeOfSingleAggregatesMean']
        m, b= np.polyfit(x, y, 1)
        ax7.plot(x, 
                 m*x+b,
                 color= '#00527C',
                 linewidth= 2.5)

        cc= round(mutant1_data.loc[mutant1_data.Stage=='Relocation & Fusion', 'AverageNumberOfAggregatesPerCellMean'].corr(mutant1_data.loc[mutant1_data.Stage=='Relocation & Fusion', 'AverageSizeOfSingleAggregatesMean'], method= 'pearson', min_periods= 2), 2)
        ax7.scatter(mutant1_data.loc[mutant1_data.Stage=='Relocation & Fusion', 'AverageNumberOfAggregatesPerCellMean'],
                    mutant1_data.loc[mutant1_data.Stage=='Relocation & Fusion', 'AverageSizeOfSingleAggregatesMean'],
                    edgecolor= 'black',
                    lw=.15,
                    color= '#FEE7CC',
                    label= f'RF, {cc}')
        x= mutant1_data.loc[mutant1_data.Stage=='Relocation & Fusion', 'AverageNumberOfAggregatesPerCellMean']
        y= mutant1_data.loc[mutant1_data.Stage=='Relocation & Fusion', 'AverageSizeOfSingleAggregatesMean']
        m, b= np.polyfit(x, y, 1)
        ax7.plot(x, 
                 m*x+b,
                 color= '#FF781F',
                 linewidth= 2.5)

        cc= round(mutant1_data.loc[mutant1_data.Stage=='Clearance', 'AverageNumberOfAggregatesPerCellMean'].corr(mutant1_data.loc[mutant1_data.Stage=='Clearance', 'AverageSizeOfSingleAggregatesMean'], method= 'pearson', min_periods= 2), 2)
        ax7.scatter(mutant1_data.loc[mutant1_data.Stage=='Clearance', 'AverageNumberOfAggregatesPerCellMean'],
                    mutant1_data.loc[mutant1_data.Stage=='Clearance', 'AverageSizeOfSingleAggregatesMean'],
                    edgecolor= 'black',
                    lw=.15,
                    color= '#D6FED2',
                    label= f'CL, {cc}')
        x= mutant1_data.loc[mutant1_data.Stage=='Clearance', 'AverageNumberOfAggregatesPerCellMean']
        y= mutant1_data.loc[mutant1_data.Stage=='Clearance', 'AverageSizeOfSingleAggregatesMean']
        m, b= np.polyfit(x, y, 1)
        ax7.plot(x, 
                 m*x+b,
                 color= '#0C8001',
                 linewidth= 2.5)
        mutant1_data_ctrl= hc_data.loc[(hc_data.Mutant==selected_mutants[2]) & (hc_data.Condition==control_label)]
        ax7.scatter(mutant1_data_ctrl.loc[:, 'AverageNumberOfAggregatesPerCellMean'],
                    mutant1_data_ctrl.loc[:, 'AverageSizeOfSingleAggregatesMean'],
                    edgecolor= 'black',
                    lw=.15,
                    color= '#ECECEC',
                    label= f'control')
        ax7.set_xlabel('agg. no.', weight= 'bold')
        ax7.set_ylabel('agg. size', weight= 'bold')
        ax7.legend(frameon= True, fontsize= 13, loc= 'upper right')
        ax7.set_xlim(0, 3.2)
        ax7.set_xticks([0, 0.8, 1.6, 2.4, 3.2])
        ax7.set_ylim(0.12, 0.5)
        ax7.set_yticks([0.1, 0.2, 0.3, 0.4, 0.5])
        ax7.text(0.85, 0.1, 
             f'{selected_mutants[2]}',
             transform=ax7.transAxes,
             fontweight='bold', style='italic', fontsize= 20,
             verticalalignment='center',
             horizontalalignment='center')
    else:
        raise ValueError(f"Select 2 or 3 mutants")
    
    
    #ax8: growth curves non-exposed
    bsc_data_non_exposed= bsc_data.loc[bsc_data.AsConcentration==0]
    for i, strain in enumerate(bsc_data_non_exposed.Strain.unique()):
        selected_strain= bsc_data_non_exposed.loc[bsc_data_non_exposed.Strain==strain]
        ax8.plot(selected_strain.Hours,
                 selected_strain.OD600Mean,
                 lw= 4,
                 label= f'{strain}' if strain=='WT' else f'${strain}$',
                 color= colors[i])
    ax8.legend(frameon= False, loc= 'upper right', ncol= 5)
    ax8.set_ylim(0, 2.25)
    ax8.set_xlabel('time (h)')
    ax8.set_ylabel('$\mathregular{OD_{600}}$')
    
    #ax9: growth curves exposed
    bsc_data_exposed= bsc_data.loc[bsc_data.AsConcentration==0.5]
    for i, strain in enumerate(bsc_data_exposed.Strain.unique()):
        selected_strain= bsc_data_exposed.loc[bsc_data_exposed.Strain==strain]
        ax9.plot(selected_strain.Hours,
                 selected_strain.OD600Mean,
                 lw= 4,
                 label= f'{strain}' if strain=='WT' else f'${strain}$',
                 color= colors[i])
    ax9.legend(frameon= False, loc= 'upper right', ncol= 5)
    ax9.set_ylim(0, 2.25)
    ax9.set_xlabel('time (h)')
    ax9.set_ylabel('$\mathregular{OD_{600}}$')
    
    #subfigure labels
    alphabet = list(string.ascii_uppercase)
    alphabet_labels = alphabet[alphabet.index(start_label):]
    
    if subfigures_label== False:
        pass
    elif subfigures_label== True:
        ax1.text(0.03, 0.825, 
                 f'{alphabet_labels[0]}',
                 transform=ax1.transAxes,
                 fontweight='bold', fontsize= 30,
                 verticalalignment='center',
                 horizontalalignment='center')
        ax2.text(0.03, 0.825, 
                 f'{alphabet_labels[1]}',
                 transform=ax2.transAxes,
                 fontweight='bold', fontsize= 30,
                 verticalalignment='center',
                 horizontalalignment='center')
        ax3.text(0.03, 0.825, 
                 f'{alphabet_labels[2]}',
                 transform=ax3.transAxes,
                 fontweight='bold', fontsize= 30,
                 verticalalignment='center',
                 horizontalalignment='center')
        ax4.text(0.075, 0.9, 
                 f'{alphabet_labels[3]}',
                 transform=ax4.transAxes,
                 fontweight='bold', fontsize= 30,
                 verticalalignment='center',
                 horizontalalignment='center')
        ax5.text(0.075, 0.9, 
                 f'{alphabet_labels[4]}',
                 transform=ax5.transAxes,
                 fontweight='bold', fontsize= 30,
                 verticalalignment='center',
                 horizontalalignment='center')
        ax6.text(0.075, 0.9, 
                 f'{alphabet_labels[5]}',
                 transform=ax6.transAxes,
                 fontweight='bold', fontsize= 30,
                 verticalalignment='center',
                 horizontalalignment='center')
        ax7.text(0.075, 0.9, 
                 f'{alphabet_labels[6]}',
                 transform=ax7.transAxes,
                 fontweight='bold', fontsize= 30,
                 verticalalignment='center',
                 horizontalalignment='center')
        ax8.text(0.032, 0.9, 
                 f'{alphabet_labels[7]}',
                 transform=ax8.transAxes,
                 fontweight='bold', fontsize= 30,
                 verticalalignment='center',
                 horizontalalignment='center')
        ax9.text(0.032, 0.9, 
                 f'{alphabet_labels[8]}',
                 transform=ax9.transAxes,
                 fontweight='bold', fontsize= 30,
                 verticalalignment='center',
                 horizontalalignment='center')
    else:
        raise ValueError(f"Invalid 'label' argument: '{label}'. Expected: 'True' or 'False'.")
        
    #export
    if export== True:
        plt.savefig(export_path, dpi= 1000)
    elif export== False:
        pass;
    else:
        raise ValueError(f"Invalid export argument: '{export}'. Expected: boolean ('True' or 'False').")
        

#supplementary table to bioscreen data (growth curves part of the main figure)
#based on functions: 'exp_phase_slopes' and 'stat_phase_OD'
#comparison of slopes during exponential phase (default betweeen 10-20h) and OD of a single timepoint in stationary phase (default at 48h)
#analysis of control/non-exposed cells vs. As-exposed cells at selected concentration (default 0.5 mM) 
#all argument (except data) set to default
#to export set 'export' to True and predefine the pathway (see 'inputs' section)
def bioscreen_supp_fig(data, as_concentration= 0.5, slope_timepoint1= 10, slope_timepoint2= 20, stat_phase_timepoint= 48, p_value_threshold= p_value_thr, export= False, path= supp_table_export_path):
    
    ###slopes during the exponential phase
    #slope control conditions (non-exposed cells)
    slope_ctrl= exp_phase_slopes(data= data, as_concentration= 0, timepoint1= slope_timepoint1, timepoint2= slope_timepoint2, p_value_threshold= p_value_threshold)
    slope_ctrl= slope_ctrl.assign(as_concentration= '0 mM')   
    #slope exposed cells (selected As concentration)
    slope_exposed= exp_phase_slopes(data= data, as_concentration= as_concentration, timepoint1= slope_timepoint1, timepoint2= slope_timepoint2, p_value_threshold= p_value_threshold)
    slope_exposed= slope_exposed.assign(as_concentration= f'{as_concentration} mM')  
    #concat
    slopes= pd.concat([slope_ctrl, slope_exposed], axis= 0)
    slopes= slopes.reindex(columns= ['strain', 'as_concentration', 'slope', 'p_value', 'significance'])
    
    
    #single timepoint in the stationary phase
    #OD control conditions (non-exposed cells)
    od_ctrl= stat_phase_OD(data= data, as_concentration= 0, timepoint= stat_phase_timepoint, p_value_threshold= p_value_threshold)
    od_ctrl= od_ctrl.assign(as_concentration= '0 mM')
    #OD exposed cells (selected As concentration)
    od_exposed= stat_phase_OD(data= data, as_concentration= as_concentration, timepoint= stat_phase_timepoint, p_value_threshold= p_value_threshold)
    od_exposed= od_exposed.assign(as_concentration= f'{as_concentration} mM')
    #concat
    od= pd.concat([od_ctrl, od_exposed], axis= 0)
    od= od.reindex(columns= ['strain', 'as_concentration', 'OD', 'p_value', 'significance'])
    
    
    #export (slopes to sheet 1, stationary-phae OD to sheet 2)
    if export== True:
        try:
            # exporting to xlsx, two separate sheets
            with pd.ExcelWriter(path) as writer:
                slopes.to_excel(writer, sheet_name='exp_phase_slopes', index=False)
                od.to_excel(writer, sheet_name='stat_phase_OD', index=False)
        except Exception as e:
            print(f'data NOT exported, error: {e}')
    
    elif export== False:
        return slopes, od
    
    else:
        raise ValueError(f"Invalid export argument: '{export}'. Expected: 'True' or 'False'.")


# -------------------------------------------------------------------------------------

# __data load, processing and statistical analysis__

# * __HC microscopy data__

# In[48]:


#load the data
hc_data= hc_data_load(source= 'db')


# In[49]:


# hc_data.to_csv(r"C:\Users\Jakub\Desktop\figures\Figure_X_ubiquitin_ligases\data\ubiquitin_ligases_processed_data.csv", encoding= 'utf-8', index= False)


# In[50]:


#interpolate potential missing values
hc_data= missing_values(hc_data)


# In[51]:


#define timerange (for follow-up generally limited to 0-7 hours)
hc_data= time_range_hours(hc_data, 0, 7)


# In[52]:


hc_data= repeats_group_mean_std_moe95_two_conditions(hc_data)
hc_data= t_test_two_conditions(hc_data)
hc_data= order_mutants_two_conditions(hc_data)


# In[53]:


hc_data.head()


# * __bisocreen data__

# In[55]:


bsc_data= bsc_data_load()


# In[56]:


bsc_data.head()


# * __stages__

# In[58]:


stage_bins_list= stage_bins(hc_data)
stage_bins_list


# In[59]:


hc_data= stages(hc_data, stage_bins_list)


# -----------------------------------------------------------------------------------------------------------

# __visualisation__

# * __fig. X: ubiquitin-ligases__

# In[63]:


follow_up_figure_basic(hc_data, bsc_data, selected_mutants= ['hrd1', 'ubr1', 'slx8'], subfigures_label= True, export= False)


# * __fig. SX: ubiquitin-ligases__

# In[65]:


bioscreen_supp_fig(data= bsc_data, export= False)

