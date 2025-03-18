#!/usr/bin/env python
# coding: utf-8

# * __libraries__

# In[46]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from scipy import stats
from scipy.stats import ttest_ind
import pwlf
from sqlalchemy import create_engine


# --------------------------------------------------------------------------------------------------

# * __params__

# In[49]:


plt.rcParams["legend.frameon"] = False
plt.rcParams['legend.fontsize'] = 15

plt.rcParams['axes.labelsize'] = 20
plt.rcParams['axes.labelweight'] = 'bold'

plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18  

plt.rcParams['font.size'] = 16

plt.rcParams['figure.dpi'] = 1000


# ----------------------------------------------------------------------

# * __inputs__

# In[52]:


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
path_to_raw_file= r"C:\Users\Jakub\Desktop\figures\Figure_S2\data\FigS2_data.csv"
path_to_plate_file= r"C:\Users\Jakub\Desktop\figures\Figure_S2\data\FigS2_plate_layout.xlsx"


# ------------------------------------

# * __data processing functions__

# In[55]:


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
    dataset= dataset.loc[dataset.Well.isin(['N05', 'N06', 'N07', 'N08', 'O05', 'O06', 'O07', 'O08', 'P05', 'P06', 'P07', 'P08'])]
    return dataset

#load from db
def raw_data_Load_and_processing_db(initial_timepoints_skipped):
    
    #mysql server connection
    connection_string = f"mysql+pymysql://{username}:{password}@{hostname}:{port}/hc_microscopy_data_v2"
    engine = create_engine(connection_string) 
    
    #query to obtain the desired data
    query = "call p_wt_characterisation_data (%s, %s)"
    param1= initial_timepoints_skipped
    param2= 'pretreatment'
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
    data= data.groupby(['Strain', 'PreTreatment', 'Conditions', 'Timepoint', 'TimepointHours', 'TimepointMinutes'])[['PercentageOfCellsContainingAggregates', 'AverageNumberOfAggregatesPerCell', 'AverageSizeOfSingleAggregates']].agg({'PercentageOfCellsContainingAggregates':list,
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


# * __visualisation functions__

# In[57]:


def pre_treated_wt_analysis(data, single_timepoints, export=False):
    
    #data split
    ctrl_prt_ctrl_con= data.loc[(data.PreTreatment=='control') & (data.Conditions=='control')] #pre-treatment: 60 mins YNB-C, microscopy conditions: YNB-C
    ctrl_prt_As_con= data.loc[(data.PreTreatment=='control') & (data.Conditions=='0.5 mM As')] #pre-treatment: 60 mins YNB-C, microscopy conditions: YNB-C 0.5 mM As
    As_prt_ctrl_con= data.loc[(data.PreTreatment=='0.5 mM As') & (data.Conditions=='control')] #pre-treatment: 60 mins YNB-C 0.5 mM As, microscopy conditions: YNB-C
    As_prt_As_con= data.loc[(data.PreTreatment=='0.5 mM As') & (data.Conditions=='0.5 mM As')] #pre-treatment: 60 mins YNB-C 0.5 mM As, microscopy conditions: YNB-C 0.5 mM As
    
    #plot layout
    fig= plt.figure(figsize= (19.6, 18), constrained_layout= True)
    gs= gridspec.GridSpec(5, 4, figure=fig)
    ax1= fig.add_subplot(gs[0:2, 0:2])
    ax2= fig.add_subplot(gs[0:2, 2:4])
    ax3= fig.add_subplot(gs[2:4, 0:2])
    ax4= fig.add_subplot(gs[2:4, 2:4])
    ax5= fig.add_subplot(gs[4:5, 0:1])
    ax6= fig.add_subplot(gs[4:5, 1:2])
    ax7= fig.add_subplot(gs[4:5, 2:3])
    ax8= fig.add_subplot(gs[4:5, 3:4])
    
    #subplot 1
    ax1.errorbar(ctrl_prt_ctrl_con.TimepointMinutes,
                      ctrl_prt_ctrl_con.PercentageOfCellsContainingAggregatesMean,
                      lw= 2.5,
                      label= 'control cells',
                      color= '#00527C')
    y_min_ctrl_prt_ctrl_con= ctrl_prt_ctrl_con.PercentageOfCellsContainingAggregatesMean - ctrl_prt_ctrl_con.PercentageOfCellsContainingAggregatesMOE95
    y_max_ctrl_prt_ctrl_con= ctrl_prt_ctrl_con.PercentageOfCellsContainingAggregatesMean + ctrl_prt_ctrl_con.PercentageOfCellsContainingAggregatesMOE95
    ax1.fill_between(ctrl_prt_ctrl_con.TimepointMinutes, 
                          y_min_ctrl_prt_ctrl_con, 
                          y_max_ctrl_prt_ctrl_con, 
                          color='#DFE9F5', 
                          label="margin of error area: control cells",
                          alpha= .6)
    
    ax1.errorbar(ctrl_prt_As_con.TimepointMinutes,
                      ctrl_prt_As_con.PercentageOfCellsContainingAggregatesMean,
                      lw= 2.5,
                      label= 'As-exposed cells',
                      color= '#FF781F')
    y_min_ctrl_prt_As_con= ctrl_prt_As_con.PercentageOfCellsContainingAggregatesMean - ctrl_prt_As_con.PercentageOfCellsContainingAggregatesMOE95
    y_max_ctrl_prt_As_con= ctrl_prt_As_con.PercentageOfCellsContainingAggregatesMean + ctrl_prt_As_con.PercentageOfCellsContainingAggregatesMOE95
    ax1.fill_between(ctrl_prt_As_con.TimepointMinutes, 
                          y_min_ctrl_prt_As_con, 
                          y_max_ctrl_prt_As_con, 
                          color='#FEE7CC', 
                          label="margin of error area: As-exposed cells",
                          alpha= .6)
    ax1.set_ylim(-15, 115)
    ax1.set_xlabel('time (min)', weight= 'bold')
    ax1.set_ylabel('percentage', weight= 'bold')
    ax1.set_xticks([0,60,120,180,240,300,360,420,480,540])
    ax1.legend(frameon= False,
               ncol= 2)

    
    #subplot 2
    ax2.errorbar(As_prt_ctrl_con.TimepointMinutes,
                      As_prt_ctrl_con.PercentageOfCellsContainingAggregatesMean,
                      lw= 2.5,
                      label= 'control cells',
                      color= '#00527C')
    y_min_As_prt_ctrl_con= As_prt_ctrl_con.PercentageOfCellsContainingAggregatesMean - As_prt_ctrl_con.PercentageOfCellsContainingAggregatesMOE95
    y_max_As_prt_ctrl_con= As_prt_ctrl_con.PercentageOfCellsContainingAggregatesMean + As_prt_ctrl_con.PercentageOfCellsContainingAggregatesMOE95
    ax2.fill_between(As_prt_ctrl_con.TimepointMinutes, 
                          y_min_As_prt_ctrl_con, 
                          y_max_As_prt_ctrl_con, 
                          color='#DFE9F5', 
                          label="margin of error area: control cells",
                          alpha= .6)
    
    ax2.errorbar(As_prt_As_con.TimepointMinutes,
                      As_prt_As_con.PercentageOfCellsContainingAggregatesMean,
                      lw= 2.5,
                      label= 'As-exposed cells',
                      color= '#FF781F')
    y_min_As_prt_As_con= As_prt_As_con.PercentageOfCellsContainingAggregatesMean - As_prt_As_con.PercentageOfCellsContainingAggregatesMOE95
    y_max_As_prt_As_con= As_prt_As_con.PercentageOfCellsContainingAggregatesMean + As_prt_As_con.PercentageOfCellsContainingAggregatesMOE95
    ax2.fill_between(As_prt_As_con.TimepointMinutes, 
                          y_min_As_prt_As_con, 
                          y_max_As_prt_As_con, 
                          color='#FEE7CC', 
                          label="margin of error area: As-exposed cells",
                          alpha= .6)
    ax2.set_ylim(-15, 115)
    ax2.set_xlabel('time (min)', weight= 'bold')
    ax2.set_ylabel('percentage', weight= 'bold')
    ax2.set_xticks([0,60,120,180,240,300,360,420,480,540])
    ax2.legend(frameon= False,
               ncol= 2)
    # ax2.set_title("pre-treatment: 60' YNB-C 0.5 mM As", weight= 'bold')

    
    #subplot 3
    selected_timepoints_ctrl_prt_ctrl_con= ctrl_prt_ctrl_con.loc[ctrl_prt_ctrl_con.TimepointMinutes.isin(single_timepoints), ['TimepointHours', 'TimepointMinutes', 'PercentageOfCellsContainingAggregates', 'PercentageOfCellsContainingAggregatesMean', 'PercentageOfCellsContainingAggregatesSTD', 'PercentageOfCellsContainingAggregatesMOE95']]
    selected_timepoints_ctrl_prt_ctrl_con= selected_timepoints_ctrl_prt_ctrl_con.assign(PreviousTimepoint= selected_timepoints_ctrl_prt_ctrl_con.PercentageOfCellsContainingAggregates.shift(1))
    selected_timepoints_ctrl_prt_ctrl_con= selected_timepoints_ctrl_prt_ctrl_con.assign(p_value= selected_timepoints_ctrl_prt_ctrl_con.apply(lambda x: single_t_test(x['PercentageOfCellsContainingAggregates'], x['PreviousTimepoint']), axis= 1),
                                                                                        significance= lambda x: np.where(x['p_value']<0.001, '*', ''))

    selected_timepoints_ctrl_prt_As_con= ctrl_prt_As_con.loc[ctrl_prt_As_con.TimepointMinutes.isin(single_timepoints), ['TimepointHours', 'TimepointMinutes', 'PercentageOfCellsContainingAggregates', 'PercentageOfCellsContainingAggregatesMean', 'PercentageOfCellsContainingAggregatesSTD', 'PercentageOfCellsContainingAggregatesMOE95']]
    selected_timepoints_ctrl_prt_As_con= selected_timepoints_ctrl_prt_As_con.assign(PreviousTimepoint= selected_timepoints_ctrl_prt_As_con.PercentageOfCellsContainingAggregates.shift(1))
    selected_timepoints_ctrl_prt_As_con= selected_timepoints_ctrl_prt_As_con.assign(p_value= selected_timepoints_ctrl_prt_As_con.apply(lambda x: single_t_test(x['PercentageOfCellsContainingAggregates'], x['PreviousTimepoint']), axis= 1),
                                                                                    significance= lambda x: np.where(x['p_value']<0.001, '*', ''))

    width= .4
    x= np.arange(0, len(selected_timepoints_ctrl_prt_ctrl_con))

    ax3.bar(x-width/2,
            selected_timepoints_ctrl_prt_ctrl_con.PercentageOfCellsContainingAggregatesMean,
            yerr= selected_timepoints_ctrl_prt_ctrl_con.PercentageOfCellsContainingAggregatesMOE95,
            width= width,
            color= '#DFE9F5',
            edgecolor= '#00527C',
            capsize= 2,
            error_kw= {'elinewidth':.75},
            lw= 1,
            label= 'control cells')

    ax3.bar(x+width/2,
            selected_timepoints_ctrl_prt_As_con.PercentageOfCellsContainingAggregatesMean,
            yerr= selected_timepoints_ctrl_prt_As_con.PercentageOfCellsContainingAggregatesMOE95,
            width= width,
            color= '#FEE7CC',
            edgecolor= '#FF781F',
            capsize= 2,
            error_kw= {'elinewidth':.75},
            lw= 1,
            label= 'As-exposed cells')

    x_coor_ctrl_prt_ctrl_con= x - width/2 + 0.011
    y_coor_ctrl_prt_ctrl_con= selected_timepoints_ctrl_prt_ctrl_con.PercentageOfCellsContainingAggregatesMean + selected_timepoints_ctrl_prt_ctrl_con.PercentageOfCellsContainingAggregatesMOE95 + 4
    coordinates_ctrl= [[x, y] for x, y in zip(x_coor_ctrl_prt_ctrl_con, y_coor_ctrl_prt_ctrl_con)]
    for i, c in enumerate(coordinates_ctrl):
        cx,cy = c[0], c[1]
        ax3.text(cx, cy, f'{round(selected_timepoints_ctrl_prt_ctrl_con.p_value, 4).fillna("").iloc[i]}', ha= 'center', rotation = 90)

    x_coor_ctrl_prt_As_con= x + width/2 + 0.011
    y_coor_ctrl_prt_As_con= selected_timepoints_ctrl_prt_As_con.PercentageOfCellsContainingAggregatesMean + selected_timepoints_ctrl_prt_As_con.PercentageOfCellsContainingAggregatesMOE95 + 4
    coordinates_exp= [[x, y] for x, y in zip(x_coor_ctrl_prt_As_con, y_coor_ctrl_prt_As_con)]
    for i, c in enumerate(coordinates_exp):
        cx,cy = c[0], c[1]
        ax3.text(cx, cy, f'{round(selected_timepoints_ctrl_prt_As_con.p_value, 4).fillna("").iloc[i]}', ha= 'center', rotation = 90)

    ax3.set_xticks(x)
    ax3.set_xticklabels(single_timepoints)
    ax3.axhline(0, lw=.5, ls= '--', color= 'gray')
    ax3.set_ylim(-15, 115)
    ax3.set_xlabel('timepoint (min)', weight= 'bold')
    ax3.set_ylabel('percentage', weight= 'bold')
    ax3.legend(frameon= False, ncol= 2)

    
    #subplot 4
    selected_timepoints_As_prt_ctrl_con= As_prt_ctrl_con.loc[As_prt_ctrl_con.TimepointMinutes.isin(single_timepoints), ['TimepointHours', 'TimepointMinutes', 'PercentageOfCellsContainingAggregates', 'PercentageOfCellsContainingAggregatesMean', 'PercentageOfCellsContainingAggregatesSTD', 'PercentageOfCellsContainingAggregatesMOE95']]
    selected_timepoints_As_prt_ctrl_con= selected_timepoints_As_prt_ctrl_con.assign(PreviousTimepoint= selected_timepoints_As_prt_ctrl_con.PercentageOfCellsContainingAggregates.shift(1))
    selected_timepoints_As_prt_ctrl_con= selected_timepoints_As_prt_ctrl_con.assign(p_value= selected_timepoints_As_prt_ctrl_con.apply(lambda x: single_t_test(x['PercentageOfCellsContainingAggregates'], x['PreviousTimepoint']), axis= 1),
                                                                                        significance= lambda x: np.where(x['p_value']<0.001, '*', ''))

    selected_timepoints_As_prt_As_con= As_prt_As_con.loc[As_prt_As_con.TimepointMinutes.isin(single_timepoints), ['TimepointHours', 'TimepointMinutes', 'PercentageOfCellsContainingAggregates', 'PercentageOfCellsContainingAggregatesMean', 'PercentageOfCellsContainingAggregatesSTD', 'PercentageOfCellsContainingAggregatesMOE95']]
    selected_timepoints_As_prt_As_con= selected_timepoints_As_prt_As_con.assign(PreviousTimepoint= selected_timepoints_As_prt_As_con.PercentageOfCellsContainingAggregates.shift(1))
    selected_timepoints_As_prt_As_con= selected_timepoints_As_prt_As_con.assign(p_value= selected_timepoints_As_prt_As_con.apply(lambda x: single_t_test(x['PercentageOfCellsContainingAggregates'], x['PreviousTimepoint']), axis= 1),
                                                                                significance= lambda x: np.where(x['p_value']<0.001, '*', ''))

    width= .4
    x= np.arange(0, len(selected_timepoints_As_prt_ctrl_con))

    ax4.bar(x-width/2,
            selected_timepoints_As_prt_ctrl_con.PercentageOfCellsContainingAggregatesMean,
            yerr= selected_timepoints_As_prt_ctrl_con.PercentageOfCellsContainingAggregatesMOE95,
            width= width,
            color= '#DFE9F5',
            edgecolor= '#00527C',
            capsize= 2,
            error_kw= {'elinewidth':.75},
            lw= 1,
            label= 'control cells')

    ax4.bar(x+width/2,
            selected_timepoints_As_prt_As_con.PercentageOfCellsContainingAggregatesMean,
            yerr= selected_timepoints_As_prt_As_con.PercentageOfCellsContainingAggregatesMOE95,
            width= width,
            color= '#FEE7CC',
            edgecolor= '#FF781F',
            capsize= 2,
            error_kw= {'elinewidth':.75},
            lw= 1,
            label= 'As-exposed cells')

    x_coor_As_prt_ctrl_con= x - width/2 + 0.011
    y_coor_As_prt_ctrl_con= selected_timepoints_As_prt_ctrl_con.PercentageOfCellsContainingAggregatesMean + selected_timepoints_As_prt_ctrl_con.PercentageOfCellsContainingAggregatesMOE95 + 4
    coordinates_ctrl= [[x, y] for x, y in zip(x_coor_As_prt_ctrl_con, y_coor_As_prt_ctrl_con)]
    for i, c in enumerate(coordinates_ctrl):
        cx,cy = c[0], c[1]
        ax4.text(cx, cy, f'{round(selected_timepoints_As_prt_ctrl_con.p_value, 4).fillna("").iloc[i]}', ha= 'center', rotation = 90)

    x_coor_As_prt_As_con= x + width/2 + 0.011
    y_coor_As_prt_As_con= selected_timepoints_As_prt_As_con.PercentageOfCellsContainingAggregatesMean + selected_timepoints_As_prt_As_con.PercentageOfCellsContainingAggregatesMOE95 + 4
    coordinates_exp= [[x, y] for x, y in zip(x_coor_As_prt_As_con, y_coor_As_prt_As_con)]
    for i, c in enumerate(coordinates_exp):
        cx,cy = c[0], c[1]
        ax4.text(cx, cy, f'{round(selected_timepoints_As_prt_As_con.p_value, 4).fillna("").iloc[i]}', ha= 'center', rotation = 90)

    ax4.set_xticks(x)
    ax4.set_xticklabels(single_timepoints)
    ax4.axhline(0, lw=.5, ls= '--', color= 'gray')
    ax4.set_ylim(-15, 115)
    ax4.set_xlabel('timepoint (min)', weight= 'bold')
    ax4.set_ylabel('percentage', weight= 'bold')
    ax4.legend(frameon= False, ncol= 2)

    
    #subplot 5
    ax5.scatter(ctrl_prt_ctrl_con.AverageNumberOfAggregatesPerCellMean,
                ctrl_prt_ctrl_con.AverageSizeOfSingleAggregatesMean,
                edgecolor= 'black',
                lw=.15,
                color= '#DFE9F5')
    ax5.set_ylim(0.125, 0.45)    
    ax5.set_xlim(0, 3.25)
    ax5.set_xlabel('avg. no. of agg. per cell', weight= 'bold', fontsize= 15)
    ax5.set_ylabel('avg. size of a single agg.', weight= 'bold', fontsize= 15)

    
    #subplot 6
    ax6.scatter(ctrl_prt_As_con.AverageNumberOfAggregatesPerCellMean,
                ctrl_prt_As_con.AverageSizeOfSingleAggregatesMean,
                edgecolor= 'black',
                lw=.15,
                color= '#FEE7CC')
    ax6.set_ylim(0.125, 0.45)
    ax6.set_xlim(0, 3.25)
    ax6.set_xlabel('avg. no. of agg. per cell', weight= 'bold', fontsize= 15)
    ax6.set_ylabel('avg. size of a single agg.', weight= 'bold', fontsize= 15)

    
    #subplot 7
    ax7.scatter(As_prt_ctrl_con.AverageNumberOfAggregatesPerCellMean,
                As_prt_ctrl_con.AverageSizeOfSingleAggregatesMean,
                edgecolor= 'black',
                lw=.15,
                color= '#DFE9F5')
    ax7.set_ylim(0.125, 0.45)    
    ax7.set_xlim(0, 3.25)
    ax7.set_xlabel('avg. no. of agg. per cell', weight= 'bold', fontsize= 15)
    ax7.set_ylabel('avg. size of a single agg.', weight= 'bold', fontsize= 15)

    
    #subplot 8
    ax8.scatter(As_prt_As_con.AverageNumberOfAggregatesPerCellMean,
                As_prt_As_con.AverageSizeOfSingleAggregatesMean,
                edgecolor= 'black',
                lw=.15,
                color= '#FEE7CC')
    ax8.set_ylim(0.125, 0.45)
    ax8.set_xlim(0, 3.25)
    ax8.set_xlabel('avg. no. of agg. per cell', weight= 'bold', fontsize= 15)
    ax8.set_ylabel('avg. size of a single agg.', weight= 'bold', fontsize= 15)

    
    #export
    if export== True:
        plt.savefig(r"C:\Users\Jakub\Desktop\FigureS2.png", dpi= 1000)
    elif export== False:
        pass;
    else:
        raise ValueError(f"Invalid export argument: '{export}'. Expected: boolean ('True' or 'False').")


# ----------------------------------------------------------------------------

# * __pre-treated WT analysis__

# In[60]:


_20250106= data_load('raw file')


# In[61]:


_20250106_plate= pd.read_excel(path_to_plate_file)
_20250106= _20250106.merge(_20250106_plate, how= 'left', on='Well')


# In[62]:


_20250106= missing_values(_20250106)


# In[63]:


# _20250106= time_range_hours(_20250106, 0, 9)


# In[64]:


_20250106= repeats_group_mean_std_moe95(_20250106)


# In[65]:


_20250106.head(5)


# ------------------------------------------------------------------------------------------------------------

# * __visualisation__

# In[68]:


pre_treated_wt_analysis(_20250106, [14, 42, 91, 154, 238, 308, 420, 546], export= False)


# In[ ]:




