# __READ ME__

# __data source__
# - __relational database__
#     - relational database and needs to be created and data uploaded prior to data load (see: https://github.com/JakubMasaryk/HC_microscopy_database)
#     - fill in the mysql server parameters ('username', 'password', 'hostname', 'port') in the section bellow
#     - use argument 'db' for 'data_load function' bellow
# - __raw file__
#     - define the pathway to the raw file ('path_to_raw_file') 
#     - define the pathway to the lookup table with descriptions ('path_to_plate_file')
# - __processed file (prefered method)__
#     - define the pathway to the processed file ('path_to_processed_file')
#     - processed file part of published supp. material

# __export__
# - __figure__
#     - exported as .PNG with 1000 dpi
#     - define the pathway for export ('path_for_export')
#     - set the 'export' argument (in the visualisation function) to 'True'

# -----------------------------------------------------------------------------------
# __libraries__
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
# __params__
plt.rcParams["legend.frameon"] = False
plt.rcParams['legend.fontsize'] = 12
plt.rcParams['axes.labelsize'] = 19
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['xtick.labelsize'] = 16
plt.rcParams['ytick.labelsize'] = 16  
plt.rcParams['font.size'] = 12
plt.rcParams['figure.dpi'] = 1000
plt.rcParams['font.family'] = 'Calibri'


# ---------------------------------------------------------------------------------
# __inputs__
###microscopy parameters
initital_timepoints_skippped= 1
microscopy_interval= 3.5
microscopy_initital_delay= 7

###mysql server connection parameters
username= ''
password= ''
hostname= ''
port= ''

###paths to files
#processed file
path_to_processed_file= r"...\Fig1_S1_processed_data.csv"

#raw file 
path_to_raw_file= r"...\Fig1_S1_raw_data.csv"
#plate layout lookup table
path_to_plate_file= r"...\Fig1_S1_plate_layout.xlsx"

###path for export
path_for_export= r"...\Figure1.png"


# ------------------------------------------------------------------
# __data processiong functions__
#load from a raw file, data processing
def raw_data_Load_and_processing_file(path, initial_delay, frequency, skipped_tmpts):
    dataset= pd.read_csv(path,
                         usecols= ['WELL LABEL', 'T', 'Cells Count wv1', 'Granules Cells with Org wv2', 'Granules Org per Cell wv2', 'Granules Area wv2'],
                         converters= {'WELL LABEL':lambda x: x.replace(' - ', '0') if len(x) == 5 else x.replace(' - ', '')})
    dataset= dataset.assign(Timepoint_minutes=dataset['T'] * frequency - (frequency - initial_delay),
                            Timepoint_hours= lambda x: x['Timepoint_minutes']/60,
                            Percentage= (dataset['Granules Cells with Org wv2']/dataset['Cells Count wv1'])*100)
    dataset.columns= ['Well', 'Timepoint', 'NumberOfCells', 'NumberOfCellsContainingAggregates', 'AverageNumberOfAggregatesPerCell', 'AverageSizeOfSingleAggregates', 'TimepointMinutes', 'TimepointHours', 'PercentageOfCellsContainingAggregates']    
    dataset= dataset.reindex(columns= ['Well', 'Timepoint', 'TimepointHours', 'TimepointMinutes', 'NumberOfCells', 'NumberOfCellsContainingAggregates','PercentageOfCellsContainingAggregates', 'AverageNumberOfAggregatesPerCell', 'AverageSizeOfSingleAggregates'])    
    dataset=dataset.loc[dataset.Timepoint > skipped_tmpts]
    dataset= dataset.loc[dataset.Well.isin(['N03', 'N04', 'O03', 'O04', 'P03', 'P04'])]
    return dataset

#load from a processed-data file
def procesed_data_Load_file(path, skipped_tmpts):
    dataset= pd.read_csv(path)    
    dataset=dataset.loc[dataset.Timepoint > skipped_tmpts]
    return dataset

#load from db
def raw_data_Load_and_processing_db(skipped_tmpts):
    
    #mysql server connection
    connection_string = f"mysql+pymysql://{username}:{password}@{hostname}:{port}/hc_microscopy_data_v2"
    engine = create_engine(connection_string) 
    
    #query to obtain the desired data
    query = "call p_wt_characterisation_data_basic (%s)"
    param1= skipped_tmpts
    data= pd.read_sql(query, engine, params= (param1,))
    
    #unifying NaNs (with 'file load')
    data= data.assign(AverageSizeOfSingleAggregates= np.where((data.NumberOfCellsContainingAggregates==0)&(data.PercentageOfCellsContainingAggregates==0)&(data.AverageNumberOfAggregatesPerCell==0), np.NaN, data.AverageSizeOfSingleAggregates))
    
    return data

# 'db' to load from mysql database, 'raw file' to load from a file, 'processed file' to load the processed file
#datasets merge to 'plate' lookup table 
#take pre-defined parameters from the section 'Inputs' as arguments
def data_load(source):
    if source=='db':
        data= raw_data_Load_and_processing_db(initital_timepoints_skippped)
        return data
    elif source=='raw file':
        data= raw_data_Load_and_processing_file(path_to_raw_file, microscopy_initital_delay, microscopy_interval, initital_timepoints_skippped)
        plate= pd.read_excel(path_to_plate_file)
        data= data.merge(plate, how= 'left', on='Well')
        return data
    elif source=='processed file':
        data= procesed_data_Load_file(path_to_processed_file, initital_timepoints_skippped)
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
    
    #only exposed data
    data= data.loc[data.Conditions=='0.5 mM As']
    
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


# __visualisation functions__
def Figure_1(data, single_timepoints, stage_bins, export= False):
    fig= plt.figure(figsize= (12.5, 10), constrained_layout= True)
    gs= gridspec.GridSpec(2,2, figure=fig)
    
    #gridspec
    ax= fig.add_subplot(gs[0:1, :])
    ax1= fig.add_subplot(gs[1:2, 0:1])
    ax2= fig.add_subplot(gs[1:2, 1:2])

    #data split
    control_data= data.loc[data.Conditions=='control']
    exposed_data= data.loc[data.Conditions!='control']

    #assign stages
    exposed_data= exposed_data.assign(Stage= pd.cut(exposed_data.TimepointMinutes, 
                                                    bins= stage_bins,
                                                    labels= ['Formation', 'Relocation & Fusion', 'Clearance']))

    #subplot 1
    ax.errorbar(control_data.TimepointMinutes,
                control_data.PercentageOfCellsContainingAggregatesMean,
                lw= 2.5,
                label= 'control cells',
                color= '#00527C')
    y_min_ctrl= control_data.PercentageOfCellsContainingAggregatesMean - control_data.PercentageOfCellsContainingAggregatesMOE95
    y_max_ctrl= control_data.PercentageOfCellsContainingAggregatesMean + control_data.PercentageOfCellsContainingAggregatesMOE95
    ax.fill_between(control_data.TimepointMinutes, 
                    y_min_ctrl, 
                    y_max_ctrl, 
                    color='#DFE9F5', 
                    label="margin of error area: control cells")

    ax.errorbar(exposed_data.TimepointMinutes,
                exposed_data.PercentageOfCellsContainingAggregatesMean,
                lw= 2.5,
                label= 'As-exposed cells',
                color= '#FF781F')
    y_min_ex= exposed_data.PercentageOfCellsContainingAggregatesMean - exposed_data.PercentageOfCellsContainingAggregatesMOE95
    y_max_ex= exposed_data.PercentageOfCellsContainingAggregatesMean + exposed_data.PercentageOfCellsContainingAggregatesMOE95
    ax.fill_between(exposed_data.TimepointMinutes, 
                    y_min_ex, 
                    y_max_ex, 
                    color='#FEE7CC', 
                    label="margin of error area: exposed cells")

    ax.set_ylim(-15, 115)
    ax.set_xlabel('time (min)', weight= 'bold')
    ax.set_ylabel('percentage', weight= 'bold')
    ax.legend(frameon= False,
              ncol= 2)

    # ax.set_xticks([0,1,2,3,4,5,6,7,8,9])
    ax.set_xticks([0,60,120,180,240,300,360,420,480,540])
    

    #subplot 2
    selected_timepoints_control= control_data.loc[control_data.TimepointMinutes.isin(single_timepoints), ['TimepointHours', 'TimepointMinutes', 'PercentageOfCellsContainingAggregates', 'PercentageOfCellsContainingAggregatesMean', 'PercentageOfCellsContainingAggregatesSTD', 'PercentageOfCellsContainingAggregatesMOE95']]
    selected_timepoints_control= selected_timepoints_control.assign(PreviousTimepoint= selected_timepoints_control.PercentageOfCellsContainingAggregates.shift(1))
    selected_timepoints_control= selected_timepoints_control.assign(p_value= selected_timepoints_control.apply(lambda x: single_t_test(x['PercentageOfCellsContainingAggregates'], x['PreviousTimepoint']), axis= 1),
                                                                    significance= lambda x: np.where(x['p_value']<0.001, '*', ''))

    selected_timepoints_exposed= exposed_data.loc[exposed_data.TimepointMinutes.isin(single_timepoints), ['TimepointHours', 'TimepointMinutes', 'PercentageOfCellsContainingAggregates', 'PercentageOfCellsContainingAggregatesMean', 'PercentageOfCellsContainingAggregatesSTD', 'PercentageOfCellsContainingAggregatesMOE95']]
    selected_timepoints_exposed= selected_timepoints_exposed.assign(PreviousTimepoint= selected_timepoints_exposed.PercentageOfCellsContainingAggregates.shift(1))
    selected_timepoints_exposed= selected_timepoints_exposed.assign(p_value= selected_timepoints_exposed.apply(lambda x: single_t_test(x['PercentageOfCellsContainingAggregates'], x['PreviousTimepoint']), axis= 1),
                                                                    significance= lambda x: np.where(x['p_value']<0.001, '*', ''))

    width= .4
    x= np.arange(0, len(selected_timepoints_control))

    ax1.bar(x-width/2,
            selected_timepoints_control.PercentageOfCellsContainingAggregatesMean,
            yerr= selected_timepoints_control.PercentageOfCellsContainingAggregatesMOE95,
            width= width,
            color= '#DFE9F5',
            edgecolor= '#00527C',
            capsize= 2,
            error_kw= {'elinewidth':.75},
            lw= 1,
            label= 'control cells')

    ax1.bar(x+width/2,
            selected_timepoints_exposed.PercentageOfCellsContainingAggregatesMean,
            yerr= selected_timepoints_exposed.PercentageOfCellsContainingAggregatesMOE95,
            width= width,
            color= '#FEE7CC',
            edgecolor= '#FF781F',
            capsize= 2,
            error_kw= {'elinewidth':.75},
            lw= 1,
            label= 'As-exposed cells')

    x_coor_ctrl= x - width/2 + 0.02 ###
    y_coor_ctrl= selected_timepoints_control.PercentageOfCellsContainingAggregatesMean + selected_timepoints_control.PercentageOfCellsContainingAggregatesMOE95 + 4
    coordinates_ctrl= [[x, y] for x, y in zip(x_coor_ctrl, y_coor_ctrl)]
    for i, c in enumerate(coordinates_ctrl):
        cx,cy = c[0], c[1]
        ax1.text(cx, cy, f'{round(selected_timepoints_control.p_value, 4).fillna("").iloc[i]}', ha= 'center', rotation = 90)

    x_coor_exp= x + width/2 + 0.02 ###353♣
    y_coor_exp= selected_timepoints_exposed.PercentageOfCellsContainingAggregatesMean + selected_timepoints_exposed.PercentageOfCellsContainingAggregatesMOE95 + 4
    coordinates_exp= [[x, y] for x, y in zip(x_coor_exp, y_coor_exp)]
    for i, c in enumerate(coordinates_exp):
        cx,cy = c[0], c[1]
        ax1.text(cx, cy, f'{round(selected_timepoints_exposed.p_value, 4).fillna("").iloc[i]}', ha= 'center', rotation = 90)

    ax1.set_xticks(x)
    ax1.set_xticklabels(single_timepoints)
    ax1.axhline(0, lw=.5, ls= '--', color= 'gray')
    ax1.set_ylim(-15, 115)
    ax1.set_xlabel('timepoint (min)', weight= 'bold')
    ax1.set_ylabel('percentage', weight= 'bold')
    ax1.legend(frameon= False, ncol= 1)


    #subplot 3
    cc= round(exposed_data.loc[exposed_data.Stage=='Formation', 'AverageNumberOfAggregatesPerCellMean'].corr(exposed_data.loc[exposed_data.Stage=='Formation', 'AverageSizeOfSingleAggregatesMean'], method= 'pearson', min_periods= 2), 2)
    ax2.scatter(exposed_data.loc[exposed_data.Stage=='Formation', 'AverageNumberOfAggregatesPerCellMean'],
                exposed_data.loc[exposed_data.Stage=='Formation', 'AverageSizeOfSingleAggregatesMean'],
                edgecolor= 'black',
                lw=.15,
                color= '#DFE9F5',
                label= f'formation, cc: {cc}')
    x= exposed_data.loc[exposed_data.Stage=='Formation', 'AverageNumberOfAggregatesPerCellMean']
    y= exposed_data.loc[exposed_data.Stage=='Formation', 'AverageSizeOfSingleAggregatesMean']
    m, b= np.polyfit(x, y, 1)
    ax2.plot(x, m*x+b,
             color= '#00527C',
             linewidth= 2.5)

    cc= round(exposed_data.loc[exposed_data.Stage=='Relocation & Fusion', 'AverageNumberOfAggregatesPerCellMean'].corr(exposed_data.loc[exposed_data.Stage=='Relocation & Fusion', 'AverageSizeOfSingleAggregatesMean'], method= 'pearson', min_periods= 2), 2)
    ax2.scatter(exposed_data.loc[exposed_data.Stage=='Relocation & Fusion', 'AverageNumberOfAggregatesPerCellMean'],
                exposed_data.loc[exposed_data.Stage=='Relocation & Fusion', 'AverageSizeOfSingleAggregatesMean'],
                edgecolor= 'black',
                lw=.15,
                color= '#FEE7CC',
                label= f'relocation & fusion, cc: {cc}')
    x= exposed_data.loc[exposed_data.Stage=='Relocation & Fusion', 'AverageNumberOfAggregatesPerCellMean']
    y= exposed_data.loc[exposed_data.Stage=='Relocation & Fusion', 'AverageSizeOfSingleAggregatesMean']
    m, b= np.polyfit(x, y, 1)
    ax2.plot(x, m*x+b,
             color= '#FF781F',
             linewidth= 2.5)

    cc= round(exposed_data.loc[exposed_data.Stage=='Clearance', 'AverageNumberOfAggregatesPerCellMean'].corr(exposed_data.loc[exposed_data.Stage=='Clearance', 'AverageSizeOfSingleAggregatesMean'], method= 'pearson', min_periods= 2), 2)
    ax2.scatter(exposed_data.loc[exposed_data.Stage=='Clearance', 'AverageNumberOfAggregatesPerCellMean'],
                exposed_data.loc[exposed_data.Stage=='Clearance', 'AverageSizeOfSingleAggregatesMean'],
                edgecolor= 'black',
                lw=.15,
                color= '#D6FED2',
                label= f'clearance, cc: {cc}')
    x= exposed_data.loc[exposed_data.Stage=='Clearance', 'AverageNumberOfAggregatesPerCellMean']
    y= exposed_data.loc[exposed_data.Stage=='Clearance', 'AverageSizeOfSingleAggregatesMean']
    m, b= np.polyfit(x, y, 1)
    ax2.plot(x, m*x+b,
             color= '#0C8001',
             linewidth= 2.5)
    ax2.set_xlabel('avg. no. of agg. per cell', weight= 'bold')
    ax2.set_ylabel('avg. size of a single agg.', weight= 'bold')
    ax2.legend(frameon= True,
               loc= 'lower right')
    
    
    #export
    if export== True:
        plt.savefig(path_for_export, dpi= 1000)
    elif export== False:
        pass;
    else:
        raise ValueError(f"Invalid export argument: '{export}'. Expected: boolean ('True' or 'False').")


# ------------------------------------------------------------
# __Figure 1: WT analysis__
_20250106= data_load('db')
_20250106= missing_values(_20250106)
_20250106= repeats_group_mean_std_moe95(_20250106)
calculated_stage_bins= stage_bins(_20250106)
# _20250106.head()

# * __Figure 1__
Figure_1(_20250106, [14, 35, 70, 119, 182, 238, 301, 357, 427, 483, 539], calculated_stage_bins, export= False)
