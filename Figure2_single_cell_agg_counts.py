# __READ ME__

# __data source__
# - __relational database__
#     - relational database and needs to be created and data uploaded prior to data load (see: https://github.com/JakubMasaryk/HC_microscopy_database)
#     - fill in the mysql server parameters ('username', 'password', 'hostname', 'port') in the section bellow
#     - use argument 'db' for 'data_load function' bellow
# - __raw file (prefered method__
#     - define the pathway to the raw file ('path_to_raw_file') 
#     - file part of published supp. material
#     - use argument 'raw file' for 'data_load function' bellow

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
import matplotlib.patches as mpatches
import seaborn as sns
from sqlalchemy import create_engine
import statsmodels.api as sm
from statsmodels.miscmodels.ordinal_model import OrderedModel


# -----------------------------------------------------------------------------------------------------------------------------------
# __params__
plt.rcParams["legend.frameon"] = False
plt.rcParams['legend.fontsize'] = 15
plt.rcParams['axes.labelsize'] = 20
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18  
plt.rcParams['figure.dpi'] = 1000


# -----------------------------------------------------------------------------------------------------------------------------------
# __inputs__
#microscopy parameters
initial_delay= 7
frequency= 3.5

#timepoints to analyse (ascending order)
selected_timepoint_1= 70 #late-formation stage
selected_timepoint_2= 147 #early-relocation & fusion stage
selected_timepoint_3= 224 #mid-relocation & fusion stage
selected_timepoint_4= 301 #late-relocation & fusion stagee

#mysql db access (fill in accordingly, or fill directly when calling 'load_from_sql_db')
mysql_username= ''
mysql_password= ''
mysql_hostname= ''
mysql_port= ''

#path to the raw file (if applicable)
path_to_raw_file= r"...\Fig2_data.csv"

#path for figure export
path_for_export= r"...\Figure2.png"


# ----------------------------------------------------------------------------------------------------------------
# __data load and process functions__

# * __from mysql database__
def load_from_sql_db(username= mysql_username, password=mysql_password, hostname= mysql_hostname, port= mysql_port, sel_tmpt1=selected_timepoint_1, sel_tmpt2= selected_timepoint_2, sel_tmpt3= selected_timepoint_3, sel_tmpt4= selected_timepoint_4):
    
    #connection
    connection_string = f"mysql+pymysql://{username}:{password}@{hostname}:{port}/hc_microscopy_data_v2"
    engine = create_engine(connection_string) 
    
    #query (stored procedure call) + parameters
    query = "call p_number_of_foci_single_cell_data_two_timepoints (%s, %s, %s, %s)"
    param1= sel_tmpt1
    param2= sel_tmpt2
    param3= sel_tmpt3
    param4= sel_tmpt4
    
    #data load
    data= pd.read_sql(query, engine, params= (param1, param2, param3, param4))
    
    #fix '\r' suffix
    data= data.assign(metal_concentration_unit= data.metal_concentration_unit.apply(lambda x: x.rstrip('\r')))
    
    #filtering the foci counts between 1-6 (incl.)
    data= data.loc[(data.number_of_foci > 0) & (data.number_of_foci <= 6)]
    
    #filtering down to As-exposed cells
    data= data.loc[data.metal_concentration>0]

    #assigning the timepoint category
    mapping_dict= {sel_tmpt1:'timepoint 1',
                   sel_tmpt2: 'timepoint 2',
                   sel_tmpt3: 'timepoint 3',
                   sel_tmpt4: 'timepoint 4'}
    data= data.assign(timepoint_category= data.timepoint_minutes.map(mapping_dict))
    
    #sorting and reindexing 
    data= data.sort_values(['timepoint_category', 'timepoint', 'experimental_well_label', 'fov_cell_id'])
    data= data.reset_index(drop= True)
    
    #unifying dtypes
    data=data.astype({'experimental_well_label':'category', 'timepoint':'int32', 'timepoint_minutes':'float16', 'number_of_foci':'int32', 'tested_metal':'category', 'metal_concentration':'float16', 'metal_concentration_unit':'category', 'timepoint_category':'category'})
    
    return data

# * __from raw file__
def load_from_file(path =path_to_raw_file, sel_tmpt1=selected_timepoint_1, sel_tmpt2= selected_timepoint_2, sel_tmpt3= selected_timepoint_3, sel_tmpt4= selected_timepoint_4, init_del= initial_delay, freq= frequency):
    
    #data load
    data= pd.read_csv(path,
                      converters= {'WELL LABEL':lambda x: x.replace(' - ', '0') if len(x) == 5 else x.replace(' - ', '')},
                      usecols= ['WELL LABEL', 'T', 'FOV', 'OBJECT ID', 'Granules Org per Cell wv2'])
    
    #new columns
    data= data.assign(fov_cell_id= data['FOV'].astype('str') + '-' + data['OBJECT ID'].astype('str'),
                      tested_metal= 'As',
                      metal_concentration= 0.5,
                      metal_concentration_unit= 'mM',
                      timepoint_minutes= data['T'] * freq - (freq - init_del))
    
    #filtering columns
    data= data.loc[:, ['WELL LABEL', 'T', 'Granules Org per Cell wv2', 'fov_cell_id', 'tested_metal', 'metal_concentration', 'metal_concentration_unit', 'timepoint_minutes']]
    
    #renaming and reordering columns
    data.columns= ['experimental_well_label', 'timepoint', 'number_of_foci', 'fov_cell_id', 'tested_metal', 'metal_concentration', 'metal_concentration_unit', 'timepoint_minutes']
    data= data.reindex(columns= ['experimental_well_label', 'timepoint', 'timepoint_minutes', 'fov_cell_id', 'number_of_foci', 'tested_metal', 'metal_concentration', 'metal_concentration_unit'])
    
    #filtering down to selected and reference timepoints
    data= data.loc[data.timepoint_minutes.isin([sel_tmpt1, sel_tmpt2, sel_tmpt3, sel_tmpt4])]
    
    #filtering the foci counts between 1-6 (incl.)
    data= data.loc[(data.number_of_foci > 0) & (data.number_of_foci <= 6)]
    
    #assigning the timepoint category
    mapping_dict= {sel_tmpt1:'timepoint 1',
                   sel_tmpt2: 'timepoint 2',
                   sel_tmpt3: 'timepoint 3',
                   sel_tmpt4: 'timepoint 4'}
    data= data.assign(timepoint_category= data.timepoint_minutes.map(mapping_dict))
    
    #sorting and reindexing 
    data= data.sort_values(['timepoint_category', 'timepoint', 'experimental_well_label', 'fov_cell_id'])
    data= data.reset_index(drop= True)
    
    #unifying dtypes
    data=data.astype({'experimental_well_label':'category', 'timepoint':'int32', 'timepoint_minutes':'float16', 'number_of_foci':'int32', 'tested_metal':'category', 'metal_concentration':'float16', 'metal_concentration_unit':'category', 'timepoint_category':'category'})
    
    return data


# -----------------------------------------------------------------------------------------------------------------------------------
# * __load__
def data_load(source):
    if source=='db':
        data= load_from_sql_db()
        return data
    elif source=='raw file':
        data= load_from_file()
        return data
    else:
        raise ValueError(f"Invalid source input: '{source}'. Expected: 'db' or 'raw file'.")


# -----------------------------------------------------------------------------------------------------------------------------------
# __data visualisation function__
#agg counts per single cells: KDE plots
def single_cell_data_foci_count_kde(dataset, sel_tmpt1=selected_timepoint_1, sel_tmpt2= selected_timepoint_2, sel_tmpt3= selected_timepoint_3, sel_tmpt4= selected_timepoint_4, export= False):

    #plot
    fig, ax= plt.subplots(2, 2, figsize= (9.8, 7.2), constrained_layout= True, sharex= 'all', sharey= 'all')
    
    lf= data.loc[data.timepoint_minutes == sel_tmpt1]
    sns.kdeplot(data=lf, 
                x='number_of_foci', 
                fill=True,  
                alpha=0.25,
                ax=ax[0][0],
                color= '#00527C')
    
    erf= data.loc[data.timepoint_minutes == sel_tmpt2]
    sns.kdeplot(data=erf, 
                x='number_of_foci', 
                fill=True,  
                alpha=0.25,
                ax=ax[0][1],
                color= '#FF781F')
    
    mrf= data.loc[data.timepoint_minutes == sel_tmpt3]
    sns.kdeplot(data=mrf, 
                x='number_of_foci', 
                fill=True,  
                alpha=0.25,
                ax=ax[1][0],
                color= '#BDE7BD')
    
    lrf= data.loc[data.timepoint_minutes == sel_tmpt4]
    sns.kdeplot(data=lrf, 
                x='number_of_foci', 
                fill=True,  
                alpha=0.25,
                ax=ax[1][1],
                color= '#FF6962',
                bw_adjust=1.25)
   
    
    #custom legend
    legend_patches_lf = [mpatches.Patch(color='#00527C', alpha= 0.25, label= f'LF ({sel_tmpt1} min)')]
    legend_patches_erf = [mpatches.Patch(color='#FF781F', alpha= 0.25, label= f'ERF ({sel_tmpt2} min)')]
    legend_patches_mrf = [mpatches.Patch(color='#BDE7BD', alpha= 0.25, label= f'MRF ({sel_tmpt3} min)')]
    legend_patches_lrf = [mpatches.Patch(color='#FF6962', alpha= 0.25, label= f'LRF ({sel_tmpt4} min)')]
    
    
    ax[0][0].legend(handles=legend_patches_lf)
    ax[0][1].legend(handles=legend_patches_erf)
    ax[1][0].legend(handles=legend_patches_mrf)
    ax[1][1].legend(handles=legend_patches_lrf)
    

    #plot parameters
    ax[0][0].set_xlabel('no. of agg. per cell')
    ax[0][0].set_ylabel('density')
    ax[0][0].set_ylim(0, 1)
    ax[0][0].set_xlim(-1, 7)
    ax[0][0].set_xticks(np.arange(0, 7))
    
    ax[1][0].set_xlabel('no. of agg. per cell')
    ax[1][0].set_ylabel('density')
    
    ax[1][1].set_xlabel('no. of agg. per cell')
    
    
    
    #export
    if export== True:
        plt.savefig(path_for_export, bbox_inches='tight')
    elif export== False:
        pass;
    else:
        raise ValueError(f"Invalid export argument: '{export}'. Expected: boolean ('True' or 'False').")
        
        
        
#agg cunts per single cells: heatmap
def single_cell_data_foci_count_heatmap(data, selected_stage_olr, sel_tmpt1=selected_timepoint_1, sel_tmpt2= selected_timepoint_2, sel_tmpt3= selected_timepoint_3, sel_tmpt4= selected_timepoint_4, export= False):
    
    #calculate the p-value using OLR (comparing ref group LF to selected stage (prefferably LRF))
    p_val= ordinal_log_regression_single_stage(data, selected_stage= selected_stage_olr)
    
    #label selected timepoints
    labels_dict= {sel_tmpt1: f'LF', # late-formation
                  sel_tmpt2: f'ERF', #early-relocation & fusion
                  sel_tmpt3: f'MRF', #mid-relocation & fusion
                  sel_tmpt4: f'LRF'} #late-relocation & fusion
    data= data.assign(timepoint_label= data.timepoint_minutes.map(labels_dict))
    
    #group and count
    data_counts= data.groupby(['timepoint_label', 'timepoint_minutes', 'number_of_foci'])[['fov_cell_id']].count().reset_index()
    
    #pivot
    data_counts= data_counts.pivot_table(index= ['timepoint_label', 'timepoint_minutes'],
                                         columns= 'number_of_foci',
                                         values= 'fov_cell_id')
    
    #calculate proportions
    data_percentage= data_counts.apply(lambda x: round(x/x.sum(),2), axis= 1).reset_index()
    
    #order
    data_percentage= data_percentage.sort_values('timepoint_minutes')
    
    #filter out sorting columns
    data_percentage= data_percentage.loc[:, ['timepoint_label', 1, 2, 3, 4, 5, 6]]
    
    #set index
    data_percentage= data_percentage.set_index('timepoint_label')
    
    # visualising
    fig, ax= plt.subplots(figsize=(6.25, 5))
    sns.heatmap(data_percentage,
                cmap= 'coolwarm',                
                annot= True,
                fmt= '.2f',
                ax= ax,
                linecolor= 'white',
                linewidths= 1,
                annot_kws={"size": 15},
                cbar= False)
    ax.set_xlabel('no. of agg.')
    ax.set_ylabel('stage')
    ax.set_yticklabels(data.timepoint_label.unique(), rotation= 0)
    
    #significance symbol
    if p_val >= 0.05:
        ax.text(1.075, 0.5, '-', va='center', ha='left', rotation= 90, transform=ax.transAxes, clip_on=False, fontsize= 20)
    elif (p_val < 0.05) & (p_val >= 0.01):
        ax.text(1.075, 0.5, '*', va='center', ha='left', rotation= 90, transform=ax.transAxes, clip_on=False, fontsize= 20)
    elif (p_val < 0.01) & (p_val >= 0.001):
        ax.text(1.075, 0.5, '**', va='center', ha='left', rotation= 90, transform=ax.transAxes, clip_on=False, fontsize= 20)
    else:
        ax.text(1.075, 0.5, '***', va='center', ha='left', rotation= 90, transform=ax.transAxes, clip_on=False, fontsize= 20)
    
    #export
    if export== True:
        plt.savefig(path_for_export, bbox_inches='tight')
    elif export== False:
        pass;
    else:
        raise ValueError(f"Invalid export argument: '{export}'. Expected: boolean ('True' or 'False').")  
        

#agg counts per single cells: barchart with cumulative percentage
def single_cell_data_foci_count_barchart(data, sel_tmpt1=selected_timepoint_1, sel_tmpt2= selected_timepoint_2, sel_tmpt3= selected_timepoint_3, sel_tmpt4= selected_timepoint_4, export= False):
    
    #label selected timepoints
    labels_dict= {sel_tmpt1: f'LF', # late-formation
                  sel_tmpt2: f'ERF', #early-relocation & fusion
                  sel_tmpt3: f'MRF', #mid-relocation & fusion
                  sel_tmpt4: f'LRF'} #late-relocation & fusion
    data= data.assign(timepoint_label= data.timepoint_minutes.map(labels_dict))
    
    #group and count
    data_counts= data.groupby(['timepoint_label', 'timepoint_minutes', 'number_of_foci'])[['fov_cell_id']].count().reset_index()
    
    #pivot
    data_counts= data_counts.pivot_table(index= ['timepoint_label', 'timepoint_minutes'],
                                         columns= 'number_of_foci',
                                         values= 'fov_cell_id')
    
    #calculate proportions
    data_percentage= data_counts.apply(lambda x: round(x/x.sum(),2), axis= 1).reset_index()
    
    #order
    data_percentage= data_percentage.sort_values('timepoint_minutes')
    
    #filter out sorting columns
    data_percentage= data_percentage.loc[:, ['timepoint_label', 1, 2, 3, 4, 5, 6]]
    
    #set index + transpose
    data_percentage= data_percentage.set_index('timepoint_label').T
    
    # visualising
    x= np.arange(1, len(data_percentage)+1)
    width= 0.15
    
    #cumulative percentage
    data_percentage= data_percentage.assign(LF_cumsum= data_percentage.LF.cumsum(axis= 0),
                                            ERF_cumsum= data_percentage.ERF.cumsum(axis= 0),
                                            MRF_cumsum= data_percentage.MRF.cumsum(axis= 0),
                                            LRF_cumsum= data_percentage.LRF.cumsum(axis= 0),)
    
    fig, ax= plt.subplots(figsize= (12, 8))
    ax1= ax.twinx()
    
    ax.bar(x-1.5*width,
           data_percentage.LF*100,
           width=width,
           color= '#00527C',
           alpha= .5,
           edgecolor= 'black',
           linewidth= .75,
           label= f'LF ({sel_tmpt1} min)')
    ax1.plot(x,
            data_percentage.LF_cumsum*100,
            color= '#00527C',
            linewidth= 2.5,
            alpha= .5,
            label= 'LF')
    ax1.plot(x,
            data_percentage.LF_cumsum*100,
            linestyle= 'None',
            marker= 'v',
            markerfacecolor='white',       
            markeredgecolor='#00527C',        
            markeredgewidth=2,
            zorder=2)
    ax.bar(x-0.5*width,
           data_percentage.ERF*100,
           width=width,
           color= '#FF781F',
           alpha= .5,
           edgecolor= 'black',
           linewidth= .75,
           label= f'ERF ({sel_tmpt2} min)')
    ax.bar(x+0.5*width,
           data_percentage.MRF*100,
           width=width,
           color= '#BDE7BD',
           alpha= .5,
           edgecolor= 'black',
           linewidth= .75,
           label= f'MRF ({sel_tmpt3} min)')
    ax.bar(x+1.5*width,
           data_percentage.LRF*100,
           width=width,
           color= '#FF6962',
           alpha= .5,
           edgecolor= 'black',
           linewidth= .75,
           label= f'LRF ({sel_tmpt4} min)')
    ax1.plot(x,
            data_percentage.LRF_cumsum*100,
            color= '#FF6962',
            linewidth= 2.5,
            alpha= .5,
            label= 'LRF')
    ax1.plot(x,
            data_percentage.LRF_cumsum*100,
            linestyle= 'None',
            marker= 'o',
            markerfacecolor='white',       
            markeredgecolor='#FF6962',        
            markeredgewidth=2,
            zorder=2)
    
    ax.set_ylim(0, 115)
    ax1.set_ylim(0, 115)
    ax.set_xticks(x)
    ax.set_xticklabels(data_percentage.index)
    ax.legend(frameon= False, ncol= 2, loc= 'upper left')
    ax.set_xlabel('no. of agg. per cell')
    ax.set_ylabel('perc. of cells')
    ax1.set_ylabel('cumulative perc. of cells')
    ax1.legend(frameon= False, loc= 'upper right')
    
    #export
    if export== True:
        plt.savefig(path_for_export, bbox_inches='tight')
    elif export== False:
        pass;
    else:
        raise ValueError(f"Invalid export argument: '{export}'. Expected: boolean ('True' or 'False').") 


# ----------------------------------------------------------------------------------
# __statistics functions__
#ordinal logistic regression, comparing each stage (ERF, MRF and LRF) with reference group (LF)
#return predictor coefficients and p-values
def ordinal_log_regression(data, sel_tmpt1=selected_timepoint_1, sel_tmpt2= selected_timepoint_2, sel_tmpt3= selected_timepoint_3, sel_tmpt4= selected_timepoint_4):
    
    #assigning stages
    stage_dict= {sel_tmpt1: 'LF',
                 sel_tmpt2: 'ERF',
                 sel_tmpt3: 'MRF',
                 sel_tmpt4: 'LRF'}
    data= data.assign(stage= data.timepoint_minutes.map(stage_dict))
    
    #removing unnecessary columns
    data= data.loc[:, ['stage', 'number_of_foci']]
    
    #converting stage into category
    data= data.astype({'stage':'category'})
    
    #ordering categories: putting reference (LF) on top- used as a reference group
    data= data.assign(stage= data.stage.cat.reorder_categories(['LF', 'ERF', 'MRF', 'LRF'], ordered=False))
    
    #defining the model
    model = OrderedModel(data['number_of_foci'],
                         pd.get_dummies(data['stage'], drop_first=True),  # stages as dummy variables
                         distr='logit')
    
    #fitting the model
    res = model.fit(method='bfgs')
    
    #outputs
    p_values= pd.DataFrame(res.pvalues).iloc[:3,:]
    coeff= pd.DataFrame(res.params).iloc[:3,:]
    
    outputs= coeff.merge(p_values, how= 'left', left_index= True, right_index= True)
    outputs.columns= ['coefficient', 'p_value']
    outputs= outputs.assign(coefficient= round(outputs.coefficient, 4),
                            p_value= round(outputs.p_value, 4))
    
    
    return outputs


# --------------------------------------------------------------------------------------------

# __data load__
#specify data source: 'db' or 'raw file'
data= data_load(source= 'raw file')
data.head()


# ------------------------------------------------------------------------
# __Figure 2__
single_cell_data_foci_count_kde(data)

single_cell_data_foci_count_heatmap(data)

single_cell_data_foci_count_barchart(data)


# -----------------------------------------------------------------------------------------
# __ordinal logistic_regression__
ordinal_log_regression(data)

