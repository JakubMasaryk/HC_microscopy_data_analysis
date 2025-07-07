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
from scipy.stats import chi2_contingency
import statsmodels.api as sm
from statsmodels.miscmodels.ordinal_model import OrderedModel
from sklearn.model_selection import train_test_split
from sklearn.utils import resample as data_resample
from matplotlib.colors import ListedColormap
from scipy.stats import ttest_ind


# -----------------------------------------------------------------------------------------------------
# __params__
plt.rcParams["legend.frameon"] = False
plt.rcParams['legend.fontsize'] = 15
plt.rcParams['axes.labelsize'] = 20
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18  
plt.rcParams['font.size'] = 16
plt.rcParams['figure.dpi'] = 1000


# -----------------------------------------------------------------------------------------------------
# __mysql server connection__
#mysql server connection parameters
username= ''
password= ''
hostname= ''
port= ''

#mysql server connection
connection_string = f"mysql+pymysql://{username}:{password}@{hostname}:{port}/hc_microscopy_data_v2"
engine = create_engine(connection_string) 


# -----------------------------------------------------------------------------------------------------
# __candidate genes identification function__
def candidate_multiallele_hits(min_no_of_alleles, min_no_of_stages):
    
    #data load
    query_multiallele_hits = "call p_unique_hit_stage_effect_allele_count ()"
    data_multiallele_hits= pd.read_sql(query_multiallele_hits, engine)
    
    #process gene labels
    data_multiallele_hits= data_multiallele_hits.assign(hit_standard_name= data_multiallele_hits.hit_standard_name.apply(lambda x: x.rstrip('\r')))
    
    ##filter and sort by allele count
    data_multiallele_hits= data_multiallele_hits.loc[data_multiallele_hits.allele_count >= min_no_of_alleles]
    data_multiallele_hits= data_multiallele_hits.sort_values('allele_count', ascending= False)
    
    #filter by no of stages in which given gene is a hit
    data_multiallele_hits= data_multiallele_hits.assign(no_of_stages= data_multiallele_hits.loc[:, ['decreased_formation', 'disrupted_relocation_and_fusion', 'slower_clearance']].sum(axis= 1))
    data_multiallele_hits= data_multiallele_hits.loc[data_multiallele_hits.no_of_stages >= min_no_of_stages]
    
    return data_multiallele_hits.iloc[:, 1:-1].reset_index(drop= True)

# __resampling functions__
def no_of_entries_per_strain(data):
    
    data= data.groupby('mutation')[['fov_cell_id']].count().reset_index()
    data.columns= ['strain', 'no_of_entries']
    data= data.sort_values('no_of_entries', ascending= False, ignore_index= True)
    
    return data
def avg_no_of_entries_mutant_alleles(data):
    data= data.loc[data.strain != 'wt control']
    return round(data.no_of_entries.mean())
def strains_resample(data, target_sample_size= 2000):
    
    if len(data) > target_sample_size:
        downsampled_dataset, a = train_test_split(data, 
                                                  train_size=target_sample_size, 
                                                  stratify=data.iloc[:, 0], 
                                                  random_state=123)
        return downsampled_dataset.reset_index(drop= True)
        
    elif len(data) < target_sample_size:
        upsampled_dataset = data_resample(data,
                                          replace= True,
                                          n_samples= target_sample_size,
                                          random_state= 123)
        return upsampled_dataset

    else: 
        pass


# __data analysis and visualisation functions__
def data_analysis_and_visualisation(selected_gene, min_cell_count, starting_timepoint_min, ending_timepoint_min, resample_agg_counts):
    
    #query params (common for both agg. size and no. data)
    param1= selected_gene
    param2= min_cell_count
    param3= starting_timepoint_min
    param4= ending_timepoint_min
    
    #load data on aggregate number (single cell)
    query_agg_no = "call p_single_cell_data_agg_no_all_alleles_selected_gene (%s, %s, %s, %s)"
    data_agg_no= pd.read_sql(query_agg_no, engine, params= (param1,param2,param3,param4))
    #load data on aggregate size (single cell)
    query_agg_sz = "call p_single_cell_data_agg_size_all_alleles_selected_gene (%s, %s, %s, %s)"
    data_agg_sz= pd.read_sql(query_agg_sz, engine, params= (param1,param2,param3,param4))
    
    #resampling agg. no.
    if resample_agg_counts== True:
        target_no_of_entries= avg_no_of_entries_mutant_alleles(no_of_entries_per_strain(data_agg_no))
        data_agg_no= data_agg_no.loc[:, ['mutation', 'number_of_foci']]
        data_agg_no= data_agg_no.groupby('mutation').apply(lambda x: strains_resample(x, target_no_of_entries)).reset_index(drop= True)
    
    elif resample_agg_counts== False:
        data_agg_no= data_agg_no.loc[:, ['mutation', 'number_of_foci']]
        
    else:
        raise ValueError(f"Invalid input argument: '{resample_agg_counts}'. Expected: boolean ('True' or 'False').")
    
    #calculate percentage of cell containing (+) or not containing (- aggregates)
    #based on the data for agg. no
    data_percentages= data_agg_no.assign(aggregates= np.where(data_agg_no.number_of_foci==0, '-', '+'))
    data_percentages= data_percentages.pivot_table(index= 'mutation',
                                                   columns= 'aggregates',
                                                   values= 'number_of_foci',
                                                   aggfunc= 'count')
    data_percentages= data_percentages.apply(lambda x: round(x/x.sum(),2), axis= 1)
    data_percentages= data_percentages.sort_index(ascending= False)
    
    #percentage of cells containing a specific number of aggregates, only cells containing at least 1 agg. (+)
    #based on the data for agg. no
    data_agg_counts= data_agg_no.assign(cell_count= 1)
    data_agg_counts= data_agg_counts.groupby(['mutation', 'number_of_foci'])[['cell_count']].count().reset_index()
    data_agg_counts= data_agg_counts.loc[(data_agg_counts.number_of_foci > 0) & (data_agg_counts.number_of_foci <= 6)]
    data_agg_counts= data_agg_counts.pivot_table(index= 'mutation',
                                                 columns= 'number_of_foci',
                                                 values= 'cell_count').fillna(0)
    data_agg_counts= data_agg_counts.apply(lambda x: round(x/x.sum(),2), axis= 1).sort_index(ascending= False)
    
    #statistical tests
    #chi2 for percentage of cells containing aggregates
    chi2_data= chi2_test(data_agg_no)
    #ordinal logistic regression for agg counts
    orl_data= ordinal_log_regression(data_agg_no)
    #t-test for agg. size
    data_agg_sz= data_agg_sz.loc[:, ['mutation', 'avg_focus_size']]
    data_agg_sz_control_data= data_agg_sz.loc[data_agg_sz.mutation=='wt control', 'avg_focus_size']
    data_agg_sz_mutant_data= data_agg_sz.loc[data_agg_sz.mutation!='wt control']

    data_agg_sz_p_values= pd.DataFrame(data_agg_sz_mutant_data.groupby('mutation')[['avg_focus_size']].apply(lambda x: t_test(data_agg_sz_control_data, x['avg_focus_size']))).reset_index()
    data_agg_sz_p_values.columns= ['allele', 'p_value']
    data_agg_sz_p_values= data_agg_sz_p_values.assign(p_value= data_agg_sz_p_values.p_value.apply(lambda x: round(x, 4)))
    
    #visuals
    fig= plt.figure(figsize= (17.5, 1 * len(data_agg_counts)), constrained_layout= False)
    gs= gridspec.GridSpec(len(data_agg_counts.index), 15, figure=fig)
    
    #percentage of cells containing agg.
    ax= fig.add_subplot(gs[:, 0:1])
    ax1= fig.add_subplot(gs[:, 1:2])
    ax_ch2= fig.add_subplot(gs[:, 2:3])
    #agg. counts
    ax2= fig.add_subplot(gs[:, 4:10])
    ax_orl= fig.add_subplot(gs[:, 10:12])
    
    light_gray = ListedColormap(['#F2F2F2']) 

    sns.heatmap(pd.DataFrame(data_percentages.loc[:, '-']),
                cmap= 'Greens',                
                annot= True,
                fmt= '.2f',
                ax= ax,
                linecolor= 'white',
                linewidths= 1,
                annot_kws={"size": 15},
                cbar= False)
    ax.set_xlabel('')
    ax.set_ylabel('allele')

    sns.heatmap(pd.DataFrame(data_percentages.loc[:, '+']),
                cmap= 'Reds',                
                annot= True,
                fmt= '.2f',
                ax= ax1,
                linecolor= 'white',
                linewidths= 1,
                annot_kws={"size": 15},
                cbar= False,
                yticklabels= False)
    ax1.set_xlabel('')
    ax1.set_ylabel('')
    
    sns.heatmap(chi2_data,
                cmap= light_gray,                
                annot= True,
                fmt= '.4f',
                ax= ax_ch2,
                linecolor= 'white',
                linewidths= 1,
                annot_kws={"size": 15},
                cbar= False,
                yticklabels= False)
    ax_ch2.set_ylabel('')

    sns.heatmap(data_agg_counts,
                cmap= 'Reds',                
                annot= True,
                fmt= '.2f',
                ax= ax2,
                linecolor= 'white',
                linewidths= 1,
                annot_kws={"size": 15},
                cbar= False,
                yticklabels= False)
    ax2.set_xlabel('no. of foci per cell')
    ax2.set_ylabel('')
    
    sns.heatmap(orl_data,
                cmap= light_gray,                
                annot= True,
                fmt= '.4f',
                ax= ax_orl,
                linecolor= 'white',
                linewidths= 1,
                annot_kws={"size": 15},
                cbar= False,
                yticklabels= False)
    ax_orl.set_ylabel('');

    #size dist.
    for i, strain in enumerate(sorted(list(data_agg_sz.mutation.unique()), reverse= True)):
        ax_sz= fig.add_subplot(gs[i:i+1, 13:15])
        
        color= 'red' if strain== 'wt control' else 'green'
        
        data_to_plot= data_agg_sz.loc[data_agg_sz.mutation.isin([strain]), ['avg_focus_size']]
        ax_sz.hist(data_to_plot,
                   bins= 40,
                   alpha= .5,
                   color= color,
                   edgecolor= 'black',
                   linewidth= .15,
                   density= True)
        ax_sz.set_xlabel('')
        ax_sz.set_ylabel(strain, fontsize= 9) 
        ax_sz.set_xticks([])           
        ax_sz.set_xticklabels([])
        ax_sz.set_yticks([])           
        ax_sz.set_yticklabels([])
        ax_sz.set_xlim(-0.2, 1.2)
        ax_sz.set_ylim(0, 5)
        ax_sz.spines['top'].set_visible(False)
        ax_sz.spines['right'].set_visible(False)
        
        if strain == 'wt control':
            ax_sz.text(0.75, 3.75, f'-', fontsize= 9, weight= 'bold')
        else:
            p_value= data_agg_sz_p_values.loc[data_agg_sz_p_values.allele==strain, 'p_value'].iloc[0]
            ax_sz.text(0.75, 3.75, f'{p_value}', fontsize= 9, weight= 'bold')
            ax_sz.hist(data_agg_sz.loc[data_agg_sz.mutation== 'wt control', ['avg_focus_size']],
                       bins= 40,
                       alpha= .25,
                       color= 'red',
                       edgecolor= 'black',
                       linewidth= .15,
                       density= True)


# __stats functions__
def chi2_test(data):
    
    chi2_data= data.assign(aggregates= np.where(data.number_of_foci > 0, '+', '-'))
    chi2_data= chi2_data.pivot_table(index= 'mutation',
                                     columns= 'aggregates', 
                                     values= 'number_of_foci',
                                     aggfunc= 'count').sort_index(ascending= False)

    strains= ['wt_control']
    p_values= [0]

    for act_allele in list(chi2_data.loc[chi2_data.index!= 'wt control'].index):

        partial_chi2_data= chi2_data.loc[(chi2_data.index== 'wt control') | (chi2_data.index== act_allele)]
        chi2, p, dof, expected = chi2_contingency(partial_chi2_data)

        strains.append(act_allele)
        p_values.append(p)

    stats= pd.DataFrame({'mutation':strains, 'χ2\np-value':p_values})
    stats= stats.set_index('mutation').replace(0, np.NaN)
    stats['χ2\np-value']= round(stats['χ2\np-value'] , 4)

    return stats  

#ordinal logistic regression, comparing each allel with reference group (wt control)
#return predictor coefficients and p-values
def ordinal_log_regression(data):
    
    #removing unnecessary columns
    data= data.loc[:, ['mutation', 'number_of_foci']]
    
    #converting stage into category
    data= data.astype({'mutation':'category'})
    
    #ordered list of strains
    strain_list= list(data.mutation.unique())
    strain_list.remove('wt control')
    strain_list.sort(reverse= True)
    strain_list.insert(0, 'wt control')
    
    #ordering categories: putting reference (LF) on top- used as a reference group
    data= data.assign(mutation= data.mutation.cat.reorder_categories(strain_list, ordered=False))
    
    #defining the model
    model = OrderedModel(data['number_of_foci'],
                         pd.get_dummies(data['mutation'], drop_first=True),  # stages as dummy variables
                         distr='logit')
    
    #fitting the model
    res = model.fit(method='bfgs')
    
    #outputs
    p_values= pd.DataFrame(res.pvalues).iloc[:len(data.mutation.unique())-1,:]
    coeff= pd.DataFrame(res.params).iloc[:len(data.mutation.unique())-1,:]
    
    outputs= coeff.merge(p_values, how= 'left', left_index= True, right_index= True)
    outputs.columns= ['coefficient', 'p_value']
    outputs= outputs.assign(coefficient= round(outputs.coefficient, 4),
                            p_value= round(outputs.p_value, 4))
    
    #add control (WT) row
    ctrl_row = pd.DataFrame([{'index':'wt control', 'coefficient': np.NaN, 'p_value': np.NaN}])
    outputs = pd.concat([outputs.reset_index(), ctrl_row])
    
    #rename columns
    outputs.columns= ['index', 'coeff.', 'OLR\np-value']
    
    return outputs.set_index('index').sort_index(ascending= False)

#two-sample ttest
def t_test(control_group, tested_group):
    
    t_stat, p_value = ttest_ind(control_group, tested_group, equal_var=False)
    
    return p_value


# --------------------------------------------------------------------------------------
# __candidate genes__
candidate_multiallele_hits(min_no_of_alleles= 5, min_no_of_stages= 2)
# -------------------------------------------------------------------------------------------------------------------------------------
# __ACT1 alleles__

# * __formation (full: 0-62 min)__
data_analysis_and_visualisation(selected_gene= 'ACT1', min_cell_count= 33, starting_timepoint_min= 0, ending_timepoint_min= 62, resample_agg_counts= True)

# * __relocation & fusion (full: 62-300 min)__
data_analysis_and_visualisation(selected_gene= 'ACT1', min_cell_count= 33, starting_timepoint_min= 62, ending_timepoint_min= 300, resample_agg_counts= True)

# * __relocation & fusion (late: 280-300 min)__
data_analysis_and_visualisation(selected_gene= 'ACT1', min_cell_count= 33, starting_timepoint_min= 280, ending_timepoint_min= 300, resample_agg_counts= True)

# * __clearance (full: > 300 min)__
data_analysis_and_visualisation(selected_gene= 'ACT1', min_cell_count= 33, starting_timepoint_min= 300, ending_timepoint_min= 400, resample_agg_counts= True)


# ----------------------------------------------------------
# __MOB2 alleles__

# * __formation (full: 0-62 min)__
data_analysis_and_visualisation(selected_gene= 'MOB2', min_cell_count= 33, starting_timepoint_min= 0, ending_timepoint_min= 62, resample_agg_counts= True)

# * __relocation & fusion (full: 62-300 min)__
data_analysis_and_visualisation(selected_gene= 'MOB2', min_cell_count= 33, starting_timepoint_min= 62, ending_timepoint_min= 300, resample_agg_counts= True)

# * __relocation & fusion (late: 280-300 min)__
data_analysis_and_visualisation(selected_gene= 'MOB2', min_cell_count= 33, starting_timepoint_min= 280, ending_timepoint_min= 300, resample_agg_counts= True)

# * __clearance (full: > 300 min)__
data_analysis_and_visualisation(selected_gene= 'MOB2', min_cell_count= 33, starting_timepoint_min= 300, ending_timepoint_min= 400, resample_agg_counts= True)


# ---------------------------------------------------------------------------------------
# __POL1 alleles__

# * __formation (full: 0-62 min)__
data_analysis_and_visualisation(selected_gene= 'POL1', min_cell_count= 33, starting_timepoint_min= 0, ending_timepoint_min= 62, resample_agg_counts= True)

# * __relocation & fusion (full: 62-300 min)__
data_analysis_and_visualisation(selected_gene= 'POL1', min_cell_count= 33, starting_timepoint_min= 62, ending_timepoint_min= 300, resample_agg_counts= True)

# * __relocation & fusion (late: 280-300 min)__
data_analysis_and_visualisation(selected_gene= 'POL1', min_cell_count= 33, starting_timepoint_min= 280, ending_timepoint_min= 300, resample_agg_counts= True)

# * __clearance (full: > 300 min)__
data_analysis_and_visualisation(selected_gene= 'POL1', min_cell_count= 33, starting_timepoint_min= 300, ending_timepoint_min= 400, resample_agg_counts= True)


# --------------------------------------------------------------------------------------------------------------------
# __CDC48 alleles__

# * __formation (full: 0-62 min)__
data_analysis_and_visualisation(selected_gene= 'CDC48', min_cell_count= 33, starting_timepoint_min= 0, ending_timepoint_min= 62, resample_agg_counts= True)

# * __relocation & fusion (full: 62-300 min)__
data_analysis_and_visualisation(selected_gene= 'CDC48', min_cell_count= 33, starting_timepoint_min= 62, ending_timepoint_min= 300, resample_agg_counts= True)

# * __relocation & fusion (late: 280-300 min)__
data_analysis_and_visualisation(selected_gene= 'CDC48', min_cell_count= 33, starting_timepoint_min= 280, ending_timepoint_min= 300, resample_agg_counts= True)

# * __clearance (full: > 300 min)__
data_analysis_and_visualisation(selected_gene= 'CDC48', min_cell_count= 33, starting_timepoint_min= 300, ending_timepoint_min= 400, resample_agg_counts= True)


# --------------------------------------------------------------------------------------
# __DAM1 alleles__

# * __formation (full: 0-62 min)__
data_analysis_and_visualisation(selected_gene= 'DAM1', min_cell_count= 33, starting_timepoint_min= 0, ending_timepoint_min= 62, resample_agg_counts= True)

# * __relocation & fusion (full: 62-300 min)__
data_analysis_and_visualisation(selected_gene= 'DAM1', min_cell_count= 33, starting_timepoint_min= 62, ending_timepoint_min= 300, resample_agg_counts= True)

# * __relocation & fusion (late: 280-300 min)__
data_analysis_and_visualisation(selected_gene= 'DAM1', min_cell_count= 33, starting_timepoint_min= 280, ending_timepoint_min= 300, resample_agg_counts= True)

# * __clearance (full: > 300 min)__
data_analysis_and_visualisation(selected_gene= 'DAM1', min_cell_count= 33, starting_timepoint_min= 300, ending_timepoint_min= 400, resample_agg_counts= True)


# ---------------------------------------------------------------------------------
# __PRE2 alleles__

# * __formation (full: 0-62 min)__
data_analysis_and_visualisation(selected_gene= 'PRE2', min_cell_count= 33, starting_timepoint_min= 0, ending_timepoint_min= 62, resample_agg_counts= True)

# * __relocation & fusion (full: 62-300 min)__
data_analysis_and_visualisation(selected_gene= 'PRE2', min_cell_count= 33, starting_timepoint_min= 62, ending_timepoint_min= 300, resample_agg_counts= True)

# * __relocation & fusion (late: 280-300 min)__
data_analysis_and_visualisation(selected_gene= 'PRE2', min_cell_count= 33, starting_timepoint_min= 280, ending_timepoint_min= 300, resample_agg_counts= True)

# * __clearance (full: > 300 min)__
data_analysis_and_visualisation(selected_gene= 'PRE2', min_cell_count= 33, starting_timepoint_min= 300, ending_timepoint_min= 400, resample_agg_counts= True)


# -----------------------------------------------------------------------------------
