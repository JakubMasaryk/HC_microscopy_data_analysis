#!/usr/bin/env python
# coding: utf-8

# __libraries__

# In[3]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from sqlalchemy import create_engine


# __params__

# In[5]:


plt.rcParams["legend.frameon"] = False
plt.rcParams['legend.fontsize'] = 15

plt.rcParams['axes.labelsize'] = 20
plt.rcParams['axes.labelweight'] = 'bold'

plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18  

plt.rcParams['figure.dpi'] = 1000


# __inputs__

# In[7]:


#microscopy parameters
initial_delay= 7
frequency= 3.5

#timepoints to analyse
selected_timepont= 315 #late relocation & fusion stage
reference_timepoint= 98 #late formation stage

#mysql db access (fill in accordingly, or fill directly when calling 'load_from_sql_db')
mysql_username= 'root'
mysql_password= 'poef.qve5353'
mysql_hostname= '127.0.0.1'
mysql_port= '3306'

#path to the raw file (if applicable)
path_to_raw_file= r"C:\Users\Jakub\Desktop\figures\Figure_2\data\Fig2_data.csv"


# ----------------------------------------------------------------------------------------------------------------

# __data load and process functions__

# * __from mysql database__

# In[11]:


def load_from_sql_db(username= mysql_username, password=mysql_password, hostname= mysql_hostname, port= mysql_port, ref_tmpt=reference_timepoint, sel_tmpt= selected_timepont):
    
    #connection
    connection_string = f"mysql+pymysql://{username}:{password}@{hostname}:{port}/hc_microscopy_data_v2"
    engine = create_engine(connection_string) 
    
    #query (stored procedure call) + parameters
    query = "call p_number_of_foci_single_cell_data_two_timepoints (%s, %s)"
    param1= ref_tmpt
    param2= sel_tmpt
    
    #data load
    data= pd.read_sql(query, engine, params= (param1, param2))
    
    #fix '\r' suffix
    data= data.assign(metal_concentration_unit= data.metal_concentration_unit.apply(lambda x: x.rstrip('\r')))
    
    #removing cells with 0 foci
    data= data.loc[(data.number_of_foci > 0)]
    
    #filtering down to As-exposed cells
    data= data.loc[data.metal_concentration>0]

    #assigning the timepoint category
    data= data.assign(timepoint_category=np.where(data.timepoint_minutes== reference_timepoint, 'reference timepoint', 'selected timepoint'))
    
    #sorting and reindexing 
    data= data.sort_values(['timepoint_category', 'timepoint', 'experimental_well_label', 'fov_cell_id'])
    data= data.reset_index(drop= True)
    
    #unifying dtypes
    data=data.astype({'experimental_well_label':'category', 'timepoint':'int32', 'timepoint_minutes':'float16', 'number_of_foci':'int32', 'tested_metal':'category', 'metal_concentration':'float16', 'metal_concentration_unit':'category', 'timepoint_category':'category'})
    
    return data


# * __from raw file__

# In[13]:


def load_from_file(path =path_to_raw_file,  ref_tmpt=reference_timepoint, sel_tmpt= selected_timepont, init_del= initial_delay, freq= frequency):
    
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
    data= data.loc[data.timepoint_minutes.isin([ref_tmpt, sel_tmpt])]
    
    #removing entries (cells) with zero foci
    data= data.loc[data.number_of_foci>0]
    
    #assigning the timepoint category
    data= data.assign(timepoint_category=np.where(data.timepoint_minutes== ref_tmpt, 'reference timepoint', 'selected timepoint'))
    
    #sorting and reindexing 
    data= data.sort_values(['timepoint_category', 'timepoint', 'experimental_well_label', 'fov_cell_id'])
    data= data.reset_index(drop= True)
    
    #unifying dtypes
    data=data.astype({'experimental_well_label':'category', 'timepoint':'int32', 'timepoint_minutes':'float16', 'number_of_foci':'int32', 'tested_metal':'category', 'metal_concentration':'float16', 'metal_concentration_unit':'category', 'timepoint_category':'category'})
    
    return data


# In[14]:


# load_from_file().equals(load_from_sql_db())


# * __load__

# In[16]:


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

# In[19]:


def single_cell_data_foci_count_kde(dataset, ref_tmpt=reference_timepoint, sel_tmpt= selected_timepont, export= False):

    #plot
    fig, ax= plt.subplots(figsize= (9.8, 7.2))
    sns.kdeplot(data=dataset, 
                x='number_of_foci', 
                fill=True,  
                alpha=0.25,
                hue= 'timepoint_category',
                ax=ax,
                palette= ['#00527C', '#FF781F'],
                common_norm= True)
    
    #custom legend
    legend_patches = [mpatches.Patch(color='#00527C', alpha= 0.25, label= f'late formation stage ({ref_tmpt} min)'),
                      mpatches.Patch(color='#FF781F', alpha= 0.25, label= f'late relocation & fusion stage ({sel_tmpt} min)')]
    ax.legend(handles=legend_patches)

    #plot parameters
    ax.set_xlabel('no. of agg. per cell')
    ax.set_ylabel('density')
    ax.set_ylim(0, 0.35)
    ax.set_xlim(-1, 11)
    ax.set_xticks(np.arange(0, 11))
    
    #export
    if export== True:
        plt.savefig(r"C:\Users\Jakub\Desktop\Figure2.png", bbox_inches='tight')
    elif export== False:
        pass;
    else:
        raise ValueError(f"Invalid export argument: '{export}'. Expected: boolean ('True' or 'False').")


# --------------------------------------------------------------------------------------------

# __data load__

# In[22]:


#specify data source: 'db' or 'raw file'
data= data_load(source= 'db')


# In[23]:


data.head()


# ------------------------------------------------------------------------

# __Figure 2__

# In[26]:


single_cell_data_foci_count_kde(dataset= data, export= False)


# In[ ]:




