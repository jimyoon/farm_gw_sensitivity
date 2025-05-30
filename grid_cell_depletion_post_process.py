# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
from itertools import product
import matplotlib.pyplot as plt
import os 

# Pandas show all rows/coluns
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.max_colwidth', None)

os.chdir('/Users/yoon644/Library/CloudStorage/OneDrive-PNNL/Documents/PyProjects/farm_gw_post')

NLDAS_id_aquifer_id = pd.read_csv("NLDAS_ID_aquifer_ID.csv")
nldas_id_farm_id = pd.read_csv('NLDAS_id_Farm_id_CONUS.csv')
nldas_id_farm_id.farm_id = nldas_id_farm_id.farm_id - 1 # correct to zero index for farm_id

# used later in script to map grid cell depletion indexed by farm_id to both NLDAS_ID and aquifer classsifications 
farm_id_aquifer_id = NLDAS_id_aquifer_id.merge(nldas_id_farm_id , how = 'inner', left_on = 'NLDAS_ID', right_on = 'nldas') 
farm_id_aquifer_id = farm_id_aquifer_id[farm_id_aquifer_id.columns.intersection(['Aquifer', 'BroaderSys', 'farm_id', 'nldas'])]

# import aggregated grid depletion results 
scenario_depletion = pd.read_csv("depletion_results_conus.csv", index_col = 0)

# #   correctly "zero indexed" missing farms
#
# missing_farms_one_indexed = pd.read_csv("missing_farm_ids_CONUS.csv")
# missing_farms_zero_indexed = pd.read_csv("missing_farm_ids_CONUS_zero_indexed.csv")
#
# overlapping_missing_indexed =  missing_farms_zero_indexed.merge(missing_farms_one_indexed, how = 'inner', left_on = 'Farm_id', right_on = 'Farm_id')
#
# nonoverlapping_missing_farm_ids = []
# for i in range(len(missing_farms_zero_indexed.Farm_id)):
#     if (missing_farms_zero_indexed.Farm_id[i] in missing_farms_one_indexed.Farm_id.values) == False:
#         nonoverlapping_missing_farm_ids.append(missing_farms_zero_indexed.Farm_id[i])
#     else:
#         continue
#
# nonoverlapping_ids = overlapping_ids = pd.DataFrame(data = np.array(nonoverlapping_missing_farm_ids), columns = ['farm_id'])
# nonoverlapping_ids = nonoverlapping_ids.merge(nldas_id_farm_id, how = 'inner', left_on = 'farm_id', right_on = 'farm_id')
# #nonoverlapping_ids.to_csv("missing_farm_ids_for_HPC_rerun.csv")
#
# overlapping_missing_farm_ids = []
# for i in range(len(missing_farms_zero_indexed.Farm_id)):
#     if (missing_farms_zero_indexed.Farm_id[i] in missing_farms_one_indexed.Farm_id.values) == True:
#         overlapping_missing_farm_ids.append(missing_farms_zero_indexed.Farm_id[i])
#     else:
#         continue
#
# overlapping_ids = pd.DataFrame(data = np.array(overlapping_missing_farm_ids), columns = ['farm_id'])
# overlapping_ids = overlapping_ids.merge(nldas_id_farm_id, how = 'inner', left_on = 'farm_id', right_on = 'farm_id')
# #overlapping_ids.to_csv("missing_farm_ids_with_outputs.csv")

###### Generate data for figure 1 (depletion during baseline by aquifer with 3D surface examples)

# Subset data based on various scenarios
subset_df = scenario_depletion[scenario_depletion['K_scenario'] == 'int_1']
subset_df = subset_df[subset_df['hydro_ratio'] == 1]
subset_df = subset_df[subset_df['econ_ratio'] == 1]
subset_df = subset_df[subset_df['gamma_scenario'] == 1]

# import initial grid cell volumes
volume = pd.read_csv('aquifer_grid_cell_volume 1.csv')

subset_df = scenario_depletion[scenario_depletion['K_scenario'] == 'int_1']
subset_df = subset_df[subset_df['hydro_ratio'] == 1]
subset_df = subset_df[subset_df['econ_ratio'] == 1]
subset_df = subset_df[subset_df['gamma_scenario'] == 1]
merged_df = pd.merge(subset_df, farm_id_aquifer_id, left_on='farm',right_on='farm_id',how='left')
merged_df = pd.merge(merged_df, volume[['nldas_y','gw_volume']], left_on='nldas', right_on='nldas_y', how='left')

merged_df['volume_depleted'] = merged_df['percent_depleted'] * merged_df['gw_volume']

aquifer_calcs = pd.DataFrame(columns=['Aquifer','percent_depleted'])
for aquifer in merged_df['Aquifer'].unique():
    print(aquifer)
    sum_gw_vol = merged_df[(merged_df.Aquifer==aquifer)]['gw_volume'].sum()
    sum_vol_dep = merged_df[(merged_df.Aquifer==aquifer)]['volume_depleted'].sum()
    print(sum_gw_vol)
    percent_depleted = sum_vol_dep / sum_gw_vol
    new_row = {'Aquifer': aquifer, 'percent_depleted': percent_depleted}
    aquifer_calcs.loc[len(aquifer_calcs)] = new_row

aquifer_calcs.to_csv('4GIS_depletion_weighted_baseline.csv')

# Extract 3D surface for example cells

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import random
from matplotlib import cm

merged_df = pd.merge(scenario_depletion, farm_id_aquifer_id, left_on='farm',right_on='farm_id',how='left')

cell_id = 'x47y85' # Tulare Basin / CA Central Valley
cell_id = 'x48y175' # Columbia Plateau
cell_id = 'x89y142' # Central Plain / Snake River
cell_id = 'x217y127' # Northern High Plains / High Plains
cell_id = 'x249y165' # Central Minnesota
cell_id = 'x275y86' # Mississippi Embayment

cell_df = merged_df[(merged_df['nldas']==cell_id)]
cell_df = cell_df[(cell_df['K_scenario']=='int_1')]
cell_df = cell_df[(cell_df['gamma_scenario']==1)]

x = cell_df['hydro_ratio']
y = cell_df['econ_ratio']
z = cell_df['percent_depleted']

my_cmap = plt.get_cmap('Reds')

fig = plt.figure(figsize=(5,4))
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_trisurf(x, y, z, cmap= my_cmap, vmin=0, vmax=0.85, linewidth=0.1, edgecolor='grey')

ax.set_zlim((0.0,1.0))

fig.savefig(cell_id + '.svg', bbox_inches='tight')

###### Generate data for figure 2 (depletion maps under various scenarios)

# subset desired scenario

# baseline
subset_df_baseline = scenario_depletion[scenario_depletion['K_scenario'] == 'int_1']
subset_df_baseline = subset_df_baseline[subset_df_baseline['hydro_ratio'] == 1]
subset_df_baseline = subset_df_baseline[subset_df_baseline['econ_ratio'] == 1]
subset_df_baseline = subset_df_baseline[subset_df_baseline['gamma_scenario'] == 1]
# lower water availability (DONE)
subset_df = scenario_depletion[scenario_depletion['K_scenario'] == 'int_1']
subset_df = subset_df[subset_df['hydro_ratio'] == 0.75]
subset_df = subset_df[subset_df['econ_ratio'] == 1]
subset_df = subset_df[subset_df['gamma_scenario'] == 1]
# higher water availability (DONE)
subset_df = scenario_depletion[scenario_depletion['K_scenario'] == 'int_1']
subset_df = subset_df[subset_df['hydro_ratio'] == 1.25]
subset_df = subset_df[subset_df['econ_ratio'] == 1]
subset_df = subset_df[subset_df['gamma_scenario'] == 1]
# higher crop prices (DONE)
subset_df = scenario_depletion[scenario_depletion['K_scenario'] == 'int_1']
subset_df = subset_df[subset_df['hydro_ratio'] == 1]
subset_df = subset_df[subset_df['econ_ratio'] == 1.25]
subset_df = subset_df[subset_df['gamma_scenario'] == 1]
# lower crop prices (DONE)
subset_df = scenario_depletion[scenario_depletion['K_scenario'] == 'int_1']
subset_df = subset_df[subset_df['hydro_ratio'] == 1]
subset_df = subset_df[subset_df['econ_ratio'] == 0.75]
subset_df = subset_df[subset_df['gamma_scenario'] == 1]
# lower water availability, higher crop prices (DONE)
subset_df = scenario_depletion[scenario_depletion['K_scenario'] == 'int_1']
subset_df = subset_df[subset_df['hydro_ratio'] == 0.75]
subset_df = subset_df[subset_df['econ_ratio'] == 1.25]
subset_df = subset_df[subset_df['gamma_scenario'] == 1]
# higher water availability, lower crop prices (DONE)
subset_df = scenario_depletion[scenario_depletion['K_scenario'] == 'int_1']
subset_df = subset_df[subset_df['hydro_ratio'] == 1.25]
subset_df = subset_df[subset_df['econ_ratio'] == 0.75]
subset_df = subset_df[subset_df['gamma_scenario'] == 1]

# calculate difference from baseline
subset_df_baseline.rename(columns={'percent_depleted': 'percent_depleted_baseline'}, inplace=True)
subset_df = pd.merge(subset_df, subset_df_baseline[['farm','percent_depleted_baseline']], on='farm', how='left')
subset_df['perc_diff'] = subset_df['percent_depleted'] - subset_df['percent_depleted_baseline']

# import initial grid cell volumes
volume = pd.read_csv('aquifer_grid_cell_volume 1.csv')

# merge scenario depletion results to volumes and calculate aquifer-level weighted depletion
merged_df = pd.merge(subset_df, farm_id_aquifer_id, left_on='farm',right_on='farm_id',how='left')
merged_df = pd.merge(merged_df, volume[['nldas_y','gw_volume']], left_on='nldas', right_on='nldas_y', how='left')

merged_df['volume_depleted_diff'] = merged_df['perc_diff'] * merged_df['gw_volume']
merged_df['volume_depleted'] = merged_df['percent_depleted'] * merged_df['gw_volume']

aquifer_calcs = pd.DataFrame(columns=['Aquifer','percent_depleted_diff','percent_depleted'])
for aquifer in merged_df['Aquifer'].unique():
    sum_gw_vol = merged_df[(merged_df.Aquifer==aquifer)]['gw_volume'].sum()
    sum_vol_dep_diff = merged_df[(merged_df.Aquifer==aquifer)]['volume_depleted_diff'].sum()
    sum_vol_dep = merged_df[(merged_df.Aquifer == aquifer)]['volume_depleted'].sum()
    percent_depleted_diff = sum_vol_dep_diff / sum_gw_vol
    percent_depleted = sum_vol_dep / sum_gw_vol
    new_row = {'Aquifer': aquifer, 'percent_depleted_diff': percent_depleted_diff, 'percent_depleted': percent_depleted}
    aquifer_calcs.loc[len(aquifer_calcs)] = new_row

# export results to csv
aquifer_calcs.to_csv('4GIS_depletiondiff_lowercrop.csv')

###### Generate Figure 3 (farm-gw classifications across all scenarios)

#%% Tally number of scenarios that meet depletion severity threshold (in percent)
#   for all econ and hydro scenarios and also within each econ/hydro quadrant

depletion_threshold = 0.50 # 0.2

# For testing, pretend that 'count depletion 20m' > 0 represents when a scenario 
# produces > 20% depletion (or whatever the defined depletion threshold is)
gamma_lower = 0.9 # 0.4
gamma_upper = 1.1 # 1.6

K = 'int_1'

scenario_results = scenario_depletion[(scenario_depletion.gamma_scenario < gamma_upper) & (scenario_depletion.gamma_scenario > gamma_lower)] #  
scenario_results = scenario_results.where(scenario_results.K_scenario == K)
scenario_results = scenario_results.dropna() 
scenario_results['count_depletion_20%'] = scenario_results['percent_depleted'] > depletion_threshold
scenario_statistics = scenario_results.groupby(['farm'])[['count_depletion_20%']].sum()

# For each quadrant, filter by corresponding econ and hydro ratio bounds 

# Quadrant 1 
scenario_results = scenario_depletion[(scenario_depletion.gamma_scenario < gamma_upper) & (scenario_depletion.gamma_scenario > gamma_lower)] #    
scenario_results = scenario_results.where(scenario_results.econ_ratio > 1)
scenario_results = scenario_results.where(scenario_results.hydro_ratio < 1)
scenario_results = scenario_results.where(scenario_results.K_scenario == K)
scenario_results = scenario_results.dropna() 
scenario_results['count_depletion_20%_Q1'] = scenario_results['percent_depleted'] > depletion_threshold
scenario_statistics_Q1 = scenario_results.groupby(['farm'])[['count_depletion_20%_Q1']].sum()

# Quadrant 2
scenario_results = scenario_depletion[(scenario_depletion.gamma_scenario < gamma_upper) & (scenario_depletion.gamma_scenario > gamma_lower)] #     
scenario_results = scenario_results.where(scenario_results.econ_ratio > 1)
scenario_results = scenario_results.where(scenario_results.hydro_ratio > 1)
scenario_results = scenario_results.where(scenario_results.K_scenario == K)
scenario_results = scenario_results.dropna() 
scenario_results['count_depletion_20%_Q2'] = scenario_results['percent_depleted'] > depletion_threshold
scenario_statistics_Q2 = scenario_results.groupby(['farm'])[['count_depletion_20%_Q2']].sum()

# Quadrant 3
scenario_results = scenario_depletion[(scenario_depletion.gamma_scenario < gamma_upper) & (scenario_depletion.gamma_scenario > gamma_lower)] #    
scenario_results = scenario_results.where(scenario_results.econ_ratio < 1)
scenario_results = scenario_results.where(scenario_results.hydro_ratio < 1)
scenario_results = scenario_results.where(scenario_results.K_scenario == K)
scenario_results = scenario_results.dropna() 
scenario_results['count_depletion_20%_Q3'] = scenario_results['percent_depleted'] > depletion_threshold
scenario_statistics_Q3 = scenario_results.groupby(['farm'])[['count_depletion_20%_Q3']].sum()

# Quadrant 4
scenario_results = scenario_depletion[(scenario_depletion.gamma_scenario < gamma_upper) & (scenario_depletion.gamma_scenario > gamma_lower)] #     
scenario_results = scenario_results.where(scenario_results.econ_ratio < 1)
scenario_results = scenario_results.where(scenario_results.hydro_ratio > 1)
scenario_results = scenario_results.where(scenario_results.K_scenario == K)
scenario_results = scenario_results.dropna() 
scenario_results['count_depletion_20%_Q4'] = scenario_results['percent_depleted'] > depletion_threshold
scenario_statistics_Q4 = scenario_results.groupby(['farm'])[['count_depletion_20%_Q4']].sum()

# Baseline (hydro = 1, econ = 1)
scenario_results = scenario_depletion.where(scenario_depletion.econ_ratio == 1)
scenario_results = scenario_results.where(scenario_results.hydro_ratio == 1)
scenario_results = scenario_results.where(scenario_results.gamma_scenario == 1)
scenario_results = scenario_results.where(scenario_results.K_scenario == K)
scenario_results = scenario_results.dropna() 
scenario_results['count_depletion_baseline'] = scenario_results['percent_depleted'] > depletion_threshold
scenario_statistics_baseline = scenario_results.groupby(['farm'])[['count_depletion_baseline']].sum()

# Concatenate depletion stats results 
results = [scenario_statistics, scenario_statistics_Q1, scenario_statistics_Q2,
               scenario_statistics_Q3, scenario_statistics_Q4, scenario_statistics_baseline]

depletion_counts = pd.concat(results, axis = 1)

#%% Categorize depletion outcomes 

depletion_classification = pd.DataFrame(data = np.zeros(len(depletion_counts.iloc[:,0])), 
                                        columns = ['depletion_classification'], 
                                        index = depletion_counts.index)

gamma_scenarios = 1 # 5 # number gamma scenarios
K_scenarios = 1 # number of K scenarios 
scenario_count_total = 25 * gamma_scenarios * K_scenarios # total number of scenarios (25 = # of hydro and econ scenario combinations)
scenario_count_quadrants = 4 * gamma_scenarios * K_scenarios # total number of scenarios within each quadrant 


depletion_counts_pct = depletion_counts.copy()
depletion_counts_pct.iloc[:,0] = depletion_counts.iloc[:,0]/scenario_count_total
depletion_counts_pct.iloc[:,1:5] = depletion_counts.iloc[:,1:5]/scenario_count_quadrants


# Uncomment V2 statements for alternative classification that includes cells with compound >20%, but NOT total count >20%
for i in range(len(depletion_counts.iloc[:,0])):
    
    # V2
    # if depletion_counts_pct.iloc[i,0] < 0.2 and depletion_counts_pct.iloc[i,1] < 0.2:
    #     depletion_classification.iloc[i] = 0
    
    # V1 - comment out this "if" statement for V2 classification option 
    if depletion_counts_pct.iloc[i,0] < 0.2:
        depletion_classification.iloc[i] = 0
    
    elif depletion_counts_pct.iloc[i,0] > 0.2 and depletion_counts_pct.iloc[i,5] == 1:
        depletion_classification.iloc[i] = 1
        
    elif depletion_counts_pct.iloc[i,0] > 0.2 and depletion_counts_pct.iloc[i,1] > 0.2 \
        and depletion_counts_pct.iloc[i,2] < 0.2 and depletion_counts_pct.iloc[i,3] < 0.2:
            
        depletion_classification.iloc[i] = 2
        
    elif depletion_counts_pct.iloc[i,0] > 0.2 and depletion_counts_pct.iloc[i,3] > 0.2 \
        and depletion_counts_pct.iloc[i,2] < 0.2:
        
        depletion_classification.iloc[i] = 3

    elif depletion_counts_pct.iloc[i,0] > 0.2 and depletion_counts_pct.iloc[i,2] > 0.2 \
        and depletion_counts_pct.iloc[i,3] < 0.2:
            
        depletion_classification.iloc[i] = 4
    
    #V2
    # elif depletion_counts_pct.iloc[i,1] > 0.2:
    #
    #     depletion_classification.iloc[i] = 5
        
    else:
        continue

# Add NLDAS column to relate farm_ids to NLDAS_ids to enable joining with aquifer_id data which is indexed by NLDAS_id
depletion_classification['FARM_ID'] = depletion_classification.index.astype('int')
depletion_classification = depletion_classification.merge(farm_id_aquifer_id, how = 'inner', left_on = 'FARM_ID', right_on = 'farm_id')

depletion_class_byaquifer = depletion_classification[['Aquifer','depletion_classification']].groupby(['Aquifer']).agg(lambda x:x.value_counts().index[0]).reset_index()


# Save as CSV for mapping and analysis 
#depletion_classification.to_csv("depletion_classification_V4.csv")
# depletion_classification.to_csv("4GIS_depletion_classification.csv")
depletion_classification.to_csv("4GIS_depletion_classification_gammas1_thres50.csv")
depletion_class_byaquifer.to_csv("4GIS_depletion_byaquifer_gamma1.csv")

# Additional analysis, such as aggregating by aquifer or broadersys 
aquifer_subunits = depletion_classification.Aquifer.unique() # 'aquifer' designation 
aquifer_regional = depletion_classification.BroaderSys.unique() # 'broadersys' designation 
aquifer_regional = np.delete(aquifer_regional, (0), axis = 0) # delete '-' entry which = no broadersys designation

#......Add scrip to iterate over subunit or regional aquifer IDs to generate aquifer-level metrics of grid cell classifications 

#%% Depletion Response surface for single grid cell (for Jim AGU)

farm_id = 36335 # IDs = 15557(x200y128), 36335(x32y103)


scenario_results = scenario_depletion.where(scenario_depletion.farm == farm_id) #    
scenario_results = scenario_results.where(scenario_results.gamma_scenario == 1) #    
scenario_results = scenario_results.where(scenario_results.K_scenario == K)
scenario_results = scenario_results.dropna() 

surface = scenario_results.pivot(index="econ_ratio", columns="hydro_ratio", 
                                values="percent_depleted")

# Surface with econ and hydro as axes 
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

# Make data
X = np.linspace(0.5, 1.5, 5, endpoint = True )
Y = np.linspace(0.5, 1.5, 5, endpoint = True )
X, Y = np.meshgrid(X, Y)
Z = surface.copy()

surf = ax.plot_surface(X, Y, Z, cmap= "Blues",
                        linewidth=0, antialiased=False)
ax.set_zlabel('Fraction GW depleted')
ax.set_xlabel('Hydro')
ax.set_ylabel('Econ')
ax.set_zlim([0,1])
ax.set_title('x32y103')
plt.show()




