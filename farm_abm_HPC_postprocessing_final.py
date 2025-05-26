# -*- coding: utf-8 -*-

import pandas as pd
import os 

# List of farm_ids to iterate over 
os.chdir("PATH FOR HPC OUTPUT PROCESSING")
farm_ids = pd.read_csv('nldas_farms_subset_final.csv', index_col = 0) # index is numeric farm_id, data column is NLDAS grid id
farm_ids['farm_id'] = farm_ids.index # add column with farm_id 
farm_ids = farm_ids.reset_index(drop = True) # reset index 

# To keep file sizes manageable, create multiple aggregated dateframes that
# contain ALL results for a given attribute (K, WL, etc.) or 
# metric (GW_Area, Profit, etc.)

# Make empty dataframes
scenario_attributes = pd.DataFrame(columns = ['hydro_ratio', 'econ_ratio', 'K_scenario', 'gamma_scenario', 'farm'])
depletion = pd.DataFrame(columns = ['perc_vol_depleted'])
m_flag = pd.DataFrame(columns = ['m'])

for i in range(len(farm_ids.iloc[:])):
    print(i)
    
    # try to import HPC output for farm_id
    try:
        # Import cases summary csv for individual farm 
        output_file = pd.read_csv("farm_" + str(int(farm_ids.farm_id[i])) + "_cases.csv")
        
        # Concatenate individual farm outputs to Dataframes 
        scenario_attributes = pd.concat([scenario_attributes, output_file.iloc[:,0:5]])
        
        depletion = pd.concat([depletion, output_file.perc_vol_depleted])

        m_flag = pd.concat([m_flag, output_file.m])

    # pass if no result output file for farm_id 
    except FileNotFoundError:
        pass 

## Save aggregated files - all files are in the "HPC_outputs" folder but we 
## processed in 4 batches (first 10,000 indexes of farm_ids, next 10,000, third 10,000, and remaining)
## to speed up processing due to memory limitations of laptop
# scenario_attributes.to_csv("scenario_attributes_batch_4.csv")
# depletion.to_csv("depletion_scenario_outputs_batch_4.csv")
# m_flag.to_csv("m_flag_scenario_outputs_batch_4.csv")
