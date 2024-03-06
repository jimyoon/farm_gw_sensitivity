##### Load all necessary modules

import pandas as pd
import os

##### Adjust pandas setting for debugging
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

os.chdir('C:\\Users\\yoon644\\OneDrive - PNNL\\Documents\\PyProjects\\farmer_gw_archetypes')

results = pd.read_csv("./farms_slope_df_profit.csv")

##### Load all relevant inputs into master table

farms_master_excel = pd.ExcelFile("./data_inputs/PMP_inputs_sensitivity_20220829.xlsx")
farms_master = farms_master_excel.parse("crop_level")
farms_id = farms_master_excel.parse("farm_level")
farms_id['farm'] = farms_id.index

farms_master_param = farms_master[['nldas','crop','gammas_total']]
farms_master_param = farms_master_param.set_index(['nldas', 'crop']).gammas_total.unstack(fill_value='').reset_index()
new_names = [(i,'gammas_' + i) for i in farms_master_param.iloc[:, 1:].columns.values]
farms_master_param.rename(columns = dict(new_names), inplace=True)

farms_master_add = farms_master[['nldas','crop','alphas_land']]
farms_master_add = farms_master_add.set_index(['nldas', 'crop']).alphas_land.unstack(fill_value='').reset_index()
new_names = [(i,'alphas_' + i) for i in farms_master_add.iloc[:, 1:].columns.values]
farms_master_add.rename(columns = dict(new_names), inplace=True)

farms_master_param = pd.merge(farms_master_param, farms_master_add, on='nldas')

nldas_gw_attributes = pd.read_csv('./NLDAS_Cost_Curve_Attributes.csv')

farms_master_param = pd.merge(farms_master_param, nldas_gw_attributes[['NLDAS_ID','K','K_hi','K_de_graaf']], left_on = 'nldas', right_on = 'NLDAS_ID')

farms_master_param = pd.merge(farms_master_param, farms_id[['nldas','farm']], on='nldas')

farms_master_param = farms_master_param.drop('NLDAS_ID', axis=1)
cols = farms_master_param.columns.tolist()

cols = ['farm',
    'nldas',
 'gammas_Corn',
 'gammas_FiberCrop',
 'gammas_FodderGrass',
 'gammas_MiscCrop',
 'gammas_OilCrop',
 'gammas_OtherGrain',
 'gammas_Rice',
 'gammas_Root_Tuber',
 'gammas_SugarCrop',
 'gammas_Wheat',
 'alphas_Corn',
 'alphas_FiberCrop',
 'alphas_FodderGrass',
 'alphas_MiscCrop',
 'alphas_OilCrop',
 'alphas_OtherGrain',
 'alphas_Rice',
 'alphas_Root_Tuber',
 'alphas_SugarCrop',
 'alphas_Wheat',
 'K',
 'K_hi',
 'K_de_graaf']

farms_master_param = farms_master_param[cols]
farms_master_param.to_csv('farms_param_inputs.csv')
