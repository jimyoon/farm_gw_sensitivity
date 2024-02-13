##### Load all necessary modules

import pandas as pd
import os

##### Adjust pandas setting for debugging
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

os.chdir('C:\\Users\\yoon644\\OneDrive - PNNL\\Documents\\PyProjects\\farmer_gw_archetypes')

sens_results_all = pd.read_csv("./sensitivity_output/farm_gw_ensemble_20231018_K_low.csv")
metric = 'cumul_gw'  # 'total_profit'
sens_results = sens_results_all[['hydro_ratio','econ_ratio','farm', metric]]

farms = sens_results.farm.unique().tolist()
hydro_ratios = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5]
econ_ratios = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5]

hydro_slope_avg = []
econ_slope_avg = []

farm_status = 0
for farm in farms:
    hydro_slope_sum = 0
    econ_slope_sum = 0
    count = 0
    print(farm_status)
    for hydro in hydro_ratios:
        if hydro == 0.5:
            continue
        else:
            for econ in econ_ratios:
                if econ == 0.5:
                    continue
                else:
                    hydro_slope_sum += sens_results[
                                           (sens_results.farm == farm) & (sens_results.hydro_ratio == hydro) & (
                                                       sens_results.econ_ratio == econ)]['cumul_gw'].values[0] - \
                                       sens_results[
                                           (sens_results.farm == farm) & (sens_results.hydro_ratio == round(hydro - 0.1,1)) & (
                                                       sens_results.econ_ratio == econ)]['cumul_gw'].values[0]
                    econ_slope_sum += sens_results[
                                           (sens_results.farm == farm) & (sens_results.hydro_ratio == hydro) & (
                                                       sens_results.econ_ratio == econ)]['cumul_gw'].values[0] - \
                                       sens_results[
                                           (sens_results.farm == farm) & (sens_results.hydro_ratio == hydro) & (
                                                       sens_results.econ_ratio == round(econ - 0.1,1))]['cumul_gw'].values[0]
                    count += 1
    hydro_slope_avg.append(hydro_slope_sum / count)
    econ_slope_avg.append(econ_slope_sum / count)
    farm_status += 1

########### VERSION 2
farm_status = 0
for farm in farms:
    hydro_slope_sum = 0
    econ_slope_sum = 0
    print(farm_status)
    count = 0
    sens_results_filtered = sens_results[(sens_results.farm == farm)]
    for hydro in hydro_ratios:
        if hydro == 0.5:
            continue
        else:
            sens_results_filtered_2 = sens_results_filtered[(sens_results_filtered.hydro_ratio == hydro) | (sens_results_filtered.hydro_ratio == round(hydro - 0.1,1))]
            for econ in econ_ratios:
                if econ == 0.5:
                    continue
                else:
                    hydro_slope_sum += sens_results_filtered_2[(sens_results_filtered_2.hydro_ratio == hydro) & (
                                                       sens_results_filtered_2.econ_ratio == econ)][metric].values[0] - \
                                       sens_results_filtered_2[(sens_results_filtered_2.hydro_ratio == round(hydro - 0.1,1)) & (
                                                       sens_results_filtered_2.econ_ratio == econ)][metric].values[0]
                    econ_slope_sum += sens_results_filtered_2[(sens_results_filtered_2.hydro_ratio == hydro) & (
                                                       sens_results_filtered_2.econ_ratio == econ)][metric].values[0] - \
                                       sens_results_filtered_2[(sens_results_filtered_2.hydro_ratio == hydro) & (
                                                       sens_results_filtered_2.econ_ratio == round(econ - 0.1,1))][metric].values[0]
                    count += 1
    hydro_slope_avg.append(hydro_slope_sum / count)
    econ_slope_avg.append(econ_slope_sum / count)
    farm_status += 1

farms_slope_df = pd.DataFrame(list(zip(farms, hydro_slope_avg, econ_slope_avg)),
              columns=['farms','hydro_slope_avg', 'econ_slope_avg'])

nldas = pd.read_csv('./data_inputs/nldas_farms_subset.csv')

farms_slope_df = pd.merge(farms_slope_df, nldas, how="left", left_on="farms", right_on="index")

# Calculate hydrologic versus economic "control" ratio (cumulative gw pumping version)

farms_slope_df["control_ratio"] = farms_slope_df["econ_slope_avg"] / farms_slope_df["hydro_slope_avg"]
farms_slope_df.loc[farms_slope_df['hydro_slope_avg'] > -.01, "control_ratio"] = 0
farms_slope_df.loc[farms_slope_df['econ_slope_avg'] < .01, "control_ratio"] = 0

farms_slope_df.to_csv('farms_slope_df_cumulgw_lowK.csv')

farms_slope_df_subset = farms_slope_df[(farms_slope_df.control_ratio != 0)]
farms_slope_df_subset.loc[farms_slope_df_subset['control_ratio'] < 0, "control_ratio"] = farms_slope_df_subset['control_ratio'] * -1.0

farms_slope_df_subset.to_csv('farms_slope_df_cumulgw_lowK_subset.csv')

# Calculate hydrologic versus economic "control" ratio (total ag profit version)

farms_slope_df["control_ratio"] = abs(farms_slope_df["econ_slope_avg"]) / abs(farms_slope_df["hydro_slope_avg"])
farms_slope_df.to_csv('farms_slope_df_profit.csv')

# Calculate change in average groundwater pumping between high and low K
results_lowk = pd.read_csv("./sensitivity_output/farm_gw_ensemble_20231018_K_low.csv")
results_highk = pd.read_csv("./sensitivity_output/farm_gw_ensemble_20230919.csv")

aggregation_functions = {'cumul_gw': 'mean', 'perc_vol_depleted': 'mean', 'total_profit':'mean'}
results_highk_summary = results_highk[['farm', 'cumul_gw','perc_vol_depleted', 'total_profit']].groupby(['farm']).aggregate(
    aggregation_functions)
results_lowk_summary = results_lowk[['farm', 'cumul_gw','perc_vol_depleted', 'total_profit']].groupby(['farm']).aggregate(
    aggregation_functions)
results_highk_summary = results_highk_summary.reset_index()
results_lowk_summary =  results_lowk_summary.reset_index()

results_highk_summary = results_highk_summary.rename(columns={'cumul_gw': 'cumul_gw_highk', 'perc_vol_depleted': 'perc_vol_depleted_highk', 'total_profit': 'total_profit_highk'})
results_lowk_summary = results_lowk_summary.rename(columns={'cumul_gw': 'cumul_gw_lowk', 'perc_vol_depleted': 'perc_vol_depleted_lowk', 'total_profit': 'total_profit_lowk'})

results_merged = pd.merge(results_highk_summary, results_lowk_summary, on='farm')

results_merged['cumul_gw_ratio'] = results_merged['cumul_gw_lowk'] / results_merged['cumul_gw_highk']
results_merged['total_profit_ratio'] = results_merged['total_profit_lowk'] / results_merged['total_profit_highk']
results_merged['percent_vol_depl_diff'] = results_merged['perc_vol_depleted_highk'] - results_merged['perc_vol_depleted_lowk']

nldas = pd.read_csv('./data_inputs/nldas_farms_subset.csv')

results_merged = pd.merge(results_merged, nldas, how="left", left_on="farm", right_on="index")

results_merged.to_csv('perc_vol_depleted_highminuslowk.csv')