pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

if first:
    df_combined = df
    first = False
else:
    df_combined = pd.concat([df_combined, df])

gw_multipliers = []
for farm in cum_gw_sum:
    i = 0
    for index, row in gw_cost_curves.iterrows():
        if farm > row['volume']:
            i = index
        else:
            break
    if gw_cost_curves.loc[i][gw_run] != 0:  # JY temp to deal with zeros in groundwater cost curves (check in with Stephen)
        gw_multiplier = gw_cost_curves.loc[i][gw_run]
    else:
        gw_multiplier = 9999999 # multiply cost by large number
    gw_multipliers.append(gw_multiplier)

import csv
import pandas as pd
df = pd.DataFrame(water_constraints_by_farm, index=[0]).T
df.to_csv('temp.csv')



net_prices_gw_subset_update = net_prices_gw_subset.copy()
data_profit_subset = data_profit[['nldas', 'farm_id']].head(53835)
cumul_gw_sum = cumul_gw_sum.merge(data_profit_subset, how='left', on='nldas')
for farm_temp in farm_ids_subset:
    for value in crop_ids_by_farm[farm_temp]:
        net_prices_gw_subset_update[value] *= cumul_gw_sum.iloc[farm_temp]['gw_multiplier']

import pickle

id_nldas = dict(zip(data_profit_short.index, data_profit_short.nldas))

with open('20221005_land_costs.p', 'wb') as handle:
    pickle.dump(dict_test, handle, protocol=pickle.HIGHEST_PROTOCOL)

# Revise to account for removal of "Fodder_Herb category"
crop_ids_by_farm_new = {}
for i in crop_ids_by_farm:
    crop_ids_by_farm_new[i] = crop_ids_by_farm[i][0:10]
crop_ids_by_farm = crop_ids_by_farm_new

with open('20221005_sobol_indices.p', 'wb') as handle:
    pickle.dump(sobol_indices, handle, protocol=pickle.HIGHEST_PROTOCOL)


with open('dict.csv', 'w') as csv_file:
    writer = csv.writer(csv_file)
    for key, value in alphas_total.items():
       writer.writerow([key, value])

(pd.DataFrame.from_dict(data=net_prices_land, orient='index')
   .to_csv('dict_file_net_prices.csv', header=False))

