from pyomo.environ import *
from pyomo.opt import SolverFactory
import os
import pandas as pd
import numpy as np
try:
    import cPickle as pickle
except ImportError:  # python 3.x
    import pickle
import shutil
import netCDF4
import logging
import sys


water_constraints_by_farm = pd.read_csv('./data_inputs/hist_avail_bias_correction_20220223.csv') # Use baseline water demand data for warmup period
water_constraints_by_farm = water_constraints_by_farm[['NLDAS_ID','sw_irrigation_vol']].reset_index()
water_constraints_by_farm = water_constraints_by_farm['sw_irrigation_vol'].to_dict()

gw_cost_curves = pd.read_csv('./data_inputs/gw_cost_curves.csv')
gw_cost_lookup = pd.read_csv('./data_inputs/Test_cost_curves_edit.csv')


## Read in PMP calibration files
data_file=pd.ExcelFile("./data_inputs/MOSART_WM_PMP_inputs_20220323_GW.xlsx")
data_profit = data_file.parse("Profit")
data_profit['farm_id'] = data_profit.index
water_nirs=data_profit["nir_corrected"]
nirs=dict(water_nirs)

# For each representative farm, estimate # of groundwater wells (taking initial total groundwater area divided by 750m x 750m)
aggregation_functions = {'area_irrigated_gw': 'sum'}
farm_gw_sum = data_profit[['nldas','area_irrigated_gw']].groupby(['nldas']).aggregate(aggregation_functions)
farm_gw_sum = farm_gw_sum.reset_index()
farm_gw_sum['no_of_wells'] = farm_gw_sum['area_irrigated_gw'] / (750.0 * 750.0 * 0.000247105)  # .000247105 conversion from square meters to acres

## C.1. Preparing model indices and constraints:
ids = range(538350) # total number of crop and nldas ID combinations
farm_ids = range(53835) # total number of farm agents / nldas IDs
with open('./data_inputs/crop_ids_by_farm.p', 'rb') as fp:
    crop_ids_by_farm = pickle.load(fp)
with open('./data_inputs/crop_ids_by_farm_and_constraint.p', 'rb') as fp:
    crop_ids_by_farm_and_constraint = pickle.load(fp)
with open('./data_inputs/max_land_constr_20220307_protocol2.p', 'rb') as fp:
    land_constraints_by_farm = pickle.load(fp)

# Revise to account for removal of "Fodder_Herb category"
crop_ids_by_farm_new = {}
for i in crop_ids_by_farm:
    crop_ids_by_farm_new[i] = crop_ids_by_farm[i][0:10]
crop_ids_by_farm = crop_ids_by_farm_new
crop_ids_by_farm_and_constraint = crop_ids_by_farm_new

# Load gammas, net prices, etc
with open('./data_inputs/gammas_total_20220323_protocol2.p', 'rb') as fp:
    gammas_total = pickle.load(fp)
with open('./data_inputs/net_prices_total_20220323_protocol2.p', 'rb') as fp:
    net_prices_total = pickle.load(fp)
with open('./data_inputs/alphas_total_20220323_protocol2.p', 'rb') as fp:
    alphas_total = pickle.load(fp)
with open('./data_inputs/net_prices_sw_20220323_protocol2.p', 'rb') as fp:
    net_prices_sw = pickle.load(fp)
with open('./data_inputs/alphas_sw_20220323_protocol2.p', 'rb') as fp:
    alphas_sw = pickle.load(fp)
with open('./data_inputs/gammas_sw_20220323_protocol2.p', 'rb') as fp:
    gammas_sw = pickle.load(fp)
with open('./data_inputs/net_prices_gw_20220323_protocol2.p', 'rb') as fp:
    net_prices_gw = pickle.load(fp)

x_start_values=dict(enumerate([0.0]*3))

# subset various dictionaries;
# keys_to_extract = [36335]  # will multiply start and end values by n+1 once integrated in loop
# keys_to_extract = [15557]
# keys_to_extract = [42988]

# gw_list = ['gw7', 'gw6', 'gw14', 'gw13', 'gw12', 'gw3', 'gw11', 'gw10', 'gw9', 'gw17', 'gw16', 'gw15']
gw_list = ['gw7', 'gw6', 'gw14', 'gw13', 'gw12', 'gw3', 'gw11', 'gw10', 'gw9', 'gw17', 'gw16', 'gw15']
first_a = True
for gw_run in gw_list:
    first_b = True
    year = 1
    cumul_gw_sum = 0
    gw_multiplier = 1
    for t in range(5):

        chunk_size = 555
        no_of_chunks = len(farm_ids) / chunk_size

        for n in range(65,66):
        # for n in range(int(no_of_chunks)):
            print(gw_run + ' ///// chunk: ' + str(n) + ' ///// timestep: ' + str(t))
            farm_ids_subset = list(range(chunk_size*n, chunk_size*(n+1)))
            crop_ids_by_farm_subset = {key: crop_ids_by_farm[key] for key in farm_ids_subset}
            ids_subset = []
            for key,list_value in crop_ids_by_farm_subset.items():
                for value in list_value:
                    ids_subset.append(value)
            ids_subset_sorted = sorted(ids_subset)

            keys_to_extract = list(range(chunk_size*n, chunk_size*(n+1)))
            net_prices_total_subset = {key: net_prices_total[key] for key in ids_subset_sorted}
            net_prices_sw_subset = {key: net_prices_sw[key] for key in ids_subset_sorted}
            net_prices_gw_subset = {key: net_prices_gw[key] for key in ids_subset_sorted}
            gammas_total_subset = {key: gammas_total[key] for key in ids_subset_sorted}
            nirs_subset = {key: nirs[key] for key in ids_subset_sorted}
            alphas_sw_subset = {key: alphas_sw[key] for key in farm_ids_subset}
            gammas_sw_subset = {key: gammas_sw[key] for key in farm_ids_subset}
            land_constraints_by_farm_subset = {key: land_constraints_by_farm[key] for key in farm_ids_subset}
            water_constraints_by_farm_subset = {key: water_constraints_by_farm[key] for key in farm_ids_subset}

            # set price to zero for gammas that are zero
            for key,value in gammas_total_subset.items():
                if value == 0:
                    net_prices_total_subset[key] = -9999999999

            fwm_s = ConcreteModel()
            fwm_s.ids = Set(initialize=ids_subset_sorted)
            if first_b:
                fwm_s.net_prices_gw = Param(fwm_s.ids, initialize=net_prices_gw_subset, mutable=True)
            else:
                net_prices_gw_subset_update = net_prices_gw_subset.copy()
                if 'farm_id' not in cumul_gw_sum.columns:
                    data_profit_subset = data_profit[['nldas', 'farm_id']].head(53835)
                    cumul_gw_sum = cumul_gw_sum.merge(data_profit_subset, how='left', on='nldas')
                for farm_temp in farm_ids_subset:
                    for value in crop_ids_by_farm[farm_temp]:
                        net_prices_gw_subset_update[value] *= cumul_gw_sum[(cumul_gw_sum.farm_id==farm_temp)]['gw_multiplier'].values[0]
                # net_prices_gw_subset_update.update((x, y * gw_multiplier) for x, y in net_prices_gw_subset_update.items())
                fwm_s.net_prices_gw = Param(fwm_s.ids, initialize=net_prices_gw_subset_update, mutable=True)
            fwm_s.farm_ids = Set(initialize=farm_ids_subset)
            fwm_s.crop_ids_by_farm = Set(fwm_s.farm_ids, initialize=crop_ids_by_farm_subset)
            fwm_s.crop_ids_by_farm_and_constraint = Set(fwm_s.farm_ids, initialize=crop_ids_by_farm_subset)
            fwm_s.net_prices_sw = Param(fwm_s.ids, initialize=net_prices_sw_subset, mutable=True)
            fwm_s.net_prices_total = Param(fwm_s.ids, initialize=net_prices_total_subset, mutable=True)
            fwm_s.gammas_total = Param(fwm_s.ids, initialize=gammas_total_subset, mutable=True)
            fwm_s.alphas_sw = Param(fwm_s.farm_ids, initialize=alphas_sw_subset, mutable=True)
            fwm_s.gammas_sw = Param(fwm_s.farm_ids, initialize=gammas_sw_subset, mutable=True)
            fwm_s.land_constraints = Param(fwm_s.farm_ids, initialize=land_constraints_by_farm_subset, mutable=True)
            fwm_s.water_constraints = Param(fwm_s.farm_ids, initialize=water_constraints_by_farm_subset,
                                            mutable=True)
            fwm_s.xs_total = Var(fwm_s.ids, domain=NonNegativeReals, initialize=x_start_values)
            fwm_s.xs_sw = Var(fwm_s.ids, domain=NonNegativeReals, initialize=x_start_values)
            fwm_s.xs_gw = Var(fwm_s.ids, domain=NonNegativeReals, initialize=x_start_values)
            # obs_lu_total = dict(data_profit["area_irrigated"])
            # obs_lu_sw = dict(area_irrigated_sw_farm["area_irrigated_sw"])
            # fwm_s.obs_lu_total = Param(fwm_s.ids, initialize=obs_lu_total, mutable=True)
            # fwm_s.obs_lu_sw = Param(fwm_s.ids, initialize=obs_lu_sw, mutable=True)
            fwm_s.nirs = Param(fwm_s.ids, initialize=nirs_subset, mutable=True)

            ## C.2. Constructing model functions:
            def obj_fun(fwm_s):
                return 0.00001 * sum(sum((fwm_s.net_prices_total[h] * fwm_s.xs_total[h] - 0.5 * fwm_s.gammas_total[
                    h] * fwm_s.xs_total[h] * fwm_s.xs_total[h]) for h in fwm_s.crop_ids_by_farm[f]) +
                                     sum((fwm_s.net_prices_sw[i] * fwm_s.xs_sw[i]) for i in
                                         fwm_s.crop_ids_by_farm[f]) +
                                     sum((fwm_s.net_prices_gw[g] * fwm_s.xs_gw[g]) for g in
                                         fwm_s.crop_ids_by_farm[f]) -
                                     (fwm_s.alphas_sw[f] * sum(fwm_s.xs_sw[s] for s in fwm_s.crop_ids_by_farm[f])) -
                                     (0.5 * fwm_s.gammas_sw[f] * sum(
                                         fwm_s.xs_sw[t] for t in fwm_s.crop_ids_by_farm[f])) * sum(
                    fwm_s.xs_sw[u] for u in fwm_s.crop_ids_by_farm[f]) for f in
                                     fwm_s.farm_ids)  # JY double check this!

            fwm_s.obj_f = Objective(rule=obj_fun, sense=maximize)

            def land_constraint(fwm_s, ff):
                return sum(fwm_s.xs_total[i] for i in fwm_s.crop_ids_by_farm_and_constraint[ff]) <= \
                       fwm_s.land_constraints[ff]

            fwm_s.c1 = Constraint(fwm_s.farm_ids, rule=land_constraint)

            def obs_lu_constraint_sum(fwm_s, i):
                return fwm_s.xs_sw[i] + fwm_s.xs_gw[i] == fwm_s.xs_total[i]

            fwm_s.c5 = Constraint(fwm_s.ids, rule=obs_lu_constraint_sum)
            #
            # def water_constraint(fwm_s, ff):
            #     return sum(fwm_s.xs_sw[i]*fwm_s.nirs[i]*1000 for i in fwm_s.crop_ids_by_farm_and_constraint[ff]) <= fwm_s.water_constraints[ff]  # JY multiplication by 1,000 gets us back to non-scaled NIR values
            # fwm_s.c2 = Constraint(fwm_s.farm_ids, rule=water_constraint)


            opt = SolverFactory("ipopt", solver_io='nl')
            results = opt.solve(fwm_s, keepfiles=False, tee=True)
            print(results.solver.termination_condition)

            result_xs_sw = dict(fwm_s.xs_sw.get_values())
            result_xs_gw = dict(fwm_s.xs_gw.get_values())
            result_xs_total = dict(fwm_s.xs_total.get_values())

            # farms_list = ['x32y103']
            # farms_list = ['x200y128']
            # farms_list = ['x38y97']


            results_pd = data_profit
            results_pd['xs_gw'] = 0
            results_pd['xs_sw'] = 0
            results_pd['xs_total'] = 0
            results_pd['id'] = results_pd['index']
            results_xs_sw_pd = pd.DataFrame.from_dict(result_xs_sw, orient='index')
            results_xs_sw_pd['id'] = results_xs_sw_pd.index + 1
            results_xs_sw_pd = results_xs_sw_pd.rename(columns={0: "xs_sw_temp"})
            results_pd = results_pd.merge(results_xs_sw_pd[['id', 'xs_sw_temp']], how='left', on=['id'])

            results_pd = results_pd[(results_pd['xs_sw_temp'].notnull())]

            # results_pd.loc[results_pd['xs_sw_temp'].notnull(), 'xs_sw'] = results_pd['xs_sw_temp']

            results_pd['xs_sw'] = results_pd['xs_sw_temp']

            results_xs_gw_pd = pd.DataFrame.from_dict(result_xs_gw, orient='index')
            results_xs_gw_pd['id'] = results_xs_gw_pd.index + 1
            results_xs_gw_pd = results_xs_gw_pd.rename(columns={0: "xs_gw_temp"})
            results_pd = results_pd.merge(results_xs_gw_pd[['id', 'xs_gw_temp']], how='left', on=['id'])

            # results_pd.loc[results_pd['xs_gw_temp'].notnull(), 'xs_gw'] = results_pd['xs_gw_temp']

            results_pd['xs_gw'] = results_pd['xs_gw_temp']

            results_xs_total_pd = pd.DataFrame.from_dict(result_xs_total, orient='index')
            results_xs_total_pd['id'] = results_xs_total_pd.index + 1
            results_xs_total_pd = results_xs_total_pd.rename(columns={0: "xs_total_temp"})
            results_pd = results_pd.merge(results_xs_total_pd[['id', 'xs_total_temp']], how='left', on=['id'])

            #results_pd.loc[results_pd['xs_total_temp'].notnull(), 'xs_total'] = results_pd['xs_total_temp']

            results_pd['xs_total'] = results_pd['xs_total_temp']

            results_pd = results_pd.drop(['xs_gw_temp', 'xs_sw_temp', 'xs_total_temp'], axis=1)
            results_pd['year'] = year
            results_pd['gw_vol'] = results_pd['xs_gw'] * results_pd['nir_corrected']

            results_pd['gw_run'] = gw_run
            results_pd['trans'] = gw_cost_lookup[(gw_cost_lookup.row_indx == gw_run)]['Trans'].values[0]

            aggregation_functions = {'gw_vol': 'sum', 'xs_gw': 'sum'}
            gw_sum = results_pd[['nldas','gw_vol','xs_gw']].groupby(['nldas']).aggregate(aggregation_functions)
            gw_sum = gw_sum.reset_index()
            gw_sum = gw_sum.merge(farm_gw_sum[['nldas', 'no_of_wells']], how='left', on=['nldas'])
            farms_per_grid = 4412 # (50 km x 50 km) / 140 acres = 4412
            gw_sum['gw_vol_well'] = gw_sum['gw_vol'] * farms_per_grid / gw_sum['no_of_wells']
            gw_sum['gw_vol_well_km3'] = gw_sum['gw_vol_well'] * 1.23348e-6
            gw_sum['gw_vol_well'] = gw_sum['gw_vol_well'].fillna(0)
            gw_sum['gw_vol_well_km3'] = gw_sum['gw_vol_well_km3'].fillna(0)
            if first_b:
                # cum_gw_sum = gw_sum['gw_vol_well_km3'].values[0]
                cumul_gw_sum = gw_sum[['nldas','gw_vol_well_km3']]  # JY PD series rather than single value
                cumul_gw_sum['gw_multiplier'] = 0.0
            else:
                # cum_gw_sum += gw_sum['gw_vol_well_km3'].values[0]
                cumul_gw_sum['gw_vol_well_km3'] += gw_sum['gw_vol_well_km3']  # JY PD series rather than single value

            results_pd = results_pd.merge(cumul_gw_sum[['nldas', 'gw_multiplier']], how='left', on='nldas')
            results_pd = results_pd.merge(cumul_gw_sum[['nldas', 'gw_vol_well_km3']], how='left', on='nldas')

            if first_a:
                results_combined = results_pd
                first_a = False
            else:
                results_combined = pd.concat([results_combined, results_pd])

            gw_multipliers = []

            for index1, farm in cumul_gw_sum.iterrows():
                i = 0
                for index2, row in gw_cost_curves.iterrows():
                    if farm['gw_vol_well_km3'] > row['volume']:
                        i = index2
                    else:
                        break
                if gw_cost_curves.loc[i][gw_run] != 0:  # JY temp to deal with zeros in groundwater cost curves (check in with Stephen)
                    gw_multiplier = gw_cost_curves.loc[i][gw_run]
                else:
                    gw_multiplier = 2 # multiply cost by large number
                cumul_gw_sum.at[index1,'gw_multiplier'] = gw_multiplier

            # i = 0
            # for index, row in gw_cost_curves.iterrows():
            #     if cum_gw_sum > row['volume']:
            #         i = index
            #     else:
            #         break
            # if gw_cost_curves.loc[i][gw_run] != 0:  # JY temp to deal with zeros in groundwater cost curves (check in with Stephen)
            #     gw_multiplier = gw_cost_curves.loc[i][gw_run]

            year += 1
            first_b = False


# Post-processing for GIS visualizations
year_start = results_combined['year'].min()
year_end = results_combined['year'].max()
results_start = results_combined[(results_combined.year == year_start)]
results_end = results_combined[(results_combined.year == year_end)]
