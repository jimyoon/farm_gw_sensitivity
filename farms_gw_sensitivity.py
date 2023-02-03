##### Load all necessary modules
from pyomo.environ import *
from pyomo.opt import SolverFactory
import pandas as pd
try:
    import cPickle as pickle
except ImportError:  # python 3.x
    import pickle
import numpy as np
from SALib.sample import saltelli
from SALib.analyze import sobol
import time

##### Load Theis drawdown response module
import Theis_pumping_with_deepening

time_start = time.time()

##### Adjust pandas setting for debugging
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

##### Initiate list of all farm ids across CONUS
farm_ids = range(53835)  # total number of farm agents / nldas IDs

##### Load in master farms table (note: this takes a couple minutes as we are loading a large table over all crops/farms,
# need to figure out how to speed up for ensemble runs, e.g., run once as part of SALib script and remove here)
farms_master = pd.ExcelFile("./data_inputs/PMP_inputs_sensitivity_20220829.xlsx")
farms_master = farms_master.parse("crop_level")

##### Load external pickle files for PMP/agents (for master reference, see excel file)
with open('./data_inputs/20220830_id_nldas.p', 'rb') as fp:  # dictionary that includes nldas id keyed on farm id [0-53834]
    id_nldas = pickle.load(fp)
with open('./data_inputs/20220830_crop_ids_by_farm.p', 'rb') as fp:  # dictionary that includes crop ids [0-538349] keyed on farm id [0-53834]
    crop_ids_by_farm = pickle.load(fp)
with open('./data_inputs/20220830_max_land_constr.p', 'rb') as fp:  # dictionary that includes maximum land constraint keyed on farm id [0-53834]
    max_land_constr = pickle.load(fp)
with open('./data_inputs/20220830_sw_calib_constraints.p', 'rb') as fp:  # dictionary that includes sw constraint (used during calibration) keyed on farm id [0-53834]
    sw_calib_constr = pickle.load(fp)
with open('./data_inputs/20220830_gw_calib_constraints.p', 'rb') as fp:  # dictionary that includes gw constraint (used during calibration) keyed on farm id [0-53834]
    gw_calib_constr = pickle.load(fp)
with open('./data_inputs/20220830_alphas_total.p', 'rb') as fp:  # dictionary that includes PMP land alphas [0-538349] keyed on crop id [0-538349]
    alphas_total = pickle.load(fp)
with open('./data_inputs/20220830_gammas_total.p', 'rb') as fp:  # dictionary that includes PMP land gammmas [0-538349] keyed on crop id [0-538349]
    gammas_total = pickle.load(fp)
with open('./data_inputs/20220830_net_prices_land.p', 'rb') as fp:  # dictionary that includes net land prices [0-538349] keyed on crop id [0-538349]
    net_prices_land = pickle.load(fp)
with open('./data_inputs/20220830_net_prices_sw.p', 'rb') as fp:  # dictionary that includes net sw prices [0-538349] keyed on crop id [0-538349]
    net_prices_sw = pickle.load(fp)
with open('./data_inputs/20220830_net_prices_gw.p', 'rb') as fp:  # dictionary that includes net gw prices [0-538349] keyed on crop id [0-538349]
    net_prices_gw = pickle.load(fp)
with open('./data_inputs/20220830_nirs.p', 'rb') as fp:  # dictionary that includes net irrigation requirements [0-538349] keyed on crop id [0-538349]
    nirs = pickle.load(fp)
with open('./data_inputs/20221005_yields.p', 'rb') as fp:  # dictionary that includes crop yields [0-538349] keyed on crop id [0-538349]
    yields = pickle.load(fp)
with open('./data_inputs/20221005_prices.p', 'rb') as fp:  # dictionary that includes crop prices [0-538349] keyed on crop id [0-538349]
    prices = pickle.load(fp)
with open('./data_inputs/20221005_land_costs.p', 'rb') as fp:  # dictionary that includes crop land costs [0-538349] keyed on crop id [0-538349]
    land_costs = pickle.load(fp)

##### Load external files for GW cost curves
gw_cost_curves = pd.read_csv('./data_inputs/20221005_cost_curves/20221005_gw_cost_curves.csv')
gw_cost_lookup = pd.read_csv('./data_inputs/20221005_cost_curves/20221005_gw_cost_lookup.csv')

##### Select farm/gw run options (these are currently implemented as higher-level options than the sensitivity inputs)
farm_id = [15557]  # Pick the farm agent [0-53834] to use as a basis for the sensitivity analysis
gw_cost_id = 'gw4'  # Pick the gw cost curve to use as a basis for the sensitivity analysis
# gw_area_well_sqm = 750 * 750  # Assumed groundwater area irrigated with a well in square meters
# gw_area_well_acres = gw_area_well_sqm * 0.000247105  # Assumed groundwater area irrigated with a well in acres
no_of_years = 50  # The number of years (annual timesteps) to run simulation
# farms_per_grid = 4412  # Assumed number of farms per NLDAS grid cell (50 km x 50 km) / 140 acres = 4412

##### Select groundwater cost curve generation options
gw_curve_option = 'internal'  # GW cost curve option ("external" - pre-processed and pulled from external file, "internal" - generate on the fly)
S = 0.2
m = 50
# K = 100 define in sensitivity analysis
WL = 5
R = 0
IrrDepth = 10
years = 80

##### For each representative farm, estimate # of groundwater wells (taking initial total groundwater area divided by 750m x 750m)
aggregation_functions = {'area_irrigated_gw': 'sum'}
farm_gw_sum = farms_master[['nldas','area_irrigated_gw']].groupby(['nldas']).aggregate(aggregation_functions)
farm_gw_sum = farm_gw_sum.reset_index()
# farm_gw_sum['no_of_wells'] = farm_gw_sum['area_irrigated_gw'] / (gw_area_well_acres)  # .000247105 conversion from square meters to acres

##### Subset relevant data inputs for the specific farm/gw run option
# subset crop ids
crop_ids_by_farm_subset = {key: crop_ids_by_farm[key] for key in farm_id}
ids_subset = []
for key, list_value in crop_ids_by_farm_subset.items():
    for value in list_value:
        ids_subset.append(value)
ids_subset_sorted = sorted(ids_subset)

# subset various dictionaries;
net_prices_land_subset_og = {key: net_prices_land[key] for key in ids_subset_sorted}
net_prices_sw_subset = {key: net_prices_sw[key] for key in ids_subset_sorted}
net_prices_gw_subset = {key: net_prices_gw[key] for key in ids_subset_sorted}
alphas_total_subset = {key: alphas_total[key] for key in ids_subset_sorted}
gammas_total_subset_og = {key: gammas_total[key] for key in ids_subset_sorted}
nirs_subset_og = {key: nirs[key] for key in ids_subset_sorted}
yields_subset = {key: yields[key] for key in ids_subset_sorted}
prices_subset = {key: prices[key] for key in ids_subset_sorted}
land_costs_subset = {key: land_costs[key] for key in ids_subset_sorted}
max_land_constr_subset = {key: max_land_constr[key] for key in farm_id}
sw_calib_constr_subset_og = {key: sw_calib_constr[key] for key in farm_id}
gw_calib_constr_subset = {key: gw_calib_constr[key] for key in farm_id}

# set land net prices to large negative # for gammas that are zero (so PMP won't produce crops, zero gamma value indicates no observed land use)
for key,value in gammas_total_subset_og.items():
    if value == 0:
        net_prices_land_subset_og[key] = -9999999999


time_load = time.time()

###### AFTER TESTING CONVERT FOLLOWING INTO A FUNCTION FOR INCORPORATION INTO SALIB

def farm_gw_model(sw_mult, price_mult, gammas_mult, K):
# def farm_gw_model(nir_mult, sw_mult, price_mult, gammas_mult, K):  # currently only picking a few variables to test (applied uniformly across all crops)

    alphas_mult = 1  # JY temp to disable alphas in sensitivity analysis
    nir_mult = 1 # JY temp to disable alphas in sensitivity analysis

    # generate gw cost curves if gw_curve_option is "internal"
    if gw_curve_option == 'internal':
        gw_cost_curve_internal = Theis_pumping_with_deepening.Analytical(S, m, K, WL, R, IrrDepth, years)

    # apply sensitivity multipliers
    gammas_total_subset = gammas_total_subset_og.copy()
    for key in gammas_total_subset:
        gammas_total_subset[key] *= gammas_mult

    nirs_subset = nirs_subset_og.copy()
    for key in nirs_subset:
        nirs_subset[key] *= nir_mult

    sw_calib_constr_subset = sw_calib_constr_subset_og.copy()
    for key in sw_calib_constr_subset:
        sw_calib_constr_subset[key] *= sw_mult

    net_prices_land_subset = net_prices_land_subset_og.copy()
    for key in net_prices_land_subset:
        net_prices_land_subset[key] = ((yields_subset[key] * prices_subset[key]) - land_costs_subset[key] - alphas_total_subset[key]) * price_mult

    # for key in net_prices_land_subset:
    #     net_prices_land_subset[key] = ((yields_subset[key] * prices_subset[key]) - land_costs_subset[key] - (alphas_total_subset[key] * alphas_mult))

    # initialize counters/trackers
    first = True
    depletion_first = True
    cumul_gw_sum = 0
    gw_multiplier = 1

    # start simulation time loop
    for t in range(no_of_years):
        print('year' + str(t))
        print('!!JY!! ' + str(nir_mult))

        # initialize pyomo model
        fwm_s = ConcreteModel()
        fwm_s.ids = Set(initialize=ids_subset_sorted)
        if first:
            fwm_s.net_prices_gw = Param(fwm_s.ids, initialize=net_prices_gw_subset, mutable=True)
        else:
            net_prices_gw_subset_update = net_prices_gw_subset.copy()
            net_prices_gw_subset_update.update((x, y * gw_multiplier) for x, y in net_prices_gw_subset_update.items())
            fwm_s.net_prices_gw = Param(fwm_s.ids, initialize=net_prices_gw_subset_update, mutable=True)
        fwm_s.farm_ids = Set(initialize=farm_id)
        fwm_s.crop_ids_by_farm = Set(fwm_s.farm_ids, initialize=crop_ids_by_farm_subset)
        fwm_s.crop_ids_by_farm_and_constraint = Set(fwm_s.farm_ids, initialize=crop_ids_by_farm_subset)
        fwm_s.net_prices_sw = Param(fwm_s.ids, initialize=net_prices_sw_subset, mutable=True)
        fwm_s.net_prices_total = Param(fwm_s.ids, initialize=net_prices_land_subset, mutable=True)
        fwm_s.gammas_total = Param(fwm_s.ids, initialize=gammas_total_subset, mutable=True)
        fwm_s.land_constraints = Param(fwm_s.farm_ids, initialize=max_land_constr_subset, mutable=True)
        fwm_s.sw_constraints = Param(fwm_s.farm_ids, initialize=sw_calib_constr_subset, mutable=True)
        fwm_s.gw_constraints = Param(fwm_s.farm_ids, initialize=gw_calib_constr_subset, mutable=True)
        x_start_values = dict(enumerate([0.0] * 3))
        fwm_s.xs_total = Var(fwm_s.ids, domain=NonNegativeReals, initialize=x_start_values)
        fwm_s.xs_sw = Var(fwm_s.ids, domain=NonNegativeReals, initialize=x_start_values)
        fwm_s.xs_gw = Var(fwm_s.ids, domain=NonNegativeReals, initialize=x_start_values)
        fwm_s.nirs = Param(fwm_s.ids, initialize=nirs_subset, mutable=True)

        # define pyomo objective function
        def obj_fun(fwm_s):
            return 0.00001 * sum(sum((fwm_s.net_prices_total[h] * fwm_s.xs_total[h] - 0.5 * fwm_s.gammas_total[h] *
                                      fwm_s.xs_total[h] * fwm_s.xs_total[h]) for h in fwm_s.crop_ids_by_farm[f]) +
                                 sum((fwm_s.net_prices_sw[i] * fwm_s.xs_sw[i]) for i in fwm_s.crop_ids_by_farm[f]) +
                                 sum((fwm_s.net_prices_gw[g] * fwm_s.xs_gw[g]) for g in fwm_s.crop_ids_by_farm[f]) for f
                                 in fwm_s.farm_ids)

        fwm_s.obj_f = Objective(rule=obj_fun, sense=maximize)

        # define constraints
        def land_constraint(fwm_s, ff):
            return sum(fwm_s.xs_total[i] for i in fwm_s.crop_ids_by_farm_and_constraint[ff]) <= \
                   fwm_s.land_constraints[ff]

        fwm_s.c1 = Constraint(fwm_s.farm_ids, rule=land_constraint)

        def obs_lu_constraint_sum(fwm_s, i):
            return fwm_s.xs_sw[i] + fwm_s.xs_gw[i] == fwm_s.xs_total[i]

        fwm_s.c5 = Constraint(fwm_s.ids, rule=obs_lu_constraint_sum)

        def sw_constraint(fwm_s, ff):
            return sum(fwm_s.xs_sw[i] * fwm_s.nirs[i] for i in fwm_s.crop_ids_by_farm_and_constraint[ff]) <= \
                   fwm_s.sw_constraints[ff]

        fwm_s.c2 = Constraint(fwm_s.farm_ids, rule=sw_constraint)

        def gw_constraint(fwm_s, ff):
            return sum(fwm_s.xs_gw[i] * fwm_s.nirs[i] for i in fwm_s.crop_ids_by_farm_and_constraint[ff]) <= \
                   fwm_s.gw_constraints[ff]

        fwm_s.c6 = Constraint(fwm_s.farm_ids, rule=gw_constraint)

        fwm_s.scaling_factor = Suffix(direction=Suffix.EXPORT)
        fwm_s.scaling_factor[fwm_s.xs_total] = 0.0001
        fwm_s.scaling_factor[fwm_s.xs_sw] = 0.0001
        fwm_s.scaling_factor[fwm_s.land_constraints] = 0.000001
        fwm_s.scaling_factor[fwm_s.sw_constraints] = 0.01
        fwm_s.scaling_factor[fwm_s.gw_constraints] = 0.01
        fwm_s.scaling_factor[fwm_s.c1] = 0.000001
        fwm_s.scaling_factor[fwm_s.c5] = 0.000001
        fwm_s.scaling_factor[fwm_s.c2] = 0.01
        fwm_s.scaling_factor[fwm_s.c6] = 0.01

        # create and run the optimization solver
        # opt = SolverFactory("ipopt", solver_io='nl')
        opt = SolverFactory("gurobi", solver_io='python')
        results = opt.solve(fwm_s, keepfiles=False, tee=True)
        print(results.solver.termination_condition)

        # store main model outputs
        result_xs_sw = dict(fwm_s.xs_sw.get_values())
        result_xs_gw = dict(fwm_s.xs_gw.get_values())
        result_xs_total = dict(fwm_s.xs_total.get_values())

        farms_list = [farms_master.loc[[farm_id[0]]].nldas.values[0]]

        # process groundwater production results and
        results_pd = farms_master[farms_master['nldas'].isin(farms_list)]
        results_pd['xs_gw'] = 0
        results_pd['xs_sw'] = 0
        results_pd['xs_total'] = 0
        results_pd['id'] = results_pd['index']
        results_xs_sw_pd = pd.DataFrame.from_dict(result_xs_sw, orient='index')
        results_xs_sw_pd['id'] = results_xs_sw_pd.index + 1
        results_xs_sw_pd = results_xs_sw_pd.rename(columns={0: "xs_sw_temp"})
        results_pd = results_pd.merge(results_xs_sw_pd[['id', 'xs_sw_temp']], how='left', on=['id'])
        results_pd.loc[results_pd['xs_sw_temp'].notnull(), 'xs_sw'] = results_pd['xs_sw_temp']
        results_xs_gw_pd = pd.DataFrame.from_dict(result_xs_gw, orient='index')
        results_xs_gw_pd['id'] = results_xs_gw_pd.index + 1
        results_xs_gw_pd = results_xs_gw_pd.rename(columns={0: "xs_gw_temp"})
        results_pd = results_pd.merge(results_xs_gw_pd[['id', 'xs_gw_temp']], how='left', on=['id'])
        results_pd.loc[results_pd['xs_gw_temp'].notnull(), 'xs_gw'] = results_pd['xs_gw_temp']
        results_xs_total_pd = pd.DataFrame.from_dict(result_xs_total, orient='index')
        results_xs_total_pd['id'] = results_xs_total_pd.index + 1
        results_xs_total_pd = results_xs_total_pd.rename(columns={0: "xs_total_temp"})
        results_pd = results_pd.merge(results_xs_total_pd[['id', 'xs_total_temp']], how='left', on=['id'])
        results_pd.loc[results_pd['xs_total_temp'].notnull(), 'xs_total'] = results_pd['xs_total_temp']
        results_pd = results_pd.drop(['xs_gw_temp', 'xs_sw_temp', 'xs_total_temp'], axis=1)
        results_pd['year'] = t+1
        results_pd['gw_vol'] = results_pd['xs_gw'] * results_pd['nir_corrected']
        try:
            results_pd['gw_mult'] = gw_multiplier
        except NameError:
            results_pd['gw_mult'] = 1
        try:
            results_pd['gw_cumul_vol'] = cumul_gw_sum
        except NameError:
            results_pd['gw_cumul_vol'] = 0
        results_pd['gw_run'] = gw_cost_id
        results_pd['trans'] = gw_cost_lookup[(gw_cost_lookup.row_indx == gw_cost_id)]['Trans'].values[0]
        if first:
            results_combined = results_pd
        else:
            results_combined = pd.concat([results_combined, results_pd])
        aggregation_functions = {'gw_vol': 'sum', 'xs_gw': 'sum'}
        gw_sum = results_pd[['nldas','gw_vol','xs_gw']].groupby(['nldas']).aggregate(aggregation_functions)
        gw_sum = gw_sum.reset_index()
        # gw_sum = gw_sum.merge(farm_gw_sum[['nldas', 'no_of_wells']], how='left', on=['nldas'])
        # gw_sum['gw_vol_well'] = gw_sum['gw_vol'] * farms_per_grid / gw_sum['no_of_wells']
        gw_sum['gw_vol_km3'] = gw_sum['gw_vol'] * 1.23348e-6
        if first:
            cumul_gw_sum = gw_sum['gw_vol_km3'].values[0]
        else:
            cumul_gw_sum += gw_sum['gw_vol_km3'].values[0]
        i = 0

        if gw_curve_option == 'external':
            for index, row in gw_cost_curves.iterrows():
                if cumul_gw_sum > row['volume']:
                    i = index
                else:
                    break
            if gw_cost_curves.loc[i][gw_cost_id] != 0:  # JY temp to deal with zeros in groundwater cost curves (check in with Stephen)
                gw_multiplier = gw_cost_curves.loc[i][gw_cost_id]
            else:
                gw_multiplier = 9999999999999  # Set groundwater cost to extremely high value to reflect groundwater exhaustion

        elif gw_curve_option == 'internal':
            for index in range(gw_cost_curve_internal[0][0].size):
                if cumul_gw_sum > gw_cost_curve_internal[1][0][index]:
                    i = index
                else:
                    break
            if gw_cost_curve_internal[0][0][i] != 0:
                gw_multiplier = gw_cost_curve_internal[0][0][i] / gw_cost_curve_internal[0][0][0]
            else:
                gw_multiplier = 9999999999999
                if depletion_first:
                    time_depletion = t
                    depletion_first = False

        print(cumul_gw_sum)
        first = False

    # Calculate
    results_combined['profit_calc'] = ((price_mult * (results_combined['price'] * results_combined['yield']) -
                                             results_combined['land_cost'] - results_combined['alphas_land']) * results_combined['xs_total']) \
                                           - (0.5 * results_combined['gammas_total'] * results_combined['xs_total'] * results_combined['xs_total']) \
                                           - (results_combined['net_prices_gw'] * results_combined['xs_gw']) \
                                           - (results_combined['net_prices_sw'] * results_combined['xs_sw'])


    #  Calculate summary results (percent change in total irrigated area as an initial result)
    aggregation_functions = {'xs_total': 'sum', 'profit_calc': 'sum'}
    summary_pd = results_combined[['nldas','year','xs_total','profit_calc']].groupby(['nldas','year']).aggregate(aggregation_functions)
    summary_pd = summary_pd.reset_index()
    year_start = summary_pd.year.min()
    year_end = summary_pd.year.max()
    end_div_start_area = summary_pd[(summary_pd.year==year_end)].xs_total.values[0] / summary_pd[(summary_pd.year == year_start)].xs_total.values[0]
    total_profit = summary_pd.profit_calc.sum()

    # Volume depleted (calculated from cost curves)
    perc_vol_depleted = cumul_gw_sum / max(gw_cost_curve_internal[1][0].tolist())
    if perc_vol_depleted > 1:
        perc_vol_depleted = 1


    return [end_div_start_area, total_profit, perc_vol_depleted, time_depletion]

##### Run sensitivity analysis using SALib
#
# problem = {
#     'num_vars': 5,
#     'names': ['nir_mult','sw_mult','price_mult','alphas_mult','gammas_mult'],
#     'bounds': [[0.8, 1.2], [0.8, 1.2], [0.8, 1.2], [0.8, 1.2], [0.8, 1.2]]
# }
#
#
# problem = {
#     'num_vars': 5,
#     'names': ['nir_mult','sw_mult', 'price_mult','K','gammas_mult'],
#     'bounds': [[0.8, 1.2], [0.8, 1.2], [0.8, 1.2], [1, 100], [0.8, 1.2]]
# }

problem = {
    'num_vars': 4,
    'names': ['sw_mult','price_mult','gammas_mult','K'],
    'bounds': [[0.8, 1.2], [0.8, 1.2], [0.8, 1.2], [1, 100]]
}

# param_values = saltelli.sample(problem, 1024)
# param_values = saltelli.sample(problem, 256)
# param_values = saltelli.sample(problem, 256)
param_values = saltelli.sample(problem, 64)

# y = np.array([farm_gw_model(*params) for params in param_values])

count = 0
results_end_div_start_area = []
results_total_profit = []
results_perc_vol_depleted = []
results_time_depletion = []
for params in param_values:
    results = farm_gw_model(*params)
    results_end_div_start_area.append(results[0])
    results_total_profit.append(results[1])
    results_perc_vol_depleted.append(results[2])
    results_time_depletion.append(results[3])
    count += 1

time_ensemble = time.time()

results_array = np.array(results_total_profit)
sobol_indices = sobol.analyze(problem, results_array)