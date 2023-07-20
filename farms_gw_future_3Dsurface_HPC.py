##### Load all necessary modules #####
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
import os
import matplotlib.pyplot as plt
from itertools import product

##### Travis: Lines 19-92 all load external data and set variables that are used for all runs -- need to find strategy to do this only once when deployed on HPC
##### Travis: Need to modify this such that it can be deployed on HPC
os.chdir('C:\\Users\\yoon644\\OneDrive - PNNL\\Documents\\PyProjects\\farmer_gw_archetypes')

##### Load Theis drawdown response module
# import Theis_pumping_with_deepening_extended_extraction
import Superwell_for_ABM_on_the_fly

time_start = time.time()

##### Adjust pandas setting for debugging
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

##### Initiate list of all farm ids across CONUS
farm_ids = range(53835)  # total number of farm agents / nldas IDs

##### Load in master farms table (note: this takes a couple minutes as we are loading a large table over all crops/farms,
##### Travis: This is a large file that loads the data for all farms across the CONUS -- for ensemble experiment, we'd only want to load it once.
farms_master = pd.ExcelFile("./data_inputs/PMP_inputs_sensitivity_20220829.xlsx")
farms_master = farms_master.parse("crop_level")

##### Load external pickle files for PMP/agents (for master reference, see excel file)
##### Travis: Likewise, for the pickle files below. Though these load much faster, we'd still only want to load them once as they each contain information across all farms.
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
# with open('./data_inputs/20221005_yields.p', 'rb') as fp:  # dictionary that includes crop yields [0-538349] keyed on crop id [0-538349]
#     yields = pickle.load(fp)
# with open('./data_inputs/20221005_prices.p', 'rb') as fp:  # dictionary that includes crop prices [0-538349] keyed on crop id [0-538349]
#     prices = pickle.load(fp)
# with open('./data_inputs/20221005_land_costs.p', 'rb') as fp:  # dictionary that includes crop land costs [0-538349] keyed on crop id [0-538349]
#     land_costs = pickle.load(fp)

keys = np.arange(538350)
yields = dict(zip(keys,farms_master.iloc[:,6]))
prices = dict(zip(keys,farms_master.price[:]))
land_costs = dict(zip(keys,farms_master.land_cost[:]))

##### For each representative farm, estimate # of groundwater wells (taking initial total groundwater area divided by 750m x 750m)
farms_master['gw_irr_vol'] = farms_master['area_irrigated_gw'] * farms_master['nir_corrected'] * 1000
aggregation_functions = {'gw_irr_vol': 'sum'}
farm_gw_sum = farms_master[['nldas', 'gw_irr_vol']].groupby(['nldas']).aggregate(aggregation_functions)
farm_gw_sum = farm_gw_sum.reset_index()
# farm_gw_sum['no_of_wells'] = farm_gw_sum['area_irrigated_gw'] / (gw_area_well_acres)  # .000247105 conversion from square meters to acres

time_load = time.time()

# Universal cost curve inputs
gw_curve_option = 'internal'  # GW cost curve option ("external" - pre-processed and pulled from external file, "internal" - generate on the fly)
WL = 10
S = 0.25
NUM_YEARS = 80
ELECTRICITY_RATE = .12
m = 110
R = .15
gw_cost_curve_internal = None

##### Load external files for NLDAS cost curve attributes
nldas_gw_attributes = pd.read_csv('./NLDAS_Cost_Curve_Attributes.csv')

##### Travis: Lines 96-460 define the method for running the farm-groundwater model
# Function-based workflow for incorporation in various sensitivity/ensemble experiments (e.g., SALib Sobol Analysis)
def farm_gw_model(hydro_ratio, econ_ratio, K, farm_id, m, R, gw_cost_curve_internal):

    # ratios used for future exploratory scenario analysis
    percent_increment = (hydro_ratio - 1.0) / (hydro_ratio + 1.0)
    numer_hydro_factor = 1 + percent_increment
    denom_hydro_factor = 1 - percent_increment
    percent_increment = (econ_ratio - 1.0) / (econ_ratio + 1.0)
    numer_econ_factor = 1 + percent_increment
    denom_econ_factor = 1 - percent_increment

    farm_id = [farm_id]  # Pick the farm agent [0-53834] to use as a basis for the sensitivity analysis
    # farm_id = [12641] # use if running a single farm, enter farm_id
    # gw_cost_id = 'gw4'  # Pick the gw cost curve to use as a basis for the sensitivity analysis
    # gw_area_well_sqm = 750 * 750  # Assumed groundwater area irrigated with a well in square meters
    # gw_area_well_acres = gw_area_well_sqm * 0.000247105  # Assumed groundwater area irrigated with a well in acres
    no_of_years = 100  # The number of years (annual timesteps) to run simulation, can be toggled
    # farms_per_grid = 4412  # Assumed number of farms per NLDAS grid cell (50 km x 50 km) / 140 acres = 4412

    ##### For each representative farm, estimate # of groundwater wells (taking initial total groundwater area divided by 750m x 750m)
    # farms_master['gw_irr_vol'] = farms_master['area_irrigated_gw'] * farms_master['nir_corrected'] * 1000
    # aggregation_functions = {'gw_irr_vol': 'sum'}
    # farm_gw_sum = farms_master[['nldas', 'gw_irr_vol']].groupby(['nldas']).aggregate(aggregation_functions)
    # farm_gw_sum = farm_gw_sum.reset_index()
    # farm_gw_sum['no_of_wells'] = farm_gw_sum['area_irrigated_gw'] / (gw_area_well_acres)  # .000247105 conversion from square meters to acres

    ##### Subset relevant data inputs for the specific farm/gw run option
    # subset crop ids
    crop_ids_by_farm_subset = {key: crop_ids_by_farm[key] for key in
                               farm_id}  # farm_id = row in id_nldas, dictionary with farm_id as key and list of crops as value
    ids_subset = []
    for key, list_value in crop_ids_by_farm_subset.items():  # creates list where first entry is farm_id and each follow entry is a crop index which is
        for value in list_value:  # row in farms_master DataFrame
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
    for key, value in gammas_total_subset_og.items():
        if value == 0:
            net_prices_land_subset_og[key] = -9999999999

    # time_load = time.time()

    # alphas_mult = 1  # JY temp to disable alphas in sensitivity analysis
    # nir_mult = 1  # JY temp to disable alphas in sensitivity analysis

    # generate gw cost curves if more than one cost curve selected in sensitivity sweep
    if gw_curve_option == 'internal':

        nldas_id = farms_master.iloc[farm_id]['nldas']
        R = nldas_gw_attributes[(nldas_gw_attributes['NLDAS_ID'] == nldas_id.values[0])]['Recharge_usgs']
        R = R * numer_hydro_factor
        IRR_DEPTH_afy = farm_gw_sum.iloc[farm_id]['gw_irr_vol']
        # IRR_DEPTH_afy = farms_master.iloc[farm_id]['area_irrigated_gw'] * farms_master.iloc[farm_id]['nir_corrected'] * 1000
        IRR_DEPTH_m3 = IRR_DEPTH_afy * 1233.48  # convert from acre-ft/yr to cubic meters/yr
        IRR_DEPTH = IRR_DEPTH_m3.values[0] / nldas_gw_attributes[(nldas_gw_attributes['NLDAS_ID'] == nldas_id.values[0])]['Area'].values[0] # convert to m/year
        # gw_cost_curve_internal = Theis_pumping_with_deepening_extended_extraction.Analytical(S, m, K, WL, R, IrrDepth,
        #                                                                                      years, option,
        #                                                                                      energy_cost)  # S, m, K, WL, R, IrrDepth, years, extended, cost):
        gw_cost_curve_internal = Superwell_for_ABM_on_the_fly.cost_curve(S, m, K, WL, R, IRR_DEPTH, NUM_YEARS, ELECTRICITY_RATE)

    # if the Superwell_for_ABM_on_the_fly fails to generate a cost curve (returns None), return all nan values
    if not gw_cost_curve_internal:
        return [np.nan] * 22

    # apply sensitivity multipliers
    gammas_total_subset = gammas_total_subset_og.copy()
    # for key in gammas_total_subset:
    #     gammas_total_subset[key] *= gammas_mult

    nirs_subset = nirs_subset_og.copy()
    for key in nirs_subset:
        nirs_subset[key] *= denom_hydro_factor

    sw_calib_constr_subset = sw_calib_constr_subset_og.copy()
    for key in sw_calib_constr_subset:
        sw_calib_constr_subset[key] *= numer_hydro_factor

    net_prices_land_subset = net_prices_land_subset_og.copy()
    for key in net_prices_land_subset:
        net_prices_land_subset[key] = ((yields_subset[key] * prices_subset[key] * numer_econ_factor) - land_costs_subset[key] -
                                       alphas_total_subset[key])

    #### Added by SF
    keys = net_prices_gw_subset.keys()
    net_prices_gw_sensitivity = net_prices_gw_subset.copy()
    for key in keys:
        net_prices_gw_sensitivity[key] = net_prices_gw_subset[key] * denom_econ_factor

    # for key in net_prices_land_subset:
    #     net_prices_land_subset[key] = ((yields_subset[key] * prices_subset[key]) - land_costs_subset[key] - (alphas_total_subset[key] * alphas_mult))

    # initialize counters/trackers
    first = True
    depletion_first = True
    cumul_gw_sum = 0
    gw_multiplier = 1
    gw_cost_added = 0
    gw_availability_multiplier = 1  # SF added for reduced GW availability cost curve option
    time_depletion = 0  # SF added to reset time to depletion

    # start simulation time loop
    for t in range(no_of_years):
        print('year' + str(t))
        print('!!JY!! ' + str(denom_hydro_factor))

        # initialize pyomo model
        fwm_s = ConcreteModel()
        fwm_s.ids = Set(initialize=ids_subset_sorted)
        if first:
            fwm_s.net_prices_gw = Param(fwm_s.ids, initialize=net_prices_gw_sensitivity, mutable=True)
        else:
            net_prices_gw_subset_update = net_prices_gw_sensitivity.copy()

            for key in net_prices_gw_subset_update:
                net_prices_gw_temp = net_prices_gw_subset_update[key]
                nir_temp = nirs_subset[key]
                gw_cost_temp = (-1.0 * net_prices_gw_temp) / (nir_temp * 1000.0)
                gw_cost_updated_temp = gw_cost_temp + (gw_cost_added * 1233.48) # 1233.48 converts $/m3 from cost curves to $/acre-ft used in PMP
                net_prices_gw_updated_temp = -1.0 * gw_cost_updated_temp * nir_temp * 1000.0
                net_prices_gw_subset_update[key] = net_prices_gw_updated_temp
            # net_prices_gw_subset_update.update((x, y * gw_multiplier) for x, y in net_prices_gw_subset_update.items())
            fwm_s.net_prices_gw = Param(fwm_s.ids, initialize=net_prices_gw_subset_update, mutable=True)
        fwm_s.farm_ids = Set(initialize=farm_id)
        fwm_s.crop_ids_by_farm = Set(fwm_s.farm_ids, initialize=crop_ids_by_farm_subset)
        fwm_s.crop_ids_by_farm_and_constraint = Set(fwm_s.farm_ids, initialize=crop_ids_by_farm_subset)
        fwm_s.net_prices_sw = Param(fwm_s.ids, initialize=net_prices_sw_subset, mutable=True)
        fwm_s.net_prices_total = Param(fwm_s.ids, initialize=net_prices_land_subset, mutable=True)
        fwm_s.gammas_total = Param(fwm_s.ids, initialize=gammas_total_subset, mutable=True)
        fwm_s.land_constraints = Param(fwm_s.farm_ids, initialize=max_land_constr_subset, mutable=True)
        fwm_s.sw_constraints = Param(fwm_s.farm_ids, initialize=sw_calib_constr_subset, mutable=True)

        ########## Added by SF to apply reduced GW availability multiplier to the gw_constraint
        gw_calib_constr_subset_copy = gw_calib_constr_subset.copy()
        gw_calib_constr_subset_copy[farm_id[0]] = gw_calib_constr_subset_copy[farm_id[0]] * gw_availability_multiplier
        fwm_s.gw_constraints = Param(fwm_s.farm_ids, initialize=gw_calib_constr_subset_copy,
                                     mutable=True)  # check that this is correct way to implement
        #########

        x_start_values = dict(enumerate([0.0] * 3))  # JY version from GitHub initialized with [0,0,0]
        fwm_s.xs_total = Var(fwm_s.ids, domain=NonNegativeReals, initialize=x_start_values[0])
        fwm_s.xs_sw = Var(fwm_s.ids, domain=NonNegativeReals, initialize=x_start_values[1])
        fwm_s.xs_gw = Var(fwm_s.ids, domain=NonNegativeReals, initialize=x_start_values[2])
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
        # try:
        results = opt.solve(fwm_s, keepfiles=False, tee=True)
        print(results.solver.termination_condition)

        # except ValueError:
        #   continue

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
        results_pd['year'] = t + 1

        ################# line for gw volume pumped in timestep ##################
        results_pd['gw_vol'] = results_pd['xs_gw'] * results_pd['nir_corrected'] * denom_hydro_factor  # area * nir (depth?) = volume
        ##########################################################################

        try:
            results_pd['gw_mult'] = gw_multiplier
            results_pd['gw_cost_added'] = gw_cost_added
        except NameError:
            # results_pd['gw_mult'] = 1
            results_pd['gw_cost_added'] = 0

        try:
            results_pd['gw_available'] = gw_availability_multiplier
        except NameError:
            results_pd['gw_available'] = 1

        try:
            results_pd['gw_cumul_vol'] = cumul_gw_sum
        except NameError:
            results_pd['gw_cumul_vol'] = 0
        # results_pd['gw_run'] = gw_cost_id

        # Creates and appends annual results
        if first:
            results_combined = results_pd
        else:
            results_combined = pd.concat([results_combined, results_pd])

        aggregation_functions = {'gw_vol': 'sum', 'xs_gw': 'sum'}
        gw_sum = results_pd[['nldas', 'gw_vol', 'xs_gw']].groupby(['nldas']).aggregate(aggregation_functions)
        gw_sum = gw_sum.reset_index()
        # gw_sum = gw_sum.merge(farm_gw_sum[['nldas', 'no_of_wells']], how='left', on=['nldas'])
        # gw_sum['gw_vol_well'] = gw_sum['gw_vol'] * farms_per_grid / gw_sum['no_of_wells']
        gw_sum['gw_vol_km3'] = gw_sum['gw_vol'] * 1.23348e-6
        if first:
            cumul_gw_sum = gw_sum['gw_vol_km3'].values[0]
        else:
            cumul_gw_sum += gw_sum['gw_vol_km3'].values[0]
        i = 0

        if gw_curve_option == 'internal':
            for index in range(gw_cost_curve_internal[0].size):
                if cumul_gw_sum > gw_cost_curve_internal[1][
                    index]:  # closest match for cum pumped vol cost curve bin
                    i = index
                else:
                    break
            if gw_cost_curve_internal[0][i] != 0:
                gw_cost_added = gw_cost_curve_internal[0][i] - gw_cost_curve_internal[0][
                0]  # unit cost divided by starting unit cost
                gw_multiplier = gw_cost_curve_internal[0][i] / gw_cost_curve_internal[0][
                    0]  # unit cost divided by starting unit cost
                gw_availability_multiplier = gw_cost_curve_internal[3][i]

            else:
                gw_cost_added = 9999999999999
                gw_multiplier = 9999999999999
                if depletion_first:
                    time_depletion = t
                    depletion_first = False


            # if option == 'reduced_capacity':  # SF added for reduced capacity option
            #     gw_availability_multiplier = gw_cost_curve_internal[3][i] / gw_cost_curve_internal[3][
            #         0]  # update with index for gw mult

        print(cumul_gw_sum)
        first = False

    # Calculate profit
    results_combined['net_prices_gw_updated'] = -1.0 * (results_combined['gw_cost'] + (gw_cost_added * 1233.48)) * results_combined['nir_corrected'] * 1000.0
    results_combined['profit_calc'] = ((numer_econ_factor * (results_combined['price'] * results_combined['yield']) -
                                        results_combined['land_cost'] - results_combined['alphas_land']) *
                                       results_combined['xs_total']) \
                                      - (0.5 * results_combined['gammas_total'] * results_combined['xs_total'] *
                                         results_combined['xs_total']) \
                                      - (results_combined['net_prices_gw_updated'] * results_combined['xs_gw']) \
                                      - (results_combined['net_prices_sw'] * results_combined['xs_sw'])

    #  Calculate summary results (percent change in total irrigated area as an initial result)
    aggregation_functions = {'xs_total': 'sum', 'profit_calc': 'sum', 'gw_vol':'sum','gw_cost_added':'mean', 'gw_available':'mean'}
    summary_pd = results_combined[['nldas', 'year', 'xs_total', 'profit_calc','gw_vol','gw_cost_added','gw_available']].groupby(['nldas', 'year']).aggregate(
        aggregation_functions)
    summary_pd = summary_pd.reset_index()
    year_start = summary_pd.year.min()
    year_end = summary_pd.year.max()
    acres_grown_total = summary_pd.xs_total[:].sum()
    acres_grown_mean = summary_pd.xs_total[:].mean()
    acres_grown_mean_vs_start = summary_pd.xs_total[:].mean() / \
                                summary_pd[(summary_pd.year == year_start)].xs_total.values[0]
    end_div_start_area = summary_pd[(summary_pd.year == year_end)].xs_total.values[0] / \
                         summary_pd[(summary_pd.year == year_start)].xs_total.values[0]
    min_acres_grown = summary_pd.xs_total.min()
    max_acres_grown = summary_pd.xs_total.max()
    max_minus_min_acres_grown = summary_pd.xs_total.max() - summary_pd.xs_total.min()
    med_acres_grown = summary_pd.xs_total.median()
    total_profit = summary_pd.profit_calc.sum()
    min_annual_profit = summary_pd.profit_calc.min()
    max_annual_profit = summary_pd.profit_calc.max()
    max_minus_min_annual_profit = summary_pd.profit_calc.max() - summary_pd.profit_calc.min()
    med_annual_profit = summary_pd.profit_calc.median()
    end_div_start_profit = summary_pd[(summary_pd.year == year_end)].profit_calc.values[0] / \
                           summary_pd[(summary_pd.year == year_start)].profit_calc.values[0]
    max_gw_vol_yr = summary_pd.gw_vol.max()
    min_gw_vol_yr = summary_pd.gw_vol.min()
    max_minus_min_gw_vol_yr = summary_pd.gw_vol.max() - summary_pd.gw_vol.min()
    med_gw_vol_yr = summary_pd.gw_vol.median()
    mean_excl0_gw_vol_yr = summary_pd[(summary_pd.gw_vol != 0)].gw_vol.mean()

    # Volume depleted (calculated from cost curves)
    grid_cell_area = 13.875 * 13.875 * 10 ** 6  # from Stephen's script -- eventually replace with actual cell area
    available_volume = (m - WL) *  grid_cell_area * S
    perc_vol_depleted = cumul_gw_sum * (10 ** 9) / available_volume
    # perc_vol_depleted = cumul_gw_sum * (10**9) / max(gw_cost_curve_internal[1].tolist())

    if perc_vol_depleted > 1:
        perc_vol_depleted = 1

    # return summary_pd

    return [end_div_start_area, total_profit, perc_vol_depleted, time_depletion, cumul_gw_sum, acres_grown_total,
            acres_grown_mean, acres_grown_mean_vs_start, end_div_start_profit, min_acres_grown, max_acres_grown,
            max_minus_min_acres_grown, med_acres_grown, min_annual_profit, max_annual_profit,
            max_minus_min_annual_profit, med_annual_profit, max_gw_vol_yr, min_gw_vol_yr, max_minus_min_gw_vol_yr,
            med_gw_vol_yr, mean_excl0_gw_vol_yr]


##### Travis: Lines 468-523 are used to define the ensemble information, which is consolidated in the "cases_df" pandas dataframe
# Varied parameters for ABM sensitivity
hydro_ratio = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5]
econ_ratio = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5]
# K_val = [1, 100]
K_val = [1]
m_val = [m] # default in one m from single cost curve as defined in 'Cost curve inputs'
single_cost_curve = 'false'

##### Travis - Here is where I provide a list of farms by ID that we subsequently loop through. For the comprehensive farm sweep on HPC, the idea would be to run for every farm in our database (0-53834) and strategize a way to parallelize this on HPC (each call of farm_gw_model is independent)
num_farms = 2
farms = [36335, 15557]
# farms = [15557, 36335]  #15557 OG (x200y128 - nearly all corn and fodder grass / gw cost 10% net price), 36335 CA (x32y103 - majority misc crop / gw cost 30% net price), 28380  (x273y77 - major oil crop, secondary fiber, corn / gw cost 20% net price)
# farms = [36335, 15557, 28380]

# Create array of ABM sensitivity parameter combinations using parameter lists defined in 'Varied parameters for ABM sensitivity'
combinations = len(hydro_ratio) * len(econ_ratio) * len(K_val) * len(m_val) * len(farms)
cases = list(product(hydro_ratio, econ_ratio, K_val, m_val, farms))
cases_df = pd.DataFrame(data = cases.copy(), columns = ['hydro_ratio', 'econ_ratio', 'K_val', 'm_val', 'farm'])

# Replace arbitrary K values with those from NLDAS data
cases_df['K_val_data'] = 0
for f in farms:
    nldas_id = farms_master.iloc[f]['nldas']
    K_low = nldas_gw_attributes[(nldas_gw_attributes['NLDAS_ID'] == nldas_id)]['K'].values[0]
    K_high = nldas_gw_attributes[(nldas_gw_attributes['NLDAS_ID'] == nldas_id)]['K_hi'].values[0]
    cases_df.loc[(cases_df.farm == f) & (cases_df.K_val == 1), 'K_val_data'] = K_high
    # cases_df.loc[(cases_df.farm == f) & (cases_df.K_val == 100), 'K_val_data'] = K_high

# Enforce minimum K for viable pumping
cases_df.loc[(cases_df.K_val_data < 0.10), 'K_val_data'] = 0.10

### Create Output arrays to store resulting data
results_end_div_start_area = []
results_total_profit = []
results_perc_vol_depleted = []
results_time_depletion = []
results_cum_gw = []
# New metrics added by SF
acres_grown_total = []
acres_grown_mean = []
acres_grown_mean_vs_start =  []
end_div_start_profit =  []
# New metrics added by JY
min_acres_grown = []
max_acres_grown = []
max_minus_min_acres_grown = []
med_acres_grown = []
min_annual_profit = []
max_annual_profit = []
max_minus_min_annual_profit = []
med_annual_profit = []
max_gw_vol_yr = []
min_gw_vol_yr = []
max_minus_min_gw_vol_yr = []
med_gw_vol_yr = []
mean_excl0_gw_vol_yr = []

### Travis - this is where we start the simulations, looping through every case in the cases_df dataframe
### Simulate parameter combinations over all selected farms
for case in range(combinations):
    print('!!JY!! run #: ' + str(case))
    results = farm_gw_model(cases_df.hydro_ratio[case], cases_df.econ_ratio[case],
                            cases_df.K_val_data[case], cases_df.farm[case], cases_df.m_val[case], R, gw_cost_curve_internal)

    results_end_div_start_area.append(results[0])
    results_total_profit.append(results[1])
    results_perc_vol_depleted.append(results[2])
    results_time_depletion.append(results[3])
    results_cum_gw.append(results[4])
    acres_grown_total.append(results[5])
    acres_grown_mean.append(results[6])
    acres_grown_mean_vs_start.append(results[7])
    end_div_start_profit.append(results[8])
    min_acres_grown.append(results[9])
    max_acres_grown.append(results[10])
    max_minus_min_acres_grown.append(results[11])
    med_acres_grown.append(results[12])
    min_annual_profit.append(results[13])
    max_annual_profit.append(results[14])
    max_minus_min_annual_profit.append(results[15])
    med_annual_profit.append(results[16])
    max_gw_vol_yr.append(results[17])
    min_gw_vol_yr.append(results[18])
    max_minus_min_gw_vol_yr.append(results[19])
    med_gw_vol_yr.append(results[20])
    mean_excl0_gw_vol_yr.append(results[21])

##### Travis: Finally, we save the results as part of the cases df dataframe
cases_df['end_div_start_area'] = results_end_div_start_area
cases_df['total_profit'] = results_total_profit
cases_df['perc_vol_depleted'] = results_perc_vol_depleted
cases_df['time_depletion'] = results_time_depletion
cases_df['cumul_gw'] = results_cum_gw
cases_df['acres_grown_total'] = acres_grown_total
cases_df['acres_grown_mean'] = acres_grown_mean
cases_df['acres_grown_mean_vs_start'] = acres_grown_mean_vs_start
cases_df['end_div_start_profit'] = end_div_start_profit
cases_df['min_acres_grown'] = min_acres_grown
cases_df['max_acres_grown'] = max_acres_grown
cases_df['max_minus_min_acres_grown'] = max_minus_min_acres_grown
cases_df['med_acres_grown'] = med_acres_grown
cases_df['min_annual_profit'] = min_annual_profit
cases_df['max_annual_profit'] = max_annual_profit
cases_df['min_div_max_annual_profit'] = cases_df['min_annual_profit'] / cases_df['max_annual_profit']
cases_df['max_minus_min_annual_profit'] = max_minus_min_annual_profit
cases_df['med_annual_profit'] = med_annual_profit
cases_df['max_gw_vol_yr'] = max_gw_vol_yr
cases_df['min_gw_vol_yr'] = min_gw_vol_yr
cases_df['max_minus_min_gw_vol_yr'] = max_minus_min_gw_vol_yr
cases_df['med_gw_vol_yr'] = med_gw_vol_yr
cases_df['mean_excl0_gw_vol_yr'] = mean_excl0_gw_vol_yr

##### Travis: And export to csv. For each farm (~50k), we will have 121 scenarios, and for each farm/scenario ~30 results that are stored (~180M values)
cases_df.to_csv('20230719_temp.csv')


