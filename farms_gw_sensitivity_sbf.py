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
import os
import matplotlib.pyplot as plt

os.chdir('C:/Users/fere556/Documents/Farm_GW_Sensitivity/farm_gw_sensitivity-main')

##### Load Theis drawdown response module
import Theis_pumping_with_deepening_extended_extraction

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


##### Load external files for GW cost curves
gw_cost_curves = pd.read_csv('./data_inputs/20221005_cost_curves/20221005_gw_cost_curves.csv')
gw_cost_lookup = pd.read_csv('./data_inputs/20221005_cost_curves/20221005_gw_cost_lookup.csv')

#%%
##### Select farm/gw run options (these are currently implemented as higher-level options than the sensitivity inputs)
farm_id = [15557]  # Pick the farm agent [0-53834] to use as a basis for the sensitivity analysis
gw_cost_id = 'gw4'  # Pick the gw cost curve to use as a basis for the sensitivity analysis
# gw_area_well_sqm = 750 * 750  # Assumed groundwater area irrigated with a well in square meters
# gw_area_well_acres = gw_area_well_sqm * 0.000247105  # Assumed groundwater area irrigated with a well in acres
no_of_years = 50  # The number of years (annual timesteps) to run simulation
# farms_per_grid = 4412  # Assumed number of farms per NLDAS grid cell (50 km x 50 km) / 140 acres = 4412

##### Select groundwater cost curve generation options
# gw_curve_option = 'internal'  # GW cost curve option ("external" - pre-processed and pulled from external file, "internal" - generate on the fly)
# S = 0.25
# m = 130
# K = 1 #define in sensitivity analysis
# WL = 10
# R = 0
# IrrDepth = 12
# years = 80
# option = 'reduced_capacity'
# energy_cost = 0.2

##### For each representative farm, estimate # of groundwater wells (taking initial total groundwater area divided by 750m x 750m)
aggregation_functions = {'area_irrigated_gw': 'sum'}
farm_gw_sum = farms_master[['nldas','area_irrigated_gw']].groupby(['nldas']).aggregate(aggregation_functions)
farm_gw_sum = farm_gw_sum.reset_index()
# farm_gw_sum['no_of_wells'] = farm_gw_sum['area_irrigated_gw'] / (gw_area_well_acres)  # .000247105 conversion from square meters to acres

##### Subset relevant data inputs for the specific farm/gw run option
# subset crop ids
crop_ids_by_farm_subset = {key: crop_ids_by_farm[key] for key in farm_id} # farm_id = row in id_nldas, dictionary with farm_id as key and list of crops as value
ids_subset = []
for key, list_value in crop_ids_by_farm_subset.items(): # creates list where first entry is farm_id and each follow entry is a crop index which is 
    for value in list_value:                            # row in farms_master DataFrame 
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

#%% Define farm_gw_model function that enables varied ABM parameters and 
# aquifer properties (for now K and m (thickness)). This function can use a 
# internally generated cost curve defined by a single set of aquifer properties (block 
# of code following function definition) OR can generate a new cost curve within 
# the function if both the gw_curve_option = 'internal" AND single_cost_curve = 'false' flags
# are defined. The Sobol analysis requires generating new cost curves for each function call
# and therefore gw_curve_option = 'internal" AND single_cost_curve = 'false' should be set 
# for SALib analyses. 

#del farm_id  

def farm_gw_model(sw_mult, gw_mult, price_mult, gammas_mult, K, farm_id, m, gw_cost_curve_internal):
    farm_id = [farm_id]  # Pick the farm agent [0-53834] to use as a basis for the sensitivity analysis
    #farm_id = [12641] # use if running a single farm, enter farm_id 
    gw_cost_id = 'gw4'  # Pick the gw cost curve to use as a basis for the sensitivity analysis
    # gw_area_well_sqm = 750 * 750  # Assumed groundwater area irrigated with a well in square meters
    # gw_area_well_acres = gw_area_well_sqm * 0.000247105  # Assumed groundwater area irrigated with a well in acres
    no_of_years = 100  # The number of years (annual timesteps) to run simulation, can be toggled 
    # farms_per_grid = 4412  # Assumed number of farms per NLDAS grid cell (50 km x 50 km) / 140 acres = 4412
      
    ##### For each representative farm, estimate # of groundwater wells (taking initial total groundwater area divided by 750m x 750m)
    aggregation_functions = {'area_irrigated_gw': 'sum'}
    farm_gw_sum = farms_master[['nldas','area_irrigated_gw']].groupby(['nldas']).aggregate(aggregation_functions)
    farm_gw_sum = farm_gw_sum.reset_index()
    # farm_gw_sum['no_of_wells'] = farm_gw_sum['area_irrigated_gw'] / (gw_area_well_acres)  # .000247105 conversion from square meters to acres
    
    ##### Subset relevant data inputs for the specific farm/gw run option
    # subset crop ids
    crop_ids_by_farm_subset = {key: crop_ids_by_farm[key] for key in farm_id} # farm_id = row in id_nldas, dictionary with farm_id as key and list of crops as value
    ids_subset = []
    for key, list_value in crop_ids_by_farm_subset.items(): # creates list where first entry is farm_id and each follow entry is a crop index which is 
        for value in list_value:                            # row in farms_master DataFrame 
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
    
    #time_load = time.time()
    
    # alphas_mult = 1  # JY temp to disable alphas in sensitivity analysis
    nir_mult = 1 # JY temp to disable alphas in sensitivity analysis
    
    # generate gw cost curves if more than one cost curve selected in sensitivity sweep 
    if gw_curve_option == 'internal' and single_cost_curve != 'true':
        gw_cost_curve_internal = Theis_pumping_with_deepening_extended_extraction.Analytical(S, m, K, WL, R, IrrDepth, years, option, energy_cost) #S, m, K, WL, R, IrrDepth, years, extended, cost):
        
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

    #### Added by SF 
    keys = net_prices_gw_subset.keys()
    net_prices_gw_sensitivity = net_prices_gw_subset.copy()
    for key in keys:
        net_prices_gw_sensitivity[key] =  net_prices_gw_subset[key]*gw_mult
    
    # for key in net_prices_land_subset:
    #     net_prices_land_subset[key] = ((yields_subset[key] * prices_subset[key]) - land_costs_subset[key] - (alphas_total_subset[key] * alphas_mult))

    # initialize counters/trackers
    first = True
    depletion_first = True
    cumul_gw_sum = 0
    gw_multiplier = 1
    gw_availability_multiplier = 1 # SF added for reduced GW availability cost curve option 
    time_depletion = 0 # SF added to reset time to depletion 

    # start simulation time loop
    for t in range(no_of_years):
        print('year' + str(t))
        print('!!JY!! ' + str(nir_mult))

        # initialize pyomo model
        fwm_s = ConcreteModel()
        fwm_s.ids = Set(initialize=ids_subset_sorted)
        if first:
            fwm_s.net_prices_gw = Param(fwm_s.ids, initialize = net_prices_gw_sensitivity, mutable=True)
        else:
            net_prices_gw_subset_update = net_prices_gw_sensitivity.copy()
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
        
        ########## Added by SF to apply reduced GW availability multiplier to the gw_constraint
        gw_calib_constr_subset_copy = gw_calib_constr_subset.copy()
        gw_calib_constr_subset_copy[farm_id[0]] = gw_calib_constr_subset_copy[farm_id[0]] * gw_availability_multiplier
        fwm_s.gw_constraints = Param(fwm_s.farm_ids, initialize=gw_calib_constr_subset_copy, mutable=True) # check that this is correct way to implement 
        #########
        
        x_start_values = dict(enumerate([0.0] * 3)) # JY version from GitHub initialized with [0,0,0]
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
        #try:
        results = opt.solve(fwm_s, keepfiles=False, tee=True)
        print(results.solver.termination_condition)
            
        #except ValueError:
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
        results_pd['year'] = t+1
        
        ################# line for gw volume pumped in timestep ##################
        results_pd['gw_vol'] = results_pd['xs_gw'] * results_pd['nir_corrected'] # area * nir (depth?) = volume
        ##########################################################################
        
        try:
            results_pd['gw_mult'] = gw_multiplier
        except NameError:
            results_pd['gw_mult'] = 1
            
        try:
            results_pd['gw_available'] = gw_availability_multiplier
        except NameError:
            results_pd['gw_available'] = 1
            
        try:
            results_pd['gw_cumul_vol'] = cumul_gw_sum
        except NameError:
            results_pd['gw_cumul_vol'] = 0
        results_pd['gw_run'] = gw_cost_id
        results_pd['trans'] = gw_cost_lookup[(gw_cost_lookup.row_indx == gw_cost_id)]['Trans'].values[0]
        
        # Creates and appends annual results
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
                if cumul_gw_sum > gw_cost_curve_internal[1][0][index]: # closest match for cum pumped vol cost curve bin
                    i = index
                else:
                    break
            if gw_cost_curve_internal[0][0][i] != 0:
                gw_multiplier = gw_cost_curve_internal[0][0][i] / gw_cost_curve_internal[0][0][0] # unit cost divided by starting unit cost                 
            
            else:
                gw_multiplier = 9999999999999
                if depletion_first:
                    time_depletion = t
                    depletion_first = False
                    
            if option == 'reduced_capacity': # SF added for reduced capacity option 
                gw_availability_multiplier = gw_cost_curve_internal[3][i] / gw_cost_curve_internal[3][0] # update with index for gw mult 

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
    acres_grown_total = summary_pd.xs_total[:].sum()
    acres_grown_mean = summary_pd.xs_total[:].mean()
    acres_grown_mean_vs_start = summary_pd.xs_total[:].mean()/summary_pd[(summary_pd.year == year_start)].xs_total.values[0]
    end_div_start_area = summary_pd[(summary_pd.year==year_end)].xs_total.values[0] / summary_pd[(summary_pd.year == year_start)].xs_total.values[0]
    total_profit = summary_pd.profit_calc.sum()
    end_div_start_profit = summary_pd[(summary_pd.year==year_end)].profit_calc.values[0] / summary_pd[(summary_pd.year == year_start)].profit_calc.values[0]

    # Volume depleted (calculated from cost curves)
    perc_vol_depleted = cumul_gw_sum / max(gw_cost_curve_internal[1][0].tolist())
    if perc_vol_depleted > 1:
        perc_vol_depleted = 1
    
    
    return [end_div_start_area, total_profit, perc_vol_depleted, time_depletion, cumul_gw_sum, acres_grown_total, 
            acres_grown_mean, acres_grown_mean_vs_start, end_div_start_profit] 

#results_combined.to_csv('Farm_12641_restricted.csv')

#%% Exploratory simulations for any number of farm cellst, SF mostly used this
# block to evaluate farm sensitivity to GW cost by using the output from a 
# single cost curve generated from parameters specified in Lines 469-478
# and using the cost curve to explore GW responsiveness accross ranges of 
# sw_mult, gw_mult, price_mult, and gammas_mult defined in Lines 486-490 

# Cost curve inputs 
gw_curve_option = 'internal'  # GW cost curve option ("external" - pre-processed and pulled from external file, "internal" - generate on the fly)
S = 0.25
m = 110 
K = 10 # define in sensitivity analysis
WL = 10
R = 0
IrrDepth = 12
years = 80
option = 'reduced_capacity'
energy_cost = 0.15

single_cost_curve = 'true' # flag to skip internal generation of cost curve in farm_gw_model function 

# Generate cost curve 
gw_cost_curve_internal = Theis_pumping_with_deepening_extended_extraction.Analytical(S, m, K, WL, R, IrrDepth, years, option, energy_cost)

# Vaired parameters for ABM sensitivity 
sw_mult = [1]
gw_mult = [1,2,3]
price_mult = [1]
gammas_mult = [1]
K_val = [K] # default is one K from single cost curve as defined in 'Cost curve inputs'
m_val = [m] # default in one m from single cost curve as defined in 'Cost curve inputs'

if len(K_val)*len(m_val > 1): # if more than one set of aquifer K and m, flag set to 'false' and cost curves generated inside the farm_gw_model function 
    single_cost_curve = 'false'
    
# Create array of ABM sensitivity parameter combinations using parameter lists defined in 'Varied parameters for ABM sensitivity'
combinations = len(sw_mult) * len(gw_mult) * len(price_mult) * len(gammas_mult) * len(K_val) * len(m_val)
cases = np.zeros((combinations,6))
sw_arr = np.sort(np.tile(sw_mult, int(combinations/len(sw_mult))))
gw_arr = np.tile(gw_mult, int(combinations/len(gw_mult)))
cases[:,0] = sw_arr[:]
cases[:,1] = gw_arr[:]
price_arr = np.sort(np.tile(price_mult, len(gw_mult)))
price_arr = np.tile(price_arr, int(combinations/len(price_arr)))
cases[:,2] = price_arr[:]
gammas_arr = np.sort(np.tile(gammas_mult, int(int(combinations/len(gammas_mult))/len(gammas_mult))))
gammas_arr = np.tile(gammas_arr,int(combinations/len(gammas_arr)))
cases[:,3] = gammas_arr[:]
K_arr = np.sort(np.tile(K_val, int(combinations/len(K_val))))
cases[:,4] = K_arr[:]
m_arr = np.sort(np.tile(m_val, int(int(combinations/len(m_val))/len(m_val))))
m_arr = np.tile(m_arr, int(combinations/len(m_arr)))
cases[:,5] = m_arr[:]

cases_df = pd.DataFrame(data = cases.copy(), columns = ['sw_mult', 'gw_mult', 'price_mult', 'gammas_mult', 'K_val', 'm_val'])

### Generate list of farm IDs to evaluate 
farm_gw_sum_sorted = farm_gw_sum.copy() # create copy of farm irrigated GW area

# Sort by GW area in descending order 
farm_gw_sum_sorted = farm_gw_sum_sorted.sort_values(by = ['area_irrigated_gw'],
                                                    ascending = False)
# Copy index of sorted GW area dataframe (farm_id)
farm_indexes = farm_gw_sum_sorted.index.values

# Sample farm indexes 
sample = np.linspace(0,999,2).astype('int') # in this example first 1000 largest GW areas 
farms = farm_indexes[sample]
num_farms = len(farms) # number of farms sampled
gw_areas = farm_gw_sum_sorted.area_irrigated_gw[farms] # gw area for each selected farm 
plt.hist(gw_areas, bins = np.linspace(0,50000, 20)) # plot histogram of sample GW areas 
plt.xlabel('GW_Area (acres)')
plt.ylabel('Count')

### Create Output arrays to store data
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

count = 0 

### Simulate parameter combinations over all selected farms 
for j in range(num_farms):
    for case in range(combinations):

        results = farm_gw_model(cases_df.sw_mult[case], cases_df.gw_mult[case], 
                    cases_df.price_mult[case], cases_df.gammas_mult[case], 
                    cases_df.K_val[case], farms[j], cases_df.m_val[case], gw_cost_curve_internal)
        
        results_end_div_start_area.append(results[0])
        results_total_profit.append(results[1])
        results_perc_vol_depleted.append(results[2])
        results_time_depletion.append(results[3])
        results_cum_gw.append(results[4])
        acres_grown_total.append(results[5])
        acres_grown_mean.append(results[6])
        acres_grown_mean_vs_start.append(results[7])
        end_div_start_profit.append(results[8])
        
        count +=1
            
# Plot results - only designed to plot when GW cost is the only parameter varied
# I used this to plot the GW cost responsiveness results from my initial analysis 
fig, axs = plt.subplots(nrows = 2, ncols = 2, figsize = (7,6))
colors = plt.cm.rainbow(np.linspace(0, 1, num_farms))
for i in range(num_farms):
    
    if i == 0:
        lower = 0
    else:
        lower = i*combinations
        
    upper = (i+1)*combinations
    
    x_val = cases[:,1] 
    axs[0,0].plot(x_val, np.array(results_total_profit[lower:upper]/max(results_total_profit[lower:upper])), color = colors[i], alpha = 1)
    axs[0,0].set_ylabel('Profit/Max_Profit')
    axs[0,1].plot(x_val, results_end_div_start_area[lower:upper], color = colors[i])
    axs[0,1].set_ylabel('End/Start Area')
    axs[1,0].plot(x_val, results_perc_vol_depleted[lower:upper], color = colors[i])
    #axs[1,0].plot(x_val, np.array(results_perc_vol_depleted[lower:upper])/max(results_perc_vol_depleted[lower:upper]), color = colors[i])
    axs[1,0].set_ylabel('% Vol Depleted')
    axs[1,1].plot(x_val, results_time_depletion[lower:upper], color = colors[i])
    axs[1,1].set_ylabel('Years to Depletion')

axs[1,0].set_xlabel('Initial GW Mult')
axs[1,1].set_xlabel('Initial GW Mult')
plt.tight_layout()


### Save farm sensitivity results 

# Save crop data for each farm 
farm_props = farms_master.iloc[farms,:]
for i in range(9):
    farm_props = pd.concat([farm_props,farms_master.iloc[farms+53835*(i+1),:]])

sensitivity_params = np.tile(cases, (num_farms,1))

# make farm id labels for output table
farm_id_labels = np.zeros(combinations*num_farms)
for i in range(num_farms):
    farm_id_labels[i*combinations:(i+1)*combinations] = farms[i]
    
results = np.array([farm_id_labels[:], results_end_div_start_area, results_total_profit,
                    results_perc_vol_depleted, results_time_depletion, acres_grown_total, 
                            acres_grown_mean, acres_grown_mean_vs_start, end_div_start_profit])

results = results.T

name = 'example' # change name for saved output prefix

# Uncomment to save outputs
# np.savetxt(name + '_results.csv', results, fmt = '%10.5f', delimiter = ',') # results 
# farm_props.to_csv(name +'_farm_props.csv') # properties from sampled farms
# np.savetxt(name +'_sensitivity_params.csv', sensitivity_params, fmt = '%10.5f', delimiter = ',') # parameter values for each case

#%% Analysis of farm sensitivity results 

# Load farm sensitivity results 
results = np.loadtxt('all_top_1000_results.csv', delimiter = ',', skiprows = 1) # load results array 
farm_props = pd.read_csv('all_top_1000_farm_props.csv', delimiter = ',', index_col = 0) # load farm properties of farms sampled 


farm_num = int(farm_props.iloc[:,0].size/10) # number of farms assessed 
combinations = int(len(results[:,0])/farm_num) # of sensitivity param combinations for each farm 

gw_depletion_variation = np.zeros(farm_num) # difference in gw depletion for each farm accorss the various case combinations 
farms = np.asarray(farm_props.index[0:farm_num], 'int') # list of farm_ids that were evaluated 

# Prescribe list of farm IDs
farm_gw_sum_sorted = farm_gw_sum.copy()
farm_gw_sum_sorted = farm_gw_sum_sorted.sort_values(by = ['area_irrigated_gw'],
                                                    ascending = False)

gw_areas = farm_gw_sum_sorted.area_irrigated_gw[farms] # gw area for each selected farm 


for i in range(farm_num):
    gw_depletion_variation[i] = max(results[i*combinations:(i+1)*combinations,3]) \
                               -min(results[i*combinations:(i+1)*combinations,3])


# Dataframe with farm_id, NLDAS grid name, and GW sensitivity 
NLDAS_ids = farm_props.nldas[farms]
gw_response_NLDAS = pd.DataFrame(data = NLDAS_ids, index = farms)
gw_response_NLDAS['gw_response'] = gw_depletion_variation
gw_response_NLDAS.to_csv('top_1000_all_gw_response_GIS.csv')

avg_gw_price = np.zeros(farm_num)
avg_gw_price_weighted =  np.zeros(farm_num)
avg_sw_price = np.zeros(farm_num)
avg_sw_price_weighted =  np.zeros(farm_num)
area_ratio_gw_sw = np.zeros(farm_num)
avg_land_price_weighted = np.zeros(farm_num)
gammas_weighted = np.zeros(farm_num)
alphas_weighted =  np.zeros(farm_num)

crop_total = np.zeros((farm_num, 10))
crop_gw = np.zeros((farm_num, 10))
crop_sw = np.zeros((farm_num, 10))

indexes = 53835*np.ones(10)
for i in range(10):
    indexes[i] = 53835*i
    
for i in range(farm_num):
    locs = (farm_props.index[i]+indexes).astype(int)
    
    avg_gw_price[i] = farm_props.net_prices_gw[locs].mean()
    avg_sw_price[i] = farm_props.net_prices_sw[locs].mean()

    avg_gw_price_weighted[i] = (farm_props.net_prices_gw[locs]*farm_props.area_irrigated_gw[locs]).sum()/gw_areas[locs[0]]
    avg_sw_price_weighted[i] = (farm_props.net_prices_sw[locs]*farm_props.area_irrigated_sw[locs]).sum()/(farm_props.area_irrigated_sw[locs]).sum()

    area_ratio_gw_sw[i] = (farm_props.area_irrigated_gw[locs]).sum()/(farm_props.area_irrigated_sw[locs]).sum()

    avg_land_price_weighted[i] = (farm_props.net_prices_land[locs]*farm_props.area_irrigated[locs]).sum()/(farm_props.area_irrigated[locs]).sum()
    
    gammas_weighted[i] = (farm_props.gammas_total[locs]*farm_props.area_irrigated[locs]).sum()/(farm_props.area_irrigated[locs]).sum()
    
    alphas_weighted[i] = (farm_props.alphas_land[locs]*farm_props.area_irrigated[locs]).sum()/(farm_props.area_irrigated[locs]).sum()
    
    crop_sw[i,:] = farm_props.area_irrigated_sw[locs]
    crop_gw[i,:] = farm_props.area_irrigated_gw[locs]
    crop_total[i,:] = farm_props.area_irrigated_sw[locs] + farm_props.area_irrigated_gw[locs]


# determine mostly mono crop farms (largest crop fractional area) 
largest_crop_fraction = np.zeros(farm_num)
largest_crop_fraction_gw = np.zeros(farm_num)
largest_crop_fraction_sw = np.zeros(farm_num)
for i in range(farm_num):
    largest_crop_fraction[i] = max(crop_total[i,:])/np.sum((crop_total[i,:]))
    largest_crop_fraction_gw[i] = max(crop_gw[i,:])/np.sum((crop_gw[i,:]))
    largest_crop_fraction_sw[i] = max(crop_sw[i,:])/np.sum((crop_sw[i,:]))

# Simple x-y plots of model parameter on the x-axis and gw depletion response on the y  
plt.scatter(-1*avg_gw_price, gw_depletion_variation, alpha = 0.5)
plt.scatter(-1*avg_sw_price, gw_depletion_variation, alpha = 0.5)
plt.scatter(-1*avg_gw_price + avg_sw_price, gw_depletion_variation, alpha = 0.5)
plt.scatter(-1*avg_gw_price_weighted, gw_depletion_variation, alpha = 0.5)
plt.scatter(-1*avg_sw_price_weighted, gw_depletion_variation, alpha = 0.5)
plt.scatter(area_ratio_gw_sw, gw_depletion_variation, alpha = 0.5)
plt.xlim([0,50])
plt.scatter(avg_land_price_weighted, gw_depletion_variation, alpha = 0.5)
#plt.xlim([0,50])
plt.scatter(gammas_weighted, gw_depletion_variation, alpha = 0.1)
plt.scatter(-1*alphas_weighted, gw_depletion_variation, alpha = 0.25)
plt.scatter(crop_gw[:,3], gw_depletion_variation, alpha = 0.25)
plt.scatter(largest_crop_fraction_gw, gw_depletion_variation, alpha = 0.25)
plt.scatter(largest_crop_fraction, gw_depletion_variation, alpha = 0.25)
plt.scatter(largest_crop_fraction_sw, largest_crop_fraction_gw, alpha = 0.25)
plt.scatter(largest_crop_fraction, largest_crop_fraction_gw, alpha = 0.25)

#%%
##### Run sensitivity analysis using SALib
# problem = {
#     'num_vars': 5,
#     'names': ['nir_mult','sw_mult','price_mult','alphas_mult','gammas_mult'],
#     'bounds': [[0.8, 1.2], [0.8, 1.2], [0.8, 1.2], [0.8, 1.2], [0.8, 1.2]]
# }

# problem = {
#     'num_vars': 5,
#     'names': ['nir_mult','sw_mult', 'price_mult','K','gammas_mult'],
#     'bounds': [[0.8, 1.2], [0.8, 1.2], [0.8, 1.2], [1, 100], [0.8, 1.2]]
# }

# problem = {
#     'num_vars': 4,
#     'names': ['sw_mult','price_mult','gammas_mult','K'],
#     'bounds': [[0.8, 1.2], [0.8, 1.2], [0.8, 1.2], [1, 100]]
# }

from SALib import ProblemSpec

# If running SALib, single cost curve is false and cost curves will be generated
# using parameter values defined by the 'param_values' array
single_cost_curve = 'false'

#def crop_price_space(price_mult, gamma_mult, gw_mult, sw_mult):
    # Run crop_price_space vaiant of the farm-gw sensitivity 
    
    #return ....params
#farm_gw_model(sw_mult, gw_mult, price_mult, gammas_mult, K, farm_id):
    
problem = {'num_vars': 3, 'names':['gw_mult', 'depth', 'K'], 
                  'bounds': [[0.5, 1.5],[40, 110],[0.5, 10]]}

price_mult = np.linspace(0.6, 1.4, 5)

# Wrapper to run each price_mult scenario
param_values = saltelli.sample(problem, 14)

results_end_div_start_area = np.zeros((len(param_values[:,0]),len(price_mult)))
results_total_profit = np.zeros((len(param_values[:,0]),len(price_mult)))
results_perc_vol_depleted = np.zeros((len(param_values[:,0]),len(price_mult)))
results_time_depletion = np.zeros((len(param_values[:,0]),len(price_mult)))
results_cum_gw = np.zeros((len(param_values[:,0]),len(price_mult)))
    
for j in range(len(price_mult)):
    #K = 1
    farm_id = 44833 # farm to run Sobol
    count = 0
    sw_mult = 1 # not varied in this example of Sobol
    gammas_mult = 1 # not varied in this example of sobol
    results_end_div_start_area_list = []
    results_total_profit_list = []
    results_perc_vol_depleted_list = []
    results_time_depletion_list = []
    results_cum_gw_list = []
    
    for i in range(len(param_values[:,0])):
        results = farm_gw_model(sw_mult, param_values[i,0], price_mult[j], 
                                gammas_mult, param_values[i,2], farm_id, param_values[i,1], )
        results_end_div_start_area_list.append(results[0])
        results_total_profit_list.append(results[1])
        results_perc_vol_depleted_list.append(results[2])
        results_time_depletion_list.append(results[3])
        results_cum_gw_list.append(results[4])
        count += 1
    
    results_end_div_start_area[:,j] = np.array(results_end_div_start_area_list)
    results_total_profit[:,j] = np.array(results_total_profit_list)
    results_perc_vol_depleted[:,j] = np.array(results_perc_vol_depleted_list)
    results_time_depletion[:,j] = np.array(results_time_depletion_list)
    results_cum_gw[:,j] = np.array(results_cum_gw_list)
    
    #time_ensemble = time.time()

# Save results

# Analyze results
results_array = np.array(results_cum_gw)
S1_vol_dep = np.zeros((5,3))
ST_vol_dep = np.zeros((5,3))
for i in range(5):
    sobol_indices = sobol.analyze(problem, results_array[:,i], calc_second_order = True)
    S1_vol_dep[i,:] = sobol_indices['S1']
    ST_vol_dep[i,:] = sobol_indices['ST']
    
# Plots of S1s
plt.figure()
labels = ['GW mult','depth','K']
for i in range(3):
    plt.plot(price_mult, S1_vol_dep[:,i], label = labels[i])
plt.legend()

# Plots of STs
plt.figure()
labels = ['GW mult','depth','K']
for i in range(3):
    plt.plot(price_mult, ST_vol_dep[:,i], label = labels[i])
plt.legend()

# Analyze results
results_array = np.array(results_end_div_start_area)
S1_area_delt = np.zeros((5,3))
ST_area_delt = np.zeros((5,3))
for i in range(5):
    sobol_indices = sobol.analyze(problem, results_array[:,i], calc_second_order = True)
    S1_area_delt[i,:] = sobol_indices['S1']
    ST_area_delt[i,:] = sobol_indices['ST']
    
# Plots of S1s
plt.figure()
labels = ['GW mult','depth','K']
for i in range(3):
    plt.plot(price_mult, S1_area_delt[:,i], label = labels[i])
plt.legend()

# Plots of STs
plt.figure()
labels = ['GW mult','depth','K']
for i in range(3):
    plt.plot(price_mult, ST_area_delt[:,i], label = labels[i])
plt.legend()
    
# Analyze results
results_array = np.array(results_perc_vol_depleted)
S1_frac_gw_dep = np.zeros((5,3))
ST_frac_gw_dep = np.zeros((5,3))
for i in range(5):
    sobol_indices = sobol.analyze(problem, results_array[:,i], calc_second_order = True)
    S1_frac_gw_dep[i,:] = sobol_indices['S1']
    ST_frac_gw_dep[i,:] = sobol_indices['ST']
    
# Plots of S1s
plt.figure()
labels = ['GW mult','depth','K']
for i in range(3):
    plt.plot(price_mult, S1_frac_gw_dep[:,i], label = labels[i])
plt.legend()

# Plots of STs
plt.figure()
labels = ['GW mult','depth','K']
for i in range(3):
    plt.plot(price_mult,  ST_frac_gw_dep[:,i], label = labels[i])
plt.legend()

# Analyze results
results_array = np.array(results_total_profit)
S1_prof = np.zeros((5,3))
ST_prof = np.zeros((5,3))
for i in range(5):
    sobol_indices = sobol.analyze(problem, results_array[:,i], calc_second_order = True)
    S1_prof[i,:] = sobol_indices['S1']
    ST_prof[i,:] = sobol_indices['ST']
    
# Plots of S1s
plt.figure()
labels = ['GW mult','depth','K']
for i in range(3):
    plt.plot(price_mult, S1_prof[:,i], label = labels[i])
plt.legend()

# Plots of STs
plt.figure()
labels = ['GW mult','depth','K']
for i in range(3):
    plt.plot(price_mult,  ST_prof[:,i], label = labels[i])
plt.legend()

#%% Sobol plots 

##plt.plot(price_mult, S1_vol_dep[])


# range of dependant variable outcomes
prediction_interval = 95
y = results_end_div_start_area
plt.plot(price_mult, np.mean(y, axis=0), label="Mean", color='black')
plt.fill_between(price_mult,
                 np.percentile(y, 50 - prediction_interval/2., axis=0),
                 np.percentile(y, 50 + prediction_interval/2., axis=0),
                 alpha=0.5, color='black',
                 label=f"{prediction_interval} % prediction interval")

plt.xlabel("Price_mult")
plt.ylabel("End/Start Area")
plt.legend(loc=4)

plt.show()
