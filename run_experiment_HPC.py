##### Load all necessary modules #####
#import pdb
import traceback
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
from timeit import default_timer as timer
import os
import matplotlib.pyplot as plt
from itertools import product
from pathlib import Path
import sys
import logging
import warnings
warnings.filterwarnings("ignore")
logging.getLogger('pyomo.core').setLevel(logging.ERROR)
logging.getLogger('pyomo.repn').setLevel(logging.ERROR)
logging.getLogger('pyomo.repn.plugins').setLevel(logging.ERROR)
logging.getLogger('pyomo').setLevel(logging.ERROR)

##### Travis: Lines 19-92 all load external data and set variables that are used for all runs -- need to find strategy to do this only once when deployed on HPC
##### Travis: Need to modify this such that it can be deployed on HPC
#os.chdir('C:/Users/fere556/Desktop/Farm_GW_Analysis') # comment out for HPC implementation

##### Load Theis drawdown response module
import Superwell_for_ABM_on_the_fly


def pretty_timer(seconds: float) -> str:
    """Formats an elapsed time in a human friendly way.

    Args:
        seconds (float): a duration of time in seconds

    Returns:
        str: human friendly string representing the duration
    """
    if seconds < 1:
        return f'{round(seconds * 1.0e3, 0)} milliseconds'
    elif seconds < 60:
        return f'{round(seconds, 3)} seconds'
    elif seconds < 3600:
        return f'{int(round(seconds) // 60)} minutes and {int(round(seconds) % 60)} seconds'
    elif seconds < 86400:
        return f'{int(round(seconds) // 3600)} hours, {int((round(seconds) % 3600) // 60)} minutes, and {int(round(seconds) % 60)} seconds'
    else:
        return f'{int(round(seconds) // 86400)} days, {int((round(seconds) % 86400) // 3600)} hours, and {int((round(seconds) % 3600) // 60)} minutes'


def solve_farm(fid: int, output_path):
    
    print(f'Solving for farm {fid}...')

    time_start = timer()

    ##### Adjust pandas setting for debugging
    pd.set_option('display.max_columns', None)
    pd.set_option('display.max_rows', None)

    ##### Initiate list of all farm ids across CONUS
    farm_ids = range(53835)  # total number of farm agents (nldas IDs, 1 farm per NLDAS grid cell)

    ##### Load in master farms table (note: this takes a couple minutes as we are loading a large table over all crops/farms,
    ##### Travis: This is a large file that loads the data for all farms across the CONUS -- for ensemble experiment, we'd only want to load it once.
    farms_master = pd.ExcelFile("./data_inputs/PMP_inputs_sensitivity_20220829.xlsx")
    farms_master = farms_master.parse("crop_level")
    baseline_profit_per_acre = pd.read_csv("baseline_profit_analysis.csv", index_col = 0) # Stephen added

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

    ##### For each representative farm, estimate the starting GW demand based on GW-irrigated crop areas and net irrigation requirements (nir)
    farms_master['gw_irr_vol'] = farms_master['area_irrigated_gw'] * farms_master['nir_corrected'] * 1000
    aggregation_functions = {'gw_irr_vol': 'sum'}
    farm_gw_sum = farms_master[['nldas', 'gw_irr_vol']].groupby(['nldas']).aggregate(aggregation_functions)
    farm_gw_sum = farm_gw_sum.reset_index()
    
    # Stephen added baseline GW irrigated area calculation 
    aggregation_functions = {'area_irrigated_gw': 'sum'}
    farm_gw_area_sum = farms_master[['nldas', 'area_irrigated_gw']].groupby(['nldas']).aggregate(aggregation_functions)
    farm_gw_area_sum = farm_gw_area_sum.reset_index()
    
    # farm_gw_sum['no_of_wells'] = farm_gw_sum['area_irrigated_gw'] / (gw_area_well_acres)  # .000247105 conversion from square meters to acres

    time_load = time.time()

    # Universal cost curve inputs. Comment our for HPC implementation
    gw_curve_option = 'internal'  # GW cost curve option ("external" - pre-processed and pulled from external file, "internal" - generate on the fly)
    ELECTRICITY_RATE = .12 # Stephen - decide what electricty rate to use
    profit_threshold = 100 # profit/acre irrigated threshold for optimization constraint 
    gw_cost_curve_internal = None

    ##### Load external files for NLDAS cost curve attributes
    nldas_gw_attributes = pd.read_csv('NLDAS_Cost_Curve_Attributes.csv')

    # Scenario values For debugging 
    # hydro_ratio = 1 
    # econ_ratio = 1
    # gamma_scenario = 1 
    # farm_id = 7832
    # K_scenario = 'sigma_high'
    
    ##### Travis: Lines 130-694 define the method for running the farm-groundwater model
    # Function-based workflow for incorporation in various sensitivity/ensemble experiments (e.g., SALib Sobol Analysis)
    def farm_gw_model(hydro_ratio, econ_ratio, K_scenario, gamma_scenario, farm_id, gw_cost_curve_internal):
        
        gw_cost_curve_internal = None
        
        # ratios used for future exploratory scenario analysis
        percent_increment = (hydro_ratio - 1.0) / (hydro_ratio + 1.0)
        numer_hydro_factor = 1 + percent_increment
        denom_hydro_factor = 1 - percent_increment
        percent_increment = (econ_ratio - 1.0) / (econ_ratio + 1.0)
        numer_econ_factor = 1 + percent_increment
        denom_econ_factor = 1 - percent_increment

        farm_id = [farm_id]  # Pick the farm agent [0-53834] to use as a basis for the sensitivity analysis
        no_of_years = 100  # The number of years (annual timesteps) to run simulation, can be toggled

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

        # time_load = time.time()
        
        #### Stephen - modified this code block to load nldas aquifer properties for S, m, K, and WL
        if gw_curve_option == 'internal':

            nldas_id = farms_master.iloc[farm_id]['nldas']
            R = nldas_gw_attributes[(nldas_gw_attributes['NLDAS_ID'] == nldas_id.values[0])]['Recharge_usgs']
            R = R * numer_hydro_factor
            m_value = nldas_gw_attributes[(nldas_gw_attributes['NLDAS_ID'] == nldas_id.values[0])]['Depth'] # Stephen - assign from nldas grid properties 
            WL_value = nldas_gw_attributes[(nldas_gw_attributes['NLDAS_ID'] == nldas_id.values[0])]['DTW'] # Stephen - assign from nldas grid properties
            S_value = nldas_gw_attributes[(nldas_gw_attributes['NLDAS_ID'] == nldas_id.values[0])]['Porosity'] # Stephen - assign from nldas grid properties
            IRR_DEPTH_afy = farm_gw_sum.iloc[farm_id[0]]['gw_irr_vol']
            # IRR_DEPTH_afy = farms_master.iloc[farm_id]['area_irrigated_gw'] * farms_master.iloc[farm_id]['nir_corrected'] * 1000
            # IRR_DEPTH_m3 = IRR_DEPTH_afy * 1233.48  # convert from acre-ft/yr to cubic meters/yr
            
            # UPDATED - IRR_DEPTH CALCULATED BASED ON BASALINE GW CROP AREA NOT TOTAL GRID CELL AREA. USING TOTAL GRID AREA (above lines commented out) RESULTS
            # IN EXTREMELY SMALL IRR_DEPTHS FOR CELLS WITH SMALL BASELINE GW USE. USE: farm_gw_area_sum instead of nldas_id.values[0])]['Area']

            # IRR_DEPTH = IRR_DEPTH_m3.values[0]/(nldas_gw_attributes[(nldas_gw_attributes['NLDAS_ID'] == nldas_id.values[0])]['Area'].values[0]) # convert to m/year, Stephen - conversion is correct
            IRR_DEPTH = 1.25 * IRR_DEPTH_afy/farm_gw_area_sum.iloc[farm_id[0]]['area_irrigated_gw'] # IRR_DEPTH in feet with a 1.25 safety factor, updated IRR_DEPTH CALCULATED BASED ON BASALINE GW CROP AREA NOT TOTAL GRID CELL AREA

            # Assign K_value from K_scenario
            # Step 1: Lookup K and calculate K standard deviation using K and K_hi values from NLDAS dataset 
            K_default = nldas_gw_attributes[(nldas_gw_attributes['NLDAS_ID'] == nldas_id.values[0])]['K']
            K_high = nldas_gw_attributes[(nldas_gw_attributes['NLDAS_ID'] == nldas_id.values[0])]['K_hi']
            K_sigma_log = np.log10(K_high) - np.log10(K_default) # calculate K standard deviation 
            
            # Step 2: Calculate list of plausible K values: default K +/- 0.5 and 1 sigma
            K_range = [10**(np.log10(K_default) - K_sigma_log).values[0], 10**(np.log10(K_default) - 0.5* K_sigma_log).values[0],
                       K_default.values[0], 10**(np.log10(K_default) + 0.5 * K_sigma_log).values[0], 10**(np.log10(K_default) + K_sigma_log).values[0]]
            
            K_value = K_range[K_scenario_dict[K_scenario]] # Stephen - assign based on K_scenario, note = K below 0.1 m/d mostly non-viable
            
            # Generate Cost Curve with Superwell_for_ABM_on_the_fly module 
            gw_cost_curve_internal = Superwell_for_ABM_on_the_fly.cost_curve(S_value.values[0], m_value.values[0], K_value, WL_value.values[0], R.values[0], IRR_DEPTH * 0.3048, no_of_years, ELECTRICITY_RATE) 
                                    
        # if the Superwell_for_ABM_on_the_fly fails to generate a cost curve (returns None), return all nan values
        if not gw_cost_curve_internal:
            return [np.nan] * 60 # combined number of outputs for Total Simulation and Annual Outputs 
        
        # Stephen - Post-process unit cost and volume pumped in case of any nan or inf values 
        # fill in Nan or Inf values for annual volume pumped and unit cost
        for n in range(no_of_years):
            
            # Replace Inf value with 0
            if gw_cost_curve_internal[0][n] > 100000: 
                gw_cost_curve_internal[0][n] = 0 
                gw_cost_curve_internal[1][n] = 0 
                
            # Fill Nan value by linear interpolation
            if np.isnan(gw_cost_curve_internal[0][n]) == 1:
                
                if gw_cost_curve_internal[0][n+1] == 0:
                    gw_cost_curve_internal[0][n] = 0
                    gw_cost_curve_internal[1][n] = 0  
                    
                elif n < no_of_years-2:
                    slope = (gw_cost_curve_internal[0][n-1] - gw_cost_curve_internal[0][n-1])/2
                    gw_cost_curve_internal[0][n] = gw_cost_curve_internal[0][n-1] + slope
                        
                else:
                    gw_cost_curve_internal[0][n] = 0 
                    gw_cost_curve_internal[1][n] = 0 
                        
        # apply sensitivity multipliers
        gammas_total_subset = gammas_total_subset_og.copy()
        for key in gammas_total_subset:
            gammas_total_subset[key] *= gamma_scenario
            
        # alphas_total_subset = alphas_total_subset_og.copy() # Stephen - if we want to add Alpha to parameter sensitivity
        # for key in alphas_total_subset:
        #     alphas_total_subset[key] *= alpha_scenario

        nirs_subset = nirs_subset_og.copy()
        for key in nirs_subset:
            nirs_subset[key] *= denom_hydro_factor

        sw_calib_constr_subset = sw_calib_constr_subset_og.copy()
        for key in sw_calib_constr_subset:
            sw_calib_constr_subset[key] *= numer_hydro_factor

        net_prices_land_subset = net_prices_land_subset_og.copy()
        for key in net_prices_land_subset:
            net_prices_land_subset[key] = ((yields_subset[key] * prices_subset[key] * numer_econ_factor) - land_costs_subset[key] -
                                           alphas_total_subset[key]) # Stephen - check sign of Alpha term

        # set land net prices to large negative # for gammas that are zero
        # (so PMP won't produce crops, zero gamma value indicates no observed land use)
        for key, value in gammas_total_subset_og.items():
            if value == 0:
                net_prices_land_subset[key] = -9999999999

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
        
        # MODIFICATION: CHECK IF SW COST > GW COST & IF SW AREA > 10%. IF 
        # TRUE, THE GW CONSTAINT SET AT THE BASELINE. IF FALSE, THE GW CONSTRAINT 
        # WILL BE ALLOWED TO INCREASE AT A DEFINED RATE 
        if net_prices_sw_subset[farm_id[0]] < net_prices_gw_subset[farm_id[0]] and \
           sw_calib_constr_subset[farm_id[0]]/(sw_calib_constr_subset[farm_id[0]] + \
           gw_calib_constr_subset[farm_id[0]]) > 0.01: # gw & sw area or volume?
            gw_expansion = False 
            
        else:
            gw_expansion = True 
            
        # MODIFICATION: CHECK BASELINE PROFIT/ACRE. FLAG IF BASELINE IS LOWER 
        # THAN PROFIT/ACRE THRESHOLD. FLAG WILL SKIP THE PROFIT/AREA CONSTRAINT
        profit_constr = True
        if baseline_profit_per_acre[(baseline_profit_per_acre['NLDAS'] == nldas_id.values[0])]['profit_per_acre'].values[0] < profit_threshold:
            profit_constr = False 
        
        # start simulation time loop
        for t in range(no_of_years):

            # initialize pyomo model
            fwm_s = ConcreteModel()
            fwm_s.ids = Set(initialize=ids_subset_sorted)
            
            # allows for evolution of GW cost and availability
            if first:
                fwm_s.net_prices_gw = Param(fwm_s.ids, initialize=net_prices_gw_sensitivity, mutable=True)
            else:
                net_prices_gw_subset_update = net_prices_gw_sensitivity.copy()

                for key in net_prices_gw_subset_update:
                    net_prices_gw_temp = net_prices_gw_subset_update[key]
                    nir_temp = nirs_subset[key]
                    gw_cost_temp = (-1.0 * net_prices_gw_temp) / (nir_temp * 1000.0) # GW cost $/acre (same value for all 10 crop types)
                    gw_cost_updated_temp = gw_cost_temp + (gw_cost_added * 1233.48) # factor of 1233.48 converts $/m3 from cost curves to $/acre-ft used in PMP
                    net_prices_gw_updated_temp = -1.0 * gw_cost_updated_temp * nir_temp * 1000.0 # converts increase GW cost ($/acre-ft) to $/acre using crop-specific nir 
                    net_prices_gw_subset_update[key] = net_prices_gw_updated_temp # updates net GW price for each crop 
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
            fwm_s.nirs = Param(fwm_s.ids, initialize=nirs_subset, mutable=True)
            
            ########## Added by Stephen to apply reduced GW availability multiplier to the gw_constraint
            gw_calib_constr_subset_copy = gw_calib_constr_subset.copy()
            
            # MODIFICATION: THIS IS WHERE THE INCREASED ANNUAL GW VOLUME CONSTRAINT IS UPDATED. 
            # THIS INVOLVES MULTIPLYING THE GW CONSTRAINT BY 10% PERCENTAGE OF BASELINE OR ADDING 500 ACFT/YEAR, 
            # WHICHEVER IS GREATE. THE MAXIMUM EXPANISION OF GW CAPACITY IS CAPPED BY INITIAL NIR*CELL AREA
            
            if gw_expansion == False:
                gw_calib_constr_subset_copy[farm_id[0]] = gw_calib_constr_subset_copy[farm_id[0]] * gw_availability_multiplier
                fwm_s.gw_constraints = Param(fwm_s.farm_ids, initialize=gw_calib_constr_subset_copy, mutable=True)
                
            else: # If expansion = True: Annual GW expansion max of 1) 10% of baseline OR 2) 500 acre-feet
                if t == 0:
                    gw_constr = gw_calib_constr_subset_copy
                    gw_capacity_increase = max(0.1 * gw_constr[farm_id[0]], 500)
                    gw_constr[farm_id[0]] = gw_constr[farm_id[0]] + gw_capacity_increase 
                    
                elif gw_constr[farm_id[0]] < IRR_DEPTH * max_land_constr_subset[farm_id[0]]/1000 and gw_availability_multiplier == 1: # limit area to land constraint, use baseline IRR_DEPTH defined by nir 
                    if gw_constr[farm_id[0]] + gw_capacity_increase > IRR_DEPTH * max_land_constr_subset[farm_id[0]]/1000:
                        gw_constr[farm_id[0]] =  IRR_DEPTH * max_land_constr_subset[farm_id[0]]/1000
                    else:
                        gw_constr[farm_id[0]] = gw_constr[farm_id[0]] + gw_capacity_increase
                
                else:
                    gw_constr[farm_id[0]] = gw_constr[farm_id[0]] * gw_availability_multiplier
                    
                fwm_s.gw_constraints = Param(fwm_s.farm_ids, initialize = gw_constr, mutable=True)
 
            # Intialized crop area at 0 
            x_start_values = dict(enumerate([0.0] * 3))
            fwm_s.xs_total = Var(fwm_s.ids, domain=NonNegativeReals, initialize=x_start_values[0])
            fwm_s.xs_sw = Var(fwm_s.ids, domain=NonNegativeReals, initialize=x_start_values[1])
            fwm_s.xs_gw = Var(fwm_s.ids, domain=NonNegativeReals, initialize=x_start_values[2])
            
            
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
                return sum(fwm_s.xs_total[i] for i in fwm_s.crop_ids_by_farm_and_constraint[ff]) <= fwm_s.land_constraints[ff]

            fwm_s.c1 = Constraint(fwm_s.farm_ids, rule=land_constraint)

            def obs_lu_constraint_sum(fwm_s, i):
                return fwm_s.xs_sw[i] + fwm_s.xs_gw[i] == fwm_s.xs_total[i]

            fwm_s.c5 = Constraint(fwm_s.ids, rule=obs_lu_constraint_sum)

            def sw_constraint(fwm_s, ff):
                return sum(fwm_s.xs_sw[i] * fwm_s.nirs[i] for i in fwm_s.crop_ids_by_farm_and_constraint[ff]) <= fwm_s.sw_constraints[ff]

            fwm_s.c2 = Constraint(fwm_s.farm_ids, rule=sw_constraint)

            def gw_constraint(fwm_s, ff):
                return sum(fwm_s.xs_gw[i] * fwm_s.nirs[i] for i in fwm_s.crop_ids_by_farm_and_constraint[ff]) <= fwm_s.gw_constraints[ff]

            fwm_s.c6 = Constraint(fwm_s.farm_ids, rule=gw_constraint)
                        
            # MODIFICATION: PROFIT/ACRE CONSTRAINT       
            if profit_constr == True: 
                def profit_constraint(fwm_s, ff):
                    return sum(sum((fwm_s.net_prices_total[h] * fwm_s.xs_total[h]/1000 - 0.5 * fwm_s.gammas_total[h] *
                                              fwm_s.xs_total[h]/1000 * fwm_s.xs_total[h]/1000) for h in fwm_s.crop_ids_by_farm[f]) +
                                          sum((fwm_s.net_prices_sw[i] * fwm_s.xs_sw[i]/1000) for i in fwm_s.crop_ids_by_farm[f]) +
                                          sum((fwm_s.net_prices_gw[g] * fwm_s.xs_gw[g]/1000) for g in fwm_s.crop_ids_by_farm[f]) for f
                                          in fwm_s.farm_ids)/sum(fwm_s.xs_total[i]/1000 for i in fwm_s.crop_ids_by_farm_and_constraint[ff]) >= profit_threshold
                    
    
                fwm_s.c7 = Constraint(fwm_s.farm_ids, rule=profit_constraint)
                
            fwm_s.scaling_factor = Suffix(direction=Suffix.EXPORT)
            fwm_s.scaling_factor[fwm_s.xs_total] = 0.0001
            fwm_s.scaling_factor[fwm_s.xs_sw] = 0.0001
            fwm_s.scaling_factor[fwm_s.land_constraints] = 0.000001
            fwm_s.scaling_factor[fwm_s.sw_constraints] = 0.01
            fwm_s.scaling_factor[fwm_s.gw_constraints] = 0.01
            
            fwm_s.scaling_factor[fwm_s.c1] = 0.000001 # land_constraint
            fwm_s.scaling_factor[fwm_s.c5] = 0.000001 # obs_lu_constraint 
            fwm_s.scaling_factor[fwm_s.c2] = 0.01     # sw_constraint 
            fwm_s.scaling_factor[fwm_s.c6] = 0.01     # gw_constraint
            if profit_constr == True:
                fwm_s.scaling_factor[fwm_s.c7] = 0.01   # Stephen - choose scaling factor for profit constraint, tested 0.01, 1, 10 and no difference

            
            # create and run the optimization solver
            opt = SolverFactory("ipopt", solver_io='nl') # USE this for HPC 
            
            #opt = SolverFactory('gurobi', solver_io='python') # Stephen used for testing 
            #opt = SolverFactory('appsi_highs')#, solver_io='python')
            #opt = SolverFactory('asl:highs')#, solver_io='python')
            #opt.options["presolve"] = "on"
            #opt.options["parallel"] = "on"
            #opt.options["solver"] = "simplex"
            #opt.options["simplex_strategy"] = 1
            #pdb.set_trace()
            #try:
            results = opt.solve(fwm_s, keepfiles=False, tee=False)
            #print(results.solver.termination_condition)

            # except ValueError:
            #   continue

            # store main model outputs
            result_xs_sw = dict(fwm_s.xs_sw.get_values())
            result_xs_gw = dict(fwm_s.xs_gw.get_values())
            result_xs_total = dict(fwm_s.xs_total.get_values())

            farms_list = [farms_master.loc[[farm_id[0]]].nldas.values[0]]

            # process groundwater production results and ? 
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
            results_pd['gw_vol'] = results_pd['xs_gw'] * results_pd['nir_corrected'] * denom_hydro_factor  # Stephen - area * nir (depth) = volume, need to correct nir by factor of 1,000? Looks like the factor of 1000 for nir is accoutned for because xs_gw is increased by 1000
            
            # Create variables to track gw use, change in cost, and change in availability (annual capacity)
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


            # Create and appends annual results
            if first:
                results_combined = results_pd
            else:
                results_combined = pd.concat([results_combined, results_pd])

            aggregation_functions = {'gw_vol': 'sum', 'xs_gw': 'sum'}
            gw_sum = results_pd[['nldas', 'gw_vol', 'xs_gw']].groupby(['nldas']).aggregate(aggregation_functions)
            gw_sum = gw_sum.reset_index()
            # gw_sum = gw_sum.merge(farm_gw_sum[['nldas', 'no_of_wells']], how='left', on=['nldas'])
            # gw_sum['gw_vol_well'] = gw_sum['gw_vol'] * farms_per_grid / gw_sum['no_of_wells']
            gw_sum['gw_vol_km3'] = gw_sum['gw_vol'] * 1.23348e-6 # Stephen - check unit conversion for gw_sum (check conversion from acre-feet to km^3 using 1.23348e-6 multiplier)
            
            if first:
                GW_pumped_timeseries = [] # Stephen - added empty array track annual GW pumping
            
            GW_pumped_timeseries.append(gw_sum['gw_vol_km3'].values[0]) # Stephen - added empty array track annual GW pumping

            # Update GW price at end of each time step 
            if first:
                results_combined['net_prices_gw_updated'] = -1.0 * (results_combined['gw_cost'] + (gw_cost_added * 1233.48)) * results_combined['nir_corrected'] * 1000.0 # Stephen - seems like this should be inside the annual simulation code block since gw_cost_added is a single-value that gets overwritten 
            else:
                results_combined.net_prices_gw_updated[t*10:(t+1)*10+1] = -1.0 * (results_combined.gw_cost.iloc[0:10] + (gw_cost_added * 1233.48)) * results_combined.nir_corrected.iloc[0:10] * 1000.0 # Stephen - seems like this should be inside the annual simulation code block since gw_cost_added is a single-value that gets overwritten 


            if first:
                cumul_gw_sum = gw_sum['gw_vol_km3'].values[0]
            else:
                cumul_gw_sum += gw_sum['gw_vol_km3'].values[0]
            
            i = 0 # index to track progression along cost curve, i set to 0 before entering code block for
                  # each annual time step (t)   
                  
            if first:
                WL_timeseries = [] # Stephen - added empty array to track annual changes in WL due to depletion 
                GW_cost_timeseries = [] # Stephen - added empty array to track annual GW cost 
                GW_cost_increase = [] # Stephen - added empty array to track total increase in annual GW cost 
                
            if gw_curve_option == 'internal':
                for index in range(gw_cost_curve_internal[0].size):
                    if cumul_gw_sum > gw_cost_curve_internal[1][index]:  # closest match for cum pumped vol cost curve bin 
                        i = index
                    else:
                        break
                    
                if gw_cost_curve_internal[0][i] != 0: # if there is a unit cost at index i

                    if i == 0:
                        # absolute change in unit cost: current unit cost minus starting unit cost
                        gw_cost_added = 0
                        gw_cost = gw_cost_curve_internal[0][0]

                        # fractional change in unit cost: current unit cost divided by starting unit cost)
                        #gw_multiplier = 1

                        # change in annual pumping capacity due to reduce well pumping rates 
                        # if option == 'reduced_capacity':  # logic for optional reduced capacity 
                        gw_availability_multiplier = 1 # Stephen - do we want reduced capacity to be optional, right now it is the default 

                        WL_value = gw_cost_curve_internal[2][0]

                    else:
                        # interpolate GW unit cost
                        slope_incremental_cost = (gw_cost_curve_internal[0][i+1] - gw_cost_curve_internal[0][i])/(gw_cost_curve_internal[1][i+1] - gw_cost_curve_internal[1][i])  # new implementation with interpolation
                        dx = cumul_gw_sum - gw_cost_curve_internal[0][i-1] # excess volume of GW pumped above the previous cost curve pumped volume index
                        gw_cost_increase = dx * slope_incremental_cost # incremental cost increase compared to previous cost curve pumped volume index 
                        gw_cost = gw_cost_curve_internal[0][i] + gw_cost_increase # GW cost for current time step
                        gw_cost_added = gw_cost - gw_cost_curve_internal[0][0] # additional GW cost compared to initial GW cost for current time step

                        # interpolate WL
                        slope_incremental_WL = (gw_cost_curve_internal[2][i+1] - gw_cost_curve_internal[2][i])/(gw_cost_curve_internal[1][i+1] - gw_cost_curve_internal[1][i]) # new
                        wl_decline = dx * slope_incremental_WL # incremental increase in WL depth compared to previous cost curve pumped volume index 
                        WL_value =  gw_cost_curve_internal[2][i] + wl_decline # WL for current time step

                        # change in annual pumping capacity due to reduce well pumping rates
                        # if option == 'reduced_capacity':  # logic for optional reduced capacity
                        gw_availability_multiplier = gw_cost_curve_internal[3][i] # Stephen - do we want reduced capacity to be optional, right now it is the default

                    # update WL - Stephen added 
                    #WL_timeseries.append(gw_cost_curve_internal[2][i]) # current
                    WL_timeseries.append(WL_value) # new

                    # update GW cost - Stephen added 
                    #GW_cost_timeseries.append(gw_cost_curve_internal[0][i]) #current  - UPDATE to initial GW cost from Farm database + gw_cost_added 
                    GW_cost_timeseries.append(gw_cost_curve_internal[0][i-1] + gw_cost) # new
                    
                    # update GW added cost - Stephen added 
                    GW_cost_increase.append(gw_cost_added)

                else:
                    gw_cost_added = 9999999999999
                    gw_multiplier = 9999999999999
                    
                    # update WL - Stephen added 
                    WL_timeseries.append(max(gw_cost_curve_internal[2]))

                    # update GW cost - Stephen added - what cost to assign to output timeseries when no GW is available? 999 for now 
                    GW_cost_timeseries.append(999)
                    
                    # update GW added cost - Stephen added 
                    GW_cost_increase.append(999)
                    
                    if depletion_first:
                        time_depletion = t
                        depletion_first = False
            
            
            first = False

        # Calculate profit - calculates profit for each crop at each annual time step # check signs on alphas land, net prices gw and sw
        results_combined['profit_calc'] = ((numer_econ_factor * (results_combined['price'] * results_combined['yield']) -
                                            results_combined['land_cost'] - results_combined['alphas_land']) *
                                           results_combined['xs_total']/1000) \
                                          - (0.5 * results_combined['gammas_total'] * results_combined['xs_total']/1000 *
                                             results_combined['xs_total']/1000) \
                                          + (results_combined['net_prices_gw_updated'] * results_combined['xs_gw']/1000) \
                                          + (results_combined['net_prices_sw'] * results_combined['xs_sw']/1000)


        #### Calculate summary results (percent change in total irrigated area as an initial result)
        # In original version, these were only single-value outputs. For the updated version, 
        # annual metrics are lists of annual values 
        
        # Total Simulation or End of Simulation Metrics - single values 
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
        
        
        # Stephen - Potential new Total Simulation metrics 
        GW_cost_filtered = GW_cost_timeseries.copy()
        GW_cost_filtered = np.array(GW_cost_filtered)
        GW_cost_indexes = np.where(GW_cost_filtered < 100)
        GW_cost_filtered = GW_cost_filtered[GW_cost_indexes]
          
        Total_increase_GW_cost = max(GW_cost_filtered) - min(GW_cost_filtered)
        Final_WL = max(WL_timeseries)
        
        # Volume depleted (calculated from cost curves)
        #available_volume = (m - WL) *  grid_cell_area * S
        #perc_vol_depleted = cumul_gw_sum * (10 ** 9) / available_volume
        perc_vol_depleted = (max(WL_timeseries) - gw_cost_curve_internal[2][0])/(m_value.values[0]-gw_cost_curve_internal[2][0])
        
        if perc_vol_depleted > 1:
            perc_vol_depleted = 1
            
        # Annual Metrics - Stephen added 
        sim_years = max(results_combined.year) # years of output 
        
        Farm_id = np.tile(farm_id, sim_years).tolist()
        Hydro_ratio = np.tile(hydro_ratio, sim_years).tolist()
        Econ_ratio = np.tile(econ_ratio, sim_years).tolist()
        K_scenario_list = np.tile(K_scenario, sim_years).tolist()
        K_values = np.tile(K_value, sim_years).tolist()
        m_values = np.tile(m_value.values[0], sim_years).tolist()
        Year = np.arange(1,sim_years+1).tolist()
        GW_WL = WL_timeseries # populate using WL_timeseries
        GW_vol_pumped = GW_pumped_timeseries # populate using GW_pumped_timeseries 
        GW_fraction_depleted = ((m_values-np.array(WL_timeseries))/(m_values-min(WL_timeseries))).tolist() # caculate using WL_timeseries, aquifer depth, and initial sat thickness 
        GW_cost = GW_cost_timeseries # populate using GW_cost_timeseries 
        
        # For annual profit iterate over each year of results_combined 
        Profit = []
        for i in range(sim_years):
            annual_profit = results_combined.where(results_combined.year == i+1)
            annual_profit = annual_profit.dropna()
            Profit.append(sum(annual_profit.profit_calc))
            
        # Annual areas available from 'results_combined' 
        Area_Corn = results_combined.where(results_combined.crop == 'Corn').dropna()
        Area_GW_Corn = Area_Corn.xs_gw.tolist()
        Area_SW_Corn = Area_Corn.xs_sw.tolist()
        
        Area_FiberCrop = results_combined.where(results_combined.crop == 'FiberCrop').dropna()
        Area_GW_FiberCrop = Area_FiberCrop.xs_gw.tolist()
        Area_SW_FiberCrop = Area_FiberCrop.xs_sw.tolist()
        
        Area_FodderGrass = results_combined.where(results_combined.crop == 'FodderGrass').dropna()
        Area_GW_FodderGrass = Area_FodderGrass.xs_gw.tolist()
        Area_SW_FodderGrass = Area_FodderGrass.xs_sw.tolist()
        
        Area_MiscCrop = results_combined.where(results_combined.crop == 'MiscCrop').dropna()
        Area_GW_MiscCrop = Area_MiscCrop.xs_gw.tolist()
        Area_SW_MiscCrop = Area_MiscCrop.xs_sw.tolist()
        
        Area_OilCrop = results_combined.where(results_combined.crop == 'OilCrop').dropna()
        Area_GW_OilCrop = Area_OilCrop.xs_gw.tolist()
        Area_SW_OilCrop = Area_OilCrop.xs_sw.tolist()

        Area_OtherGrain = results_combined.where(results_combined.crop == 'OtherGrain').dropna() 
        Area_GW_OtherGrain = Area_OtherGrain.xs_gw.tolist()
        Area_SW_OtherGrain = Area_OtherGrain.xs_sw.tolist()

        Area_Rice = results_combined.where(results_combined.crop == 'Rice').dropna()
        Area_GW_Rice = Area_Rice.xs_gw.tolist()
        Area_SW_Rice = Area_Rice.xs_sw.tolist()

        Area_Root_Tuber = results_combined.where(results_combined.crop == 'Root_Tuber').dropna()      
        Area_GW_Root_Tuber = Area_Root_Tuber.xs_gw.tolist()
        Area_SW_Root_Tuber = Area_Root_Tuber.xs_sw.tolist()
        
        Area_SugarCrop = results_combined.where(results_combined.crop == 'SugarCrop').dropna()
        Area_GW_SugarCrop = Area_SugarCrop.xs_gw.tolist()
        Area_SW_SugarCrop = Area_SugarCrop.xs_sw.tolist()
        
        Area_Wheat = results_combined.where(results_combined.crop == 'Wheat').dropna()
        Area_GW_Wheat = Area_Wheat.xs_gw.tolist()
        Area_SW_Wheat = Area_Wheat.xs_sw.tolist() 
        
        # For annual total gw and sw crop area iterate over each year of results_combined 
        Area_GW_total = []
        Area_SW_total = []
        for i in range(sim_years):
            annual_profit = results_combined.where(results_combined.year == i+1)
            annual_profit = annual_profit.dropna()
            Area_GW_total.append(sum(annual_profit.xs_gw))
            Area_SW_total.append(sum(annual_profit.xs_sw))      
            
        Gamma_mult = np.tile(gamma_scenario, sim_years).tolist()
        S_values = np.tile(S_value.values[0], sim_years).tolist()
        
        
        # Stephen - Added additional output columns
        return [end_div_start_area, total_profit, perc_vol_depleted, time_depletion, cumul_gw_sum, acres_grown_total,
                acres_grown_mean, acres_grown_mean_vs_start, end_div_start_profit, min_acres_grown, max_acres_grown,
                max_minus_min_acres_grown, med_acres_grown, min_annual_profit, max_annual_profit,
                max_minus_min_annual_profit, med_annual_profit, max_gw_vol_yr, min_gw_vol_yr, max_minus_min_gw_vol_yr,
                med_gw_vol_yr, mean_excl0_gw_vol_yr, Final_WL, Total_increase_GW_cost, 
                
                Farm_id, Hydro_ratio, Econ_ratio, K_scenario_list, K_values, m_values, Year, 
                GW_WL, GW_vol_pumped, GW_fraction_depleted, GW_cost, 
                Profit, Area_GW_Corn, Area_SW_Corn, Area_GW_FiberCrop, Area_SW_FiberCrop,
                Area_GW_FodderGrass, Area_SW_FodderGrass, Area_GW_MiscCrop, Area_SW_MiscCrop,
                Area_GW_OilCrop, Area_SW_OilCrop, Area_GW_OtherGrain, Area_SW_OtherGrain, 
                Area_GW_Rice, Area_SW_Rice, Area_GW_Root_Tuber, Area_SW_Root_Tuber, 
                Area_GW_SugarCrop, Area_SW_SugarCrop, Area_GW_Wheat, Area_SW_Wheat,
                Area_GW_total, Area_SW_total, Gamma_mult, S_values]
                

    ##### Travis: Lines 468-523 are used to define the ensemble information, which is consolidated in the "cases_df" pandas dataframe
    # Varied parameters for ABM sensitivity
    hydro_ratio = [0.7, 1, 1.3] #[0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5]
    econ_ratio = [0.7, 1, 1.3] #[0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5]
    K_scenario = ['default', 'sigma_high'] #['sigma_low', 'half_sigma_low', 'default', 'half_sigma_high', 'sigma_high'] # Stephen added K scenarios & replaced K_val variable 
    K_scenario_dict = dict(sigma_low = 0, half_sigma_low = 1, default = 2, half_sigma_high = 3, sigma_high = 4)
    gamma_scenario = [0.5, 0.75, 1, 1.25, 1.5] #[0.5, 1, 1.5] #[0.8, 0.9, 1, 1.1, 1.2] # Stephen added gamma multiplier scenario values 
    # m_val = [m] # default in one m from single cost curve as defined in 'Cost curve inputs'
    single_cost_curve = 'false'

    ##### Travis - Here is where I provide a list of farms by ID that we subsequently loop through. For the comprehensive farm sweep on HPC, the idea would be to run for every farm in our database (0-53834) and strategize a way to parallelize this on HPC (each call of farm_gw_model is independent)
    farms = [fid]
    

    # Create array of ABM sensitivity parameter combinations using parameter lists defined in 'Varied parameters for ABM sensitivity'
    combinations = len(hydro_ratio) * len(econ_ratio) * len(K_scenario) * len(gamma_scenario) * len(farms)
    cases = list(product(hydro_ratio, econ_ratio, K_scenario, gamma_scenario, farms))
    cases_df = pd.DataFrame(data = cases.copy(), columns = ['hydro_ratio', 'econ_ratio', 'K_scenario', 'gamma_scenario', 'farm']) 
                    
    # Stephen - create empty DataFrame for Annual outputs, should not pre-allocate 
    Annual_df = pd.DataFrame()

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
    Final_WL = []
    Total_increase_GW_cost = []

    # Stephen - Annual metrics 
    Farm_id = []
    Hydro_ratio = [] 
    Econ_ratio = []
    K_scenario = []
    K_value = []
    m_value = [] 
    Year = []
    GW_WL = []
    GW_vol_pumped = []
    GW_fraction_depleted = []
    GW_cost = []
    Profit = []
    Area_GW_Corn = []
    Area_SW_Corn = []
    Area_GW_FiberCrop = [] 
    Area_SW_FiberCrop = []
    Area_GW_FodderGrass = [] 
    Area_SW_FodderGrass = []
    Area_GW_MiscCrop = []
    Area_SW_MiscCrop = []
    Area_GW_OilCrop = []
    Area_SW_OilCrop = []
    Area_GW_OtherGrain = [] 
    Area_SW_OtherGrain = []
    Area_GW_Rice = []
    Area_SW_Rice = []
    Area_GW_Root_Tuber = []
    Area_SW_Root_Tuber = []
    Area_GW_SugarCrop = []
    Area_SW_SugarCrop = []
    Area_GW_Wheat = []
    Area_SW_Wheat = []
    Area_GW_total = []
    Area_SW_total = []
    Gamma_mult  = []
    S_value = []
    
    ### Travis - this is where we start the simulations, looping through every case in the cases_df dataframe
    ### Simulate parameter combinations over all selected farms
    for case in range(combinations):
        #print('!!JY!! run #: ' + str(case))
        #print(f'  Case {case+1}', end="; ")
        
        # Run case and save results to Output arrays 
        results = farm_gw_model(cases_df.hydro_ratio[case], cases_df.econ_ratio[case],
                                cases_df.K_scenario[case], cases_df.gamma_scenario[case],
                                cases_df.farm[case], gw_cost_curve_internal)

        # Output arrays for Total Simulation metrics 
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
        Final_WL.append(results[22])
        Total_increase_GW_cost.append(results[23])
        
        # Stephen - Output arrays for Annual values 
        try:
            Farm_id.extend(results[24])
            Hydro_ratio.extend(results[25]) 
            Econ_ratio.extend(results[26])
            K_scenario.extend(results[27])
            K_value.extend(results[28])
            m_value.extend(results[29])
            Year.extend(results[30])
            GW_WL.extend(results[31])
            GW_vol_pumped.extend(results[32])
            GW_fraction_depleted.extend(results[33])
            GW_cost.extend(results[34])
            Profit.extend(results[35])
            Area_GW_Corn.extend(results[36])
            Area_SW_Corn.extend(results[37])
            Area_GW_FiberCrop.extend(results[38]) 
            Area_SW_FiberCrop.extend(results[39])
            Area_GW_FodderGrass.extend(results[40]) 
            Area_SW_FodderGrass.extend(results[41])
            Area_GW_MiscCrop.extend(results[42])
            Area_SW_MiscCrop.extend(results[43])
            Area_GW_OilCrop.extend(results[44])
            Area_SW_OilCrop.extend(results[45])
            Area_GW_OtherGrain.extend(results[46]) 
            Area_SW_OtherGrain.extend(results[47])
            Area_GW_Rice.extend(results[48])
            Area_SW_Rice.extend(results[49])
            Area_GW_Root_Tuber.extend(results[50])
            Area_SW_Root_Tuber.extend(results[51])
            Area_GW_SugarCrop.extend(results[52])
            Area_SW_SugarCrop.extend(results[53])
            Area_GW_Wheat.extend(results[54])
            Area_SW_Wheat.extend(results[55])
            Area_GW_total.extend(results[56])
            Area_SW_total.extend(results[57])
            Gamma_mult.extend(results[58])
            S_value.extend(results[59])
            
        except TypeError:
            continue 
        

    ##### Travis: Save the results as part of the cases_df DataFrame
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
    cases_df['final_WL'] = Final_WL
    cases_df['Total_increase_GW_cost'] = Total_increase_GW_cost

    # Stephen Added - Save results from the annual outputs to a DataFrame
    Annual_df['Farm_id'] = Farm_id
    Annual_df['Hydro_ratio'] = Hydro_ratio
    Annual_df['Econ_ratio'] = Econ_ratio
    Annual_df['K_scenario'] = K_scenario
    Annual_df['K_value'] = K_value
    Annual_df['m_value'] = m_value
    Annual_df['Year'] = Year
    Annual_df['GW_WL'] = GW_WL
    Annual_df['GW_vol_pumped'] = GW_vol_pumped
    Annual_df['GW_fraction_depleted'] = GW_fraction_depleted
    Annual_df['GW_cost'] = GW_cost
    Annual_df['Profit'] = Profit
    Annual_df['Area_GW_Corn'] = Area_GW_Corn
    Annual_df['Area_SW_Corn'] = Area_SW_Corn
    Annual_df['Area_GW_FiberCrop'] = Area_GW_FiberCrop
    Annual_df['Area_SW_FiberCrop'] = Area_SW_FiberCrop
    Annual_df['Area_GW_FodderGrass'] = Area_GW_FodderGrass
    Annual_df['Area_SW_FodderGrass'] = Area_SW_FodderGrass
    Annual_df['Area_GW_MiscCrop'] = Area_GW_MiscCrop
    Annual_df['Area_SW_MiscCrop'] = Area_SW_MiscCrop
    Annual_df['Area_GW_OilCrop'] = Area_GW_OilCrop
    Annual_df['Area_SW_OilCrop'] = Area_SW_OilCrop
    Annual_df['Area_GW_OtherGrain'] = Area_GW_OtherGrain
    Annual_df['Area_SW_OtherGrain'] = Area_SW_OtherGrain
    Annual_df['Area_GW_Rice'] =  Area_GW_Rice
    Annual_df['Area_SW_Rice'] =  Area_SW_Rice
    Annual_df['Area_GW_Root_Tuber'] = Area_GW_Root_Tuber
    Annual_df['Area_SW_Root_Tuber'] =  Area_SW_Root_Tuber
    Annual_df['Area_GW_SugarCrop'] =  Area_GW_SugarCrop
    Annual_df['Area_SW_SugarCrop'] =  Area_SW_SugarCrop
    Annual_df['Area_GW_Wheat'] =  Area_GW_Wheat
    Annual_df['Area_SW_Wheat'] =  Area_SW_Wheat
    Annual_df['Area_GW_total'] =  Area_GW_total
    Annual_df['Area_SW_total'] =  Area_SW_total
    Annual_df['Gamma_mult'] =  Gamma_mult
    Annual_df['S_value'] =  S_value
    
    ##### Travis: And export to csv. For each farm (~50k), we will have 121 scenarios, and for each farm/scenario ~30 results that are stored (~180M values)
    cases_df.to_csv(f'{output_path}/farm_{fid}_cases.csv')
    Annual_df.to_csv(f'{output_path}/farm_{fid}_annual.csv') # Stephen added 

    time_to_solve = timer() - time_start
    print(f'Farm {fid} solved in {pretty_timer(time_to_solve)}.')


if __name__ == '__main__':
    counter = int(sys.argv[3])
    output_path = Path(sys.argv[2])
    output_path.mkdir(parents=True, exist_ok=True)
    fid = pd.read_csv(sys.argv[1]).iloc[counter]['fid']
    try:
        solve_farm(fid, output_path)
    except Exception as e:
        out_path = Path(f'{output_path}/farm_{fid}_error.txt')
        with out_path.open('w') as fp:
            print(f'Farm {fid} failed to solve.')
            print(e, file=fp)
            print(traceback.format_exc())


#plt.plot(Annual_df.Profit[0:99]) 
#plt.ylabel('Annual Profit')
#plt.xlabel("year")
