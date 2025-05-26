##### Load all necessary modules #####
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


def t_depletion(depletion_ts, no_of_years, depletion_val):
    for t in range(no_of_years):
        if depletion_ts[t] > depletion_val:
            t_depletion = t + 1
            break

        else: # -999 flag if depletion threshold is not crossed
            t_depletion = -999 
            
    return t_depletion


#####  load external data and set variables that are used for all runs 

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
    with open('./data_inputs/20250214_sw_calib_constraints.p', 'rb') as fp:  # dictionary that includes sw constraint (used during calibration) keyed on farm id [0-53834]
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
    ELECTRICITY_RATE = .125 # electricty rate $/kWhr 
    profit_threshold = 100 # minimum $ profit/acre threshold used for optimization constraint 
    gw_cost_curve_internal = None

    ##### Load external files for NLDAS cost curve attributes
    nldas_gw_attributes = pd.read_csv('NLDAS_Cost_Curve_Attributes.csv')

    # Scenario values for debugging 
    # hydro_ratio = 1
    # econ_ratio = 1
    # gamma_scenario = 1
    # farm_id = 36872 
    # K_scenario = 'int_1'
    
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
        crop_ids_by_farm_subset = {key: crop_ids_by_farm[key] for key in farm_id}  # farm_id = row in id_nldas, dictionary with farm_id as key and list of crops as value
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
        
        if gw_curve_option == 'internal':
		
            # load nldas grid cell aquifer properties for S, m, K, and WL
            nldas_id = farms_master.iloc[farm_id]['nldas']
            R = nldas_gw_attributes[(nldas_gw_attributes['NLDAS_ID'] == nldas_id.values[0])]['Recharge_usgs'] # annual average recharge (m/year)
            R = R * numer_hydro_factor 
            m_value = nldas_gw_attributes[(nldas_gw_attributes['NLDAS_ID'] == nldas_id.values[0])]['Depth'] # aquifer thickness (meters)
            WL_value = nldas_gw_attributes[(nldas_gw_attributes['NLDAS_ID'] == nldas_id.values[0])]['DTW'] # initial depth to water (meters below ground surface)
            grid_cell_area = nldas_gw_attributes[(nldas_gw_attributes['NLDAS_ID'] == nldas_id.values[0])]['Area'] # grid cell area (meters^2)
            
            # If assumed initial water level is be below aquifer depth from input datasets, set initial saturated thicknes to 50.9999 meters 
            if WL_value.values[0] >= m_value.values[0]:
                m_value.values[0] = WL_value.values[0] + 50.9999
            
            S_value = nldas_gw_attributes[(nldas_gw_attributes['NLDAS_ID'] == nldas_id.values[0])]['Porosity'] # porosity (unitless)
            IRR_DEPTH_afy = farm_gw_sum.iloc[farm_id[0]]['gw_irr_vol'] # baseline dataset annual groundwater volume for irrigation (acre-feet)            
            IRR_DEPTH = 1.25 * IRR_DEPTH_afy/farm_gw_area_sum.iloc[farm_id[0]]['area_irrigated_gw'] # IRR_DEPTH in feet with a 1.25 safety factor, updated IRR_DEPTH CALCULATED BASED ON BASALINE GW CROP AREA NOT TOTAL GRID CELL AREA
            
            if farm_gw_area_sum.iloc[farm_id[0]]['area_irrigated_gw'] == 0:
                IRR_DEPTH = 1
                
            # Prevent annual recharge rate from exceeding water level decline due to pumping 
            if R.values[0] >= IRR_DEPTH * 0.3048: # factor of 0.3048 m/ft used to convert IRR_DEPTH in feet to meters/year
                R.values[0] = IRR_DEPTH * 0.3048
                
            elif np.isnan(R.values[0]):
                R.values[0] = 0 # set recharge to 0 
                
            else:
                pass
                
            # Assign K_value from K_scenario name            
            if K_scenario == 'gleeson': # use average Gleeson value
                K_value = nldas_gw_attributes[(nldas_gw_attributes['NLDAS_ID'] == nldas_id.values[0])]['K']
                K_value = K_value.values[0]
            
            else: # use K_scenario value fromK K_dict dictionary 
                K_value = K_dict[K_scenario]
                            
            # Generate Cost Curve with Superwell_for_ABM_on_the_fly module 
            try:
                gw_cost_curve_internal = Superwell_for_ABM_on_the_fly.cost_curve(S_value.values[0], m_value.values[0], K_value, WL_value.values[0],
                                                                                 R.values[0], IRR_DEPTH * 0.3048, no_of_years, ELECTRICITY_RATE, grid_cell_area) 
                cost_curve = True
                
            except ValueError: # this occurs when a cost curve fails to produce due to a math domain error which can occur when T and S are very small
                pass 
                      
        # if the Superwell_for_ABM_on_the_fly fails to generate a cost curve (returns None), return all nan values
        if not gw_cost_curve_internal:
            
            #return [np.nan] * 89 # combined number of outputs for Total Simulation and Annual Outputs 
            cost_curve = False 
            gw_cost_curve_internal = [np.zeros(100), (np.zeros(100)), np.zeros(100)]
            max_gw_capacity = 0
        
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
            
        # alphas_total_subset = alphas_total_subset_og.copy() # uncomment to add Alpha to parameter sensitivity
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
                                           alphas_total_subset[key]) 

        # set land net prices to large negative value for gammas that are zero
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
        gw_cost_added = 0
        gw_availability_multiplier = 1  # SF added for reduced GW availability cost curve option
        
        if cost_curve == False:
            gw_availability_multiplier = .0000001
            gw_cost_added = 9999
            gw_calib_constr_subset[farm_id[0]]  = 0.0000001 
            for key in keys:
                net_prices_gw_sensitivity[key] = net_prices_gw_subset[key] * 9999
                
        time_depletion = 0  # SF added to reset time to depletion
        
        # MODIFICATION: CHECK IF SW COST > GW COST & IF SW AREA > 10%. IF 
        # TRUE, THE GW CONSTRAINT SET AT THE BASELINE. IF FALSE, THE GW CONSTRAINT 
        # WILL BE ALLOWED TO INCREASE AT A DEFINED RATE 
        if net_prices_sw_subset[farm_id[0]] < net_prices_gw_subset[farm_id[0]] and \
            sw_calib_constr_subset[farm_id[0]]/(sw_calib_constr_subset[farm_id[0]] + \
            gw_calib_constr_subset[farm_id[0]]) > 0.10: 
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
                keys_list = list(keys)
                gw_cost_initial = (-1.0 * net_prices_gw_sensitivity[keys_list[0]])/(nirs_subset[keys_list[0]] * 1000.0)
            
            else:
                net_prices_gw_subset_update = net_prices_gw_sensitivity.copy()

                for key in net_prices_gw_subset_update:
                    net_prices_gw_temp = net_prices_gw_subset_update[key]
                    nir_temp = nirs_subset[key]
                    gw_cost_initial = (-1.0 * net_prices_gw_temp) / (nir_temp * 1000.0) # initial GW cost $/acre-foot (same value for all 10 crop types)
                    gw_cost_updated_temp = gw_cost_initial + (gw_cost_added * 1233.48) # factor of 1233.48 converts $/m3 from cost curves to $/acre-ft used in PMP
                    net_prices_gw_updated_temp = -1.0 * gw_cost_updated_temp * nir_temp * 1000.0 # converts increase GW cost ($/acre-ft) to $/acre using crop-specific nir 
                    net_prices_gw_subset_update[key] = net_prices_gw_updated_temp # updates net GW price for each crop 
                
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
            initial_gw_capacity = gw_calib_constr_subset_copy[farm_id[0]]
            
            # MODIFICATION: THIS IS WHERE THE INCREASED ANNUAL GW VOLUME CONSTRAINT IS UPDATED. 
            # THIS INVOLVES MULTIPLYING THE GW CONSTRAINT BY 10% PERCENTAGE OF BASELINE OR ADDING 500 ACFT/YEAR, 
            # WHICHEVER IS GREATER. THE MAXIMUM EXPANISION OF GW CAPACITY IS CAPPED BY INITIAL NIR*CELL LAND CONSTRAINT AREA
            
            if gw_expansion == False:
                gw_calib_constr_subset_copy[farm_id[0]] = initial_gw_capacity * gw_availability_multiplier
                fwm_s.gw_constraints = Param(fwm_s.farm_ids, initialize=gw_calib_constr_subset_copy, mutable=True)
                
            else: # If expansion = True: Annual GW expansion the larger of 1) 10% of baseline OR 2) 500 acre-feet
                if t == 0:
                    gw_constr = gw_calib_constr_subset_copy
                    gw_capacity_increase = max(0.1 * gw_constr[farm_id[0]], 500)
                    
                    if gw_constr[farm_id[0]] + gw_capacity_increase > IRR_DEPTH * max_land_constr_subset[farm_id[0]]/1000:
                        gw_constr[farm_id[0]] =  IRR_DEPTH * max_land_constr_subset[farm_id[0]]/1000
                        max_gw_capacity = IRR_DEPTH * max_land_constr_subset[farm_id[0]]/1000

                    else:
                        gw_constr[farm_id[0]] = gw_constr[farm_id[0]] + gw_capacity_increase
                        max_gw_capacity = gw_constr[farm_id[0]] # assigned to max gw constraint in case the if statment is false before max_gw_capacity is assigned 
                      
                elif gw_constr[farm_id[0]] < IRR_DEPTH * max_land_constr_subset[farm_id[0]]/1000 and gw_availability_multiplier == 1: # limit area to land constraint, use baseline IRR_DEPTH defined by nir 
                    if gw_constr[farm_id[0]] + gw_capacity_increase > IRR_DEPTH * max_land_constr_subset[farm_id[0]]/1000:
                        gw_constr[farm_id[0]] =  IRR_DEPTH * max_land_constr_subset[farm_id[0]]/1000
                        max_gw_capacity =  IRR_DEPTH * max_land_constr_subset[farm_id[0]]/1000

                    else:
                        gw_constr[farm_id[0]] = gw_constr[farm_id[0]] + gw_capacity_increase
                        max_gw_capacity = gw_constr[farm_id[0]] # assigned to max gw constraint in case the elif statment is false before max_gw_capacity is assigned 
                        
                else: 
                    gw_constr[farm_id[0]] = max_gw_capacity * gw_availability_multiplier 
                    
                
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
                fwm_s.scaling_factor[fwm_s.c7] = 0.01   # choose scaling factor for profit constraint, tested 0.01, 1, 10 and no difference

            
            # create and run the optimization solver
            opt = SolverFactory("ipopt", solver_io='nl') # USE this for HPC 
            
            #opt = SolverFactory('gurobi', solver_io='python') # Stephen used for testing 
            # opt = SolverFactory('appsi_highs', solver_io='python')
            # opt = SolverFactory('asl:highs')#, solver_io='python')
            #opt.options["presolve"] = "on"
            #opt.options["parallel"] = "on"
            #opt.options["solver"] = "simplex"
            #opt.options["simplex_strategy"] = 1
            #pdb.set_trace()
            
            try:
                tol = 1e-6
                results = opt.solve(fwm_s, keepfiles=False, tee=False, options={
                    'constr_viol_tol': tol,
                    'max_iter': 3000,
                })

            except:
                try:
                    tol = 1e-11
                    results = opt.solve(fwm_s, keepfiles=False, tee=False, options={
                        'constr_viol_tol': tol,
                        'max_iter': 3000,
                    })

                except: 
                    bad_res = list(np.full(89, np.nan))
                    
                    for i in np.arange(52,89):
                        bad_res[i] = [np.nan]
                    bad_res[52] = [farm_id]
                    bad_res[53] = [hydro_ratio]
                    bad_res[54] = [econ_ratio]
                    bad_res[55] = [K_scenario]
                    bad_res[86] = [gamma_scenario]
                    bad_res[88] = tol
                    
                    return bad_res 

            # store main model outputs
            result_xs_sw = dict(fwm_s.xs_sw.get_values())
            result_xs_gw = dict(fwm_s.xs_gw.get_values())
            result_xs_total = dict(fwm_s.xs_total.get_values())

            farms_list = [farms_master.loc[[farm_id[0]]].nldas.values[0]]

            # process groundwater production results 
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
                results_pd['gw_mult'] = gw_availability_multiplier
                results_pd['gw_cost_added'] = gw_cost_added
            except NameError:
                results_pd['gw_cost_added'] = 0

            try:
                results_pd['gw_available'] = gw_availability_multiplier
            except NameError:
                results_pd['gw_available'] = 1

            try:
                results_pd['gw_cumul_vol'] = cumul_gw_sum
            except NameError:
                results_pd['gw_cumul_vol'] = 0


            # Create and append annual results
            if first:
                results_combined = results_pd
            else:
                results_combined = pd.concat([results_combined, results_pd])

            aggregation_functions = {'gw_vol': 'sum', 'xs_gw': 'sum'}
            gw_sum = results_pd[['nldas', 'gw_vol', 'xs_gw']].groupby(['nldas']).aggregate(aggregation_functions)
            gw_sum = gw_sum.reset_index()
            # gw_sum = gw_sum.merge(farm_gw_sum[['nldas', 'no_of_wells']], how='left', on=['nldas'])
            # gw_sum['gw_vol_well'] = gw_sum['gw_vol'] * farms_per_grid / gw_sum['no_of_wells']
            gw_sum['gw_vol_km3'] = gw_sum['gw_vol'] * 1.23348e-6 # unit conversion for gw_sum from acre-feet to km^3 using 1.23348e-6 multiplier
            
            if first:
                GW_pumped_timeseries = [] # empty list to track annual GW pumping
            
            GW_pumped_timeseries.append(gw_sum['gw_vol_km3'].values[0]) 

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
                GW_cost_increase = [] # Stephen - added empty array to track total increase in annual GW unit cost 
                
            if gw_curve_option == 'internal':
                for index in range(gw_cost_curve_internal[0].size):
                    
                    if cumul_gw_sum > max(gw_cost_curve_internal[1]):
                        #print("gw_vol_exceeds_cost_curve") 
                        i = 99
                        break 
                    
                    if cumul_gw_sum < gw_cost_curve_internal[1][index]:  # switched from cumul_gw_sum greater than to less than cost curve pumped volume at index closest match for cumulative pumped vol cost curve bin 
                        
                        if index == 0:
                            i = 0 
                            break
                        
                        else:
                            i = index - 1
                            break
  
                if gw_cost_curve_internal[0][i] != 0 and i != no_of_years-1: # if there is a unit cost at index i

                    if i == 0: # if the cumulative pumped volume is smaller than the first cumulative pumped volume bin
                        # absolute change in unit cost: current unit cost minus starting unit cost
                        gw_cost_added = 0
                        gw_cost = gw_cost_curve_internal[0][0] # groundwater cost is set to initial unit cost 


                        # change in annual pumping capacity due to reduce well pumping rates 
                        # if option == 'reduced_capacity':  # logic for optional reduced capacity 
                        gw_availability_multiplier = 1 # Stephen - do we want reduced capacity to be optional, right now it is the default 

                        WL_value = gw_cost_curve_internal[2][0] # depth to water is set as initial depth to water

                    else:
                        # interpolate GW unit cost
                        slope_incremental_cost = (gw_cost_curve_internal[0][i+1] - gw_cost_curve_internal[0][i])/(gw_cost_curve_internal[1][i+1] - gw_cost_curve_internal[1][i])  # new implementation with interpolation
                        dx = cumul_gw_sum - gw_cost_curve_internal[1][i] # excess volume of GW pumped above the previous cost curve pumped volume index
                        gw_cost_increase = dx * slope_incremental_cost # incremental cost increase compared to previous cost curve pumped volume index 
                        gw_cost = gw_cost_curve_internal[0][i] + gw_cost_increase # GW cost for current time step
                        gw_cost_added = gw_cost - gw_cost_curve_internal[0][0] # additional GW cost compared to initial cost curve GW cost for current annual time step

                        # interpolate WL
                        slope_incremental_WL = (gw_cost_curve_internal[2][i+1] - gw_cost_curve_internal[2][i])/(gw_cost_curve_internal[1][i+1] - gw_cost_curve_internal[1][i]) # new
                        wl_decline = dx * slope_incremental_WL # incremental increase in WL depth compared to previous cost curve pumped volume index 
                        WL_value =  gw_cost_curve_internal[2][i] + wl_decline # WL for current time step

                        # change in annual pumping capacity due to reduce well pumping rates
                        # if option == 'reduced_capacity':  # logic for optional reduced capacity
                        gw_availability_multiplier = gw_cost_curve_internal[3][i] # reduced capacity to be optional, right now it is the default

                    # update WL 
                    WL_timeseries.append(WL_value) # new

                    # update GW cost  
                    GW_cost_timeseries.append(gw_cost_initial + gw_cost_added * 1233.48) # initial GW cost from Farm database ($/acre-ffot) + gw_cost_added (converted to $/acre foot)
                    
                    # update GW added cost  
                    GW_cost_increase.append(gw_cost_added)

                else:
                    gw_cost_added = 9999
                    gw_availability_multiplier = .00001
                    
                    # update WL 
                    WL_timeseries.append(max(gw_cost_curve_internal[2]))

                    # update GW cost - Stephen added - what cost to assign to output timeseries when no GW is available? 999 for now 
                    GW_cost_timeseries.append(9999)
                    
                    # update GW added cost - Stephen added 
                    GW_cost_increase.append(9999)
                    
                    if depletion_first and i < no_of_years-1:
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
        m_value = m_value.values[0]

        
        # Stephen - New Total Simulation metrics 
        GW_cost_filtered = GW_cost_timeseries.copy()
        GW_cost_filtered = np.array(GW_cost_filtered)
        GW_cost_indexes = np.where(GW_cost_filtered < 2000)
        GW_cost_filtered = GW_cost_filtered[GW_cost_indexes]
        
        if cost_curve == True:   
            Total_increase_GW_cost = max(GW_cost_filtered) - min(GW_cost_filtered)
            Initial_WL = WL_timeseries[0] 
            Final_WL = max(WL_timeseries)
        else:
            Total_increase_GW_cost = 0
            Initial_WL = WL_value.values[0]
            Final_WL = WL_value.values[0]
        
            
        Area_GW_summary = sum(results_combined.xs_gw)
        Area_SW_summary = sum(results_combined.xs_sw)
        
        Corn = results_combined.where(results_combined.crop == 'Corn').dropna()
        Area_GW_Corn_summary = sum(Corn.xs_gw)
        Area_SW_Corn_summary = sum(Corn.xs_sw)
        
        FiberCrop = results_combined.where(results_combined.crop == 'FiberCrop').dropna()
        Area_GW_FiberCrop_summary = sum(FiberCrop.xs_gw)
        Area_SW_FiberCrop_summary = sum(FiberCrop.xs_sw)
        
        FodderGrass = results_combined.where(results_combined.crop == 'FodderGrass').dropna()
        Area_GW_FodderGrass_summary = sum(FodderGrass.xs_gw)
        Area_SW_FodderGrass_summary = sum(FodderGrass.xs_sw)
        
        MiscCrop = results_combined.where(results_combined.crop == 'MiscCrop').dropna()
        Area_GW_MiscCrop_summary = sum(MiscCrop.xs_gw)
        Area_SW_MiscCrop_summary = sum(MiscCrop.xs_sw)
        
        OilCrop = results_combined.where(results_combined.crop == 'OilCrop').dropna()
        Area_GW_OilCrop_summary = sum(OilCrop.xs_gw)
        Area_SW_OilCrop_summary = sum(OilCrop.xs_sw)
        
        OtherGrain = results_combined.where(results_combined.crop == 'OtherGrain').dropna()
        Area_GW_OtherGrain_summary = sum(OtherGrain.xs_gw)
        Area_SW_OtherGrain_summary = sum(OtherGrain.xs_sw)
        
        Rice = results_combined.where(results_combined.crop == 'Rice').dropna()
        Area_GW_Rice_summary = sum(Rice.xs_gw)
        Area_SW_Rice_summary = sum(Rice.xs_sw)
        
        Root_Tuber = results_combined.where(results_combined.crop == 'Root_Tuber').dropna()
        Area_GW_Root_Tuber_summary = sum(Root_Tuber.xs_gw)
        Area_SW_Root_Tuber_summary = sum(Root_Tuber.xs_sw)
        
        SugarCrop = results_combined.where(results_combined.crop == 'SugarCrop').dropna()
        Area_GW_SugarCrop_summary = sum(SugarCrop.xs_gw)
        Area_SW_SugarCrop_summary = sum(SugarCrop.xs_sw)
        
        Wheat = results_combined.where(results_combined.crop == 'Wheat').dropna()
        Area_GW_Wheat_summary = sum(Wheat.xs_gw)
        Area_SW_Wheat_summary = sum(Wheat.xs_sw)
        
        
        if cost_curve == True:
            perc_vol_depleted = (max(WL_timeseries) - gw_cost_curve_internal[2][0])/(m_value - gw_cost_curve_internal[2][0])
            
            if perc_vol_depleted > 1:
                perc_vol_depleted = 1
                
            
            # Time (years) to depletion thresholds [10m, 20m, 30m, 40m, 50m]
            depletion_ts = WL_timeseries - WL_timeseries[0]
            
            t_depletion_10m = t_depletion(depletion_ts, no_of_years, 10)
            t_depletion_20m = t_depletion(depletion_ts, no_of_years, 20)
            t_depletion_30m = t_depletion(depletion_ts, no_of_years, 30)
            t_depletion_40m = t_depletion(depletion_ts, no_of_years, 40)
            t_depletion_50m = t_depletion(depletion_ts, no_of_years, 50)
            
        else:
            perc_vol_depleted = 0 
            
            t_depletion_10m = -999
            t_depletion_20m = -999
            t_depletion_30m = -999
            t_depletion_40m = -999
            t_depletion_50m = -999
              
        # Annual Metrics - 
        sim_years = max(results_combined.year) # years of output 
        
        Farm_id = np.tile(farm_id, sim_years).tolist()
        Hydro_ratio = np.tile(hydro_ratio, sim_years).tolist()
        Econ_ratio = np.tile(econ_ratio, sim_years).tolist()
        K_scenario_list = np.tile(K_scenario, sim_years).tolist()
        K_values = np.tile(K_value, sim_years).tolist()
        m_values_arr = np.tile(m_value, sim_years).tolist()
        Year = np.arange(1,sim_years+1).tolist()
        GW_WL = WL_timeseries # populate using WL_timeseries
        GW_vol_pumped = GW_pumped_timeseries # populate using GW_pumped_timeseries 
        
        if cost_curve == True:
            GW_fraction_depleted = ((m_values_arr-np.array(WL_timeseries))/(m_value-min(WL_timeseries))).tolist() # caculate using WL_timeseries, aquifer depth, and initial sat thickness 
        else:
            GW_fraction_depleted = np.zeros(int(no_of_years)).tolist()
        
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
                med_gw_vol_yr, mean_excl0_gw_vol_yr, Initial_WL, Final_WL, Total_increase_GW_cost, Area_GW_summary, 
                Area_SW_summary,  Area_GW_Corn_summary, Area_SW_Corn_summary, Area_GW_FiberCrop_summary, Area_SW_FiberCrop_summary, 
                Area_GW_FodderGrass_summary, Area_SW_FodderGrass_summary, Area_GW_MiscCrop_summary, Area_SW_MiscCrop_summary, 
                Area_GW_OilCrop_summary, Area_SW_OilCrop_summary, Area_GW_OtherGrain_summary, Area_SW_OtherGrain_summary, 
                Area_GW_Rice_summary, Area_SW_Rice_summary, Area_GW_Root_Tuber_summary, Area_SW_Root_Tuber_summary, 
                Area_GW_SugarCrop_summary, Area_SW_SugarCrop_summary, Area_GW_Wheat_summary, Area_SW_Wheat_summary, 
                t_depletion_10m, t_depletion_20m, t_depletion_30m, t_depletion_40m, t_depletion_50m,
                
                
                Farm_id, Hydro_ratio, Econ_ratio, K_scenario_list, K_values, m_value, Year, 
                GW_WL, GW_vol_pumped, GW_fraction_depleted, GW_cost, 
                Profit, Area_GW_Corn, Area_SW_Corn, Area_GW_FiberCrop, Area_SW_FiberCrop,
                Area_GW_FodderGrass, Area_SW_FodderGrass, Area_GW_MiscCrop, Area_SW_MiscCrop,
                Area_GW_OilCrop, Area_SW_OilCrop, Area_GW_OtherGrain, Area_SW_OtherGrain, 
                Area_GW_Rice, Area_SW_Rice, Area_GW_Root_Tuber, Area_SW_Root_Tuber, 
                Area_GW_SugarCrop, Area_SW_SugarCrop, Area_GW_Wheat, Area_SW_Wheat,
                Area_GW_total, Area_SW_total, Gamma_mult, S_values, tol]
                

    ##### Travis: The following code defines the ensemble information, which is consolidated in the "cases_df" pandas dataframe    
    # Varied parameters for ABM sensitivity
        
    hydro_ratio = [0.5, 0.75, 1, 1.25, 1.5] # 
    econ_ratio = [0.5, 0.75, 1, 1.25, 1.5] # 
    K_dict = dict(low = 0.5, int_1 = 2.5, int_2 = 10.0, high = 50.0, gleeson = 0.0) # dictionary shows values in m/d
    K_scenario = ['low', 'int_1', 'int_2', 'high', 'gleeson']
    gamma_scenario = [0.5, 0.75, 1, 1.25, 1.5] # gamma multiplier scenario values 
    
    ##### Travis - Here is where I provide a list of farms by ID that we subsequently loop through. For the comprehensive farm sweep on HPC, the idea would be to run for every farm in our database (0-53834) and strategize a way to parallelize this on HPC (each call of farm_gw_model is independent)
    farms = [fid]
    
    # Create array of ABM sensitivity parameter combinations using parameter lists defined in 'Varied parameters for ABM sensitivity'
    combinations = len(hydro_ratio) * len(econ_ratio) * len(K_scenario) * len(gamma_scenario) * len(farms)
    cases = list(product(hydro_ratio, econ_ratio, K_scenario, gamma_scenario, farms))
    cases_df = pd.DataFrame(data = cases.copy(), columns = ['hydro_ratio', 'econ_ratio', 'K_scenario', 'gamma_scenario', 'farm']) 
                    
    # Create empty DataFrame for Annual outputs, should not pre-allocate 
    # Annual_df = pd.DataFrame()

    ### Create Output arrays to store resulting data
    
    # Total 100 year summary metrics 
    results_end_div_start_area = []
    results_total_profit = []
    results_perc_vol_depleted = []
    results_time_depletion = []
    results_cum_gw = []
    acres_grown_total = []
    acres_grown_mean = []
    acres_grown_mean_vs_start =  []
    end_div_start_profit =  []
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
    Initial_WL = []
    Final_WL = []
    Total_increase_GW_cost = []
    Area_GW_summary = []
    Area_SW_summary = []
    Area_GW_Corn_summary = []
    Area_SW_Corn_summary = []
    Area_GW_FiberCrop_summary = []
    Area_SW_FiberCrop_summary = []
    Area_GW_FodderGrass_summary = []
    Area_SW_FodderGrass_summary = []
    Area_GW_MiscCrop_summary = []
    Area_SW_MiscCrop_summary = []
    Area_GW_OilCrop_summary = []
    Area_SW_OilCrop_summary = []
    Area_GW_OtherGrain_summary = []
    Area_SW_OtherGrain_summary = []
    Area_GW_Rice_summary = []
    Area_SW_Rice_summary = []
    Area_GW_Root_Tuber_summary = []
    Area_SW_Root_Tuber_summary  = []
    Area_GW_SugarCrop_summary = []
    Area_SW_SugarCrop_summary = []
    Area_GW_Wheat_summary = []
    Area_SW_Wheat_summary = []
    t_depletion_10m = []
    t_depletion_20m = []
    t_depletion_30m = []
    t_depletion_40m = []
    t_depletion_50m = [] 
    m_value = []
    tol = []

    ## Annual metrics 
    # Farm_id = []
    # Hydro_ratio = [] 
    # Econ_ratio = []
    # K_scenario = []
    # K_value = []
    # m_value_arr = [] 
    # Year = []
    # GW_WL = []
    # GW_vol_pumped = []
    # GW_fraction_depleted = []
    # GW_cost = []
    # Profit = []
    # Area_GW_Corn = []
    # Area_SW_Corn = []
    # Area_GW_FiberCrop = [] 
    # Area_SW_FiberCrop = []
    # Area_GW_FodderGrass = [] 
    # Area_SW_FodderGrass = []
    # Area_GW_MiscCrop = []
    # Area_SW_MiscCrop = []
    # Area_GW_OilCrop = []
    # Area_SW_OilCrop = []
    # Area_GW_OtherGrain = [] 
    # Area_SW_OtherGrain = []
    # Area_GW_Rice = []
    # Area_SW_Rice = []
    # Area_GW_Root_Tuber = []
    # Area_SW_Root_Tuber = []
    # Area_GW_SugarCrop = []
    # Area_SW_SugarCrop = []
    # Area_GW_Wheat = []
    # Area_SW_Wheat = []
    # Area_GW_total = []
    # Area_SW_total = []
    # Gamma_mult  = []
    # S_value = []
    
    ### Travis - this is where we start the simulations, looping through every case in the cases_df dataframe
    ### Simulate parameter combinations over all selected farms
    for case in range(combinations):

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
        Initial_WL.append(results[22])
        Final_WL.append(results[23])
        Total_increase_GW_cost.append(results[24])
        Area_GW_summary.append(results[25])
        Area_SW_summary.append(results[26])
        Area_GW_Corn_summary.append(results[27])
        Area_SW_Corn_summary.append(results[28])
        Area_GW_FiberCrop_summary.append(results[29])
        Area_SW_FiberCrop_summary.append(results[30])
        Area_GW_FodderGrass_summary.append(results[31])
        Area_SW_FodderGrass_summary.append(results[32])
        Area_GW_MiscCrop_summary.append(results[33])
        Area_SW_MiscCrop_summary.append(results[34])
        Area_GW_OilCrop_summary.append(results[35])
        Area_SW_OilCrop_summary.append(results[36])
        Area_GW_OtherGrain_summary.append(results[37])
        Area_SW_OtherGrain_summary.append(results[38])
        Area_GW_Rice_summary.append(results[39])
        Area_SW_Rice_summary.append(results[40])
        Area_GW_Root_Tuber_summary.append(results[41])
        Area_SW_Root_Tuber_summary.append(results[42])
        Area_GW_SugarCrop_summary.append(results[43])
        Area_SW_SugarCrop_summary.append(results[44])
        Area_GW_Wheat_summary.append(results[45])
        Area_SW_Wheat_summary.append(results[46])
        t_depletion_10m.append(results[47])
        t_depletion_20m.append(results[48])
        t_depletion_30m.append(results[49])
        t_depletion_40m.append(results[50])
        t_depletion_50m.append(results[51])
        m_value.append(results[57])
        tol.append(results[88])
        
        # Stephen - Output arrays for Annual values 
        # try:
        #     Farm_id.extend(results[52])
        #     Hydro_ratio.extend(results[53]) 
        #     Econ_ratio.extend(results[54])
        #     K_scenario.extend(results[55])
        #     K_value.extend(results[56])
        #     m_value_arr.extend(results[57])
        #     Year.extend(results[58])
        #     GW_WL.extend(results[59])
        #     GW_vol_pumped.extend(results[60])
        #     GW_fraction_depleted.extend(results[61])
        #     GW_cost.extend(results[62])
        #     Profit.extend(results[63])
        #     Area_GW_Corn.extend(results[64])
        #     Area_SW_Corn.extend(results[65])
        #     Area_GW_FiberCrop.extend(results[66]) 
        #     Area_SW_FiberCrop.extend(results[67])
        #     Area_GW_FodderGrass.extend(results[68]) 
        #     Area_SW_FodderGrass.extend(results[69])
        #     Area_GW_MiscCrop.extend(results[70])
        #     Area_SW_MiscCrop.extend(results[71])
        #     Area_GW_OilCrop.extend(results[72])
        #     Area_SW_OilCrop.extend(results[73])
        #     Area_GW_OtherGrain.extend(results[74]) 
        #     Area_SW_OtherGrain.extend(results[75])
        #     Area_GW_Rice.extend(results[76])
        #     Area_SW_Rice.extend(results[77])
        #     Area_GW_Root_Tuber.extend(results[78])
        #     Area_SW_Root_Tuber.extend(results[79])
        #     Area_GW_SugarCrop.extend(results[80])
        #     Area_SW_SugarCrop.extend(results[81])
        #     Area_GW_Wheat.extend(results[82])
        #     Area_SW_Wheat.extend(results[83])
        #     Area_GW_total.extend(results[84])
        #     Area_SW_total.extend(results[85])
        #     Gamma_mult.extend(results[86])
        #     S_value.extend(results[87])
            
        # except TypeError:
        #     continue 
        

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
    cases_df['area_GW'] = Area_GW_summary 
    cases_df['area_SW'] = Area_SW_summary
    cases_df['area_GW_corn'] = Area_GW_Corn_summary 
    cases_df['area_SW_corn'] = Area_SW_Corn_summary 
    cases_df['area_GW_FiberCrop'] = Area_GW_FiberCrop_summary 
    cases_df['area_SW_FiberCrop'] = Area_SW_FiberCrop_summary 
    cases_df['area_GW_FodderGrass'] = Area_GW_FodderGrass_summary 
    cases_df['area_SW_FodderGrass'] = Area_SW_FodderGrass_summary 
    cases_df['area_GW_MiscCrop'] = Area_GW_MiscCrop_summary
    cases_df['area_SW_MiscCrop'] = Area_SW_MiscCrop_summary
    cases_df['area_GW_OilCrop'] = Area_GW_OilCrop_summary 
    cases_df['area_SW_OilCrop'] = Area_SW_OilCrop_summary 
    cases_df['area_GW_OtherGrain'] = Area_GW_OtherGrain_summary
    cases_df['area_SW_OtherGrain'] = Area_SW_OtherGrain_summary
    cases_df['area_GW_Rice'] = Area_GW_Rice_summary
    cases_df['area_SW_Rice'] = Area_SW_Rice_summary
    cases_df['area_GW_Root_Tuber'] = Area_GW_Root_Tuber_summary 
    cases_df['area_SW_Root_Tuber'] = Area_SW_Root_Tuber_summary 
    cases_df['area_GW_SugarCrop'] = Area_GW_SugarCrop_summary 
    cases_df['area_SW_SugarCrop'] = Area_SW_SugarCrop_summary 
    cases_df['area_GW_Wheat'] = Area_GW_Wheat_summary
    cases_df['area_SW_Wheat'] = Area_SW_Wheat_summary
    cases_df['t_depletion_10m'] = t_depletion_10m 
    cases_df['t_depletion_20m'] = t_depletion_20m
    cases_df['t_depletion_30m'] = t_depletion_30m 
    cases_df['t_depletion_40m'] = t_depletion_40m 
    cases_df['t_depletion_50m'] = t_depletion_50m  
    cases_df['m'] = m_value 
    cases_df['constraint_tolerance'] = tol
        
        
    ## Save results from the annual outputs to a DataFrame
    # Annual_df['Farm_id'] = Farm_id
    # Annual_df['Hydro_ratio'] = Hydro_ratio
    # Annual_df['Econ_ratio'] = Econ_ratio
    # Annual_df['K_scenario'] = K_scenario
    # Annual_df['K_value'] = K_value
    # Annual_df['m_value'] = m_value_arr
    # Annual_df['Year'] = Year
    # Annual_df['GW_WL'] = GW_WL
    # Annual_df['GW_vol_pumped'] = GW_vol_pumped
    # Annual_df['GW_fraction_depleted'] = GW_fraction_depleted
    # Annual_df['GW_cost'] = GW_cost
    # Annual_df['Profit'] = Profit
    # Annual_df['Area_GW_Corn'] = Area_GW_Corn
    # Annual_df['Area_SW_Corn'] = Area_SW_Corn
    # Annual_df['Area_GW_FiberCrop'] = Area_GW_FiberCrop
    # Annual_df['Area_SW_FiberCrop'] = Area_SW_FiberCrop
    # Annual_df['Area_GW_FodderGrass'] = Area_GW_FodderGrass
    # Annual_df['Area_SW_FodderGrass'] = Area_SW_FodderGrass
    # Annual_df['Area_GW_MiscCrop'] = Area_GW_MiscCrop
    # Annual_df['Area_SW_MiscCrop'] = Area_SW_MiscCrop
    # Annual_df['Area_GW_OilCrop'] = Area_GW_OilCrop
    # Annual_df['Area_SW_OilCrop'] = Area_SW_OilCrop
    # Annual_df['Area_GW_OtherGrain'] = Area_GW_OtherGrain
    # Annual_df['Area_SW_OtherGrain'] = Area_SW_OtherGrain
    # Annual_df['Area_GW_Rice'] =  Area_GW_Rice
    # Annual_df['Area_SW_Rice'] =  Area_SW_Rice
    # Annual_df['Area_GW_Root_Tuber'] = Area_GW_Root_Tuber
    # Annual_df['Area_SW_Root_Tuber'] =  Area_SW_Root_Tuber
    # Annual_df['Area_GW_SugarCrop'] =  Area_GW_SugarCrop
    # Annual_df['Area_SW_SugarCrop'] =  Area_SW_SugarCrop
    # Annual_df['Area_GW_Wheat'] =  Area_GW_Wheat
    # Annual_df['Area_SW_Wheat'] =  Area_SW_Wheat
    # Annual_df['Area_GW_total'] =  Area_GW_total
    # Annual_df['Area_SW_total'] =  Area_SW_total
    # Annual_df['Gamma_mult'] =  Gamma_mult
    # Annual_df['S_value'] =  S_value
        
        
    ##### Travis: And export to csv. For each farm (~50k), we will have 625 scenarios
    cases_df.to_csv(f'{output_path}/farm_{fid}_cases.csv', index=False)
    #Annual_df.to_csv(f'{output_path}/farm_{fid}_annual.csv', index=False) # Stephen added 

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
