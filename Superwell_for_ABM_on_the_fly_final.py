# -*- coding: utf-8 -*-

import numpy as np 
import pandas as pd 
import math 
import os 

# load well function lookup table 
W_lookup = pd.read_csv('Theis_well_function_table.csv', header = "infer") 
lookup_idx = pd.Index(W_lookup.W)

#%% define Theis function 

def drawdown_theis(time, r, S, T, Q):
    u = r**2 * S/(4 * T * time)
    
    if u > 5.9: # for large u values, W will be insignificant and drawdown (s) will ~= 0 
        W = 0
        
    elif 5.9 > u and u > .6: # use W(u) lookup table for intermediate values where approximation is insufficient 
        lookup_idx = pd.Index(W_lookup.u)
        lookup_loc = lookup_idx.get_loc(u, method = 'nearest')
        W = W_lookup.W[lookup_loc]
        
    else: # use approximation for small u values (< 0.6)
        W = -0.57721 - math.log(u) + u - u**2/(2*2)
        
    s = W * Q / (4 * 3.1416 * T)
    
    return(s, u)


#%% Cost curve function
    
def cost_curve(S, m, K, WL, R, IRR_DEPTH, NUM_YEARS, ELECTRICITY_RATE, grid_cell_area, INCREASE_WELL_DENSITY = False,  EXPONENTIAL_K_DECLINE = False):
    
    # constants used by function 
    DAYS = 100 # days of pumping per year, code designed for increments of 10 
    #DEPLETION_LIMIT = 0.95 
    EFFICIENCY = 0.7 # well efficiency
    WELL_LIFETIME = 30 # years before well needs to be replaced 
    INTEREST_RATE = 0.1 # loan rate (0.1 = 10%)
    MAINTENANCE_RATE = 0.07 # annual cost in % of initial cost (0.1 = 10%)
    SPECIFIC_WEIGHT = 9800 # specific weight of water Newtons/m^3
    well_unit_cost = 350 # $/m
    MAXIMUM_INITIAL_SAT_THICKNESS = 100 # m 
    
    # user defined constants used by function 
    IRR_DEPTH = IRR_DEPTH
    NUM_YEARS = NUM_YEARS
    ELECTRICITY_RATE = ELECTRICITY_RATE
    RECHARGE = R # annual average recharge (m) 

    total_thickness = m # m 
    
    #grid_cell_area = 13375 * 13375 # m^2 in 1/8th deg grid cell 
    
    # depth to water 
    depth_to_peiz_surface = WL # m
    DTW_array = np.zeros(NUM_YEARS) # tracks depth to water for each year 
    DTW_array[0] = depth_to_peiz_surface # initial depth to water
    
    # saturated thickness 
    initial_sat_thickness = m - WL # m
    
    if initial_sat_thickness > MAXIMUM_INITIAL_SAT_THICKNESS:
        sat_thickness_array = np.zeros(NUM_YEARS)
        sat_thickness_array[0] = MAXIMUM_INITIAL_SAT_THICKNESS # m
        well_length_array = np.zeros(NUM_YEARS)
        well_length_array[0] = sat_thickness_array[0] + DTW_array[0] # m
    else:    
        sat_thickness_array = np.zeros(NUM_YEARS)
        sat_thickness_array[0] = initial_sat_thickness # m
        well_length_array = np.zeros(NUM_YEARS)
        well_length_array[0] = total_thickness # m
        
    # available volume (meters^3)
    available_volume = initial_sat_thickness * grid_cell_area * S
    
    # aquifer properties for Theis 
    S = S # storativity = porosity [-]
    
    if EXPONENTIAL_K_DECLINE == False: 
        
        K = K/86400 # convert K in m/d to m/s 
        T = K * sat_thickness_array[0] # m^2/s
        T_array = np.zeros(NUM_YEARS) # tracks T for each year 
        T_array[0] = T # initial T
    
    # CODE BLOCK FOR TRANSMISSIVITY DEFINED BY EXPONENTIAL K DECLINE WITH DEPTH
    else:
        alpha = 100
        depth_profile = np.linspace(1, int(np.ceil(total_thickness)), int(np.ceil(total_thickness)))
        K_profile = np.zeros(len(depth_profile))
        
        for depth_interval in range(len(depth_profile)):
            K_profile[depth_interval] = K * math.e**(-1*(depth_interval+1)/alpha)
    
        K_avg = np.mean(K_profile[DTW_array[0]:well_length_array[0]]) # m/s
        T = K_avg * sat_thickness_array[0] # m^2/s
        T_array = np.zeros(NUM_YEARS) # tracks T for each year 
        T_array[0] = T # initial T
    
    #################### determine initial well Q #############################
    
    # time and well radius for Theis solution
    time_Q = 2 * 365 * 86400 # time period used (seconds) for determining initial well Q
    well_r = 0.5 * 0.28 # 1/2 * well diameter in m  
    
    # candidate well pumping rates (gallons per minute)
    Q_array_gpm = np.linspace(50, 1000, 96)

    Q_array = np.zeros(len(Q_array_gpm))
    # convert candidate pumping rates to m^3/s
    for i, Q in enumerate(Q_array_gpm):
        Q_array[i] = Q/(60*264.17)
        
    # drawdown at t = 2 years for all candidate well Qs 
    s_array = np.zeros(len(Q_array))
    u_array = np.zeros(len(Q_array))
    for i, Q in enumerate(Q_array):
        s_array[i], u_array[i] = drawdown_theis(time_Q, well_r, S, T, Q)
    
    # find largest Q that meets screening criteria
    # screening criteria 
    max_s_frac = .35     # max drawdown as % of available sat thickness
    max_s_absolute = 50  # max drawdown in m
    
    Q_viability = np.zeros(len(Q_array))
    
    for i, s in enumerate(s_array):
        if s/initial_sat_thickness < max_s_frac and s < max_s_absolute and s != 0:
            Q_viability[i] = 1
    
    # skip grid cell if no pumping rates are viable
    if np.sum(Q_viability) == 0:
        return print("Cost curve not produced, no viable pumping rates")
    
    initial_Q_indx_arr = np.where(Q_viability == 1) 
    initial_Q_indx = np.max(initial_Q_indx_arr[:]) # index of largest viable Q
    initial_Q = Q_array[initial_Q_indx]
    Well_Q_array = np.zeros(NUM_YEARS)
    Well_Q_array[0] = initial_Q 
    
    ###################### determine initial well Area ############################
    initial_well_area = initial_Q * DAYS * 86400 / (IRR_DEPTH) # m^2
    initial_roi = (initial_well_area/math.pi) ** 0.5 # m 
    well_roi_array = np.zeros(NUM_YEARS)
    well_roi_array[0] = initial_roi
    well_area_array = np.zeros(NUM_YEARS)
    well_area_array[0] = initial_well_area
    
    ####################### annual pumping simulation loop #######################
    depleted_volume_fraction = np.zeros(NUM_YEARS) # initialize 
    
    for year in range(NUM_YEARS):
        # if depleted_volume_fraction[year-1] > DEPLETION_LIMIT:
        #     year = year - 1 
        #     break 
            
        # test viability for current year (simulate drawdown at t = 100 days of pumping)
        # initialize viability variables 
        s_theis = 0 
        u_theis = 0
        s_theis_interference = 0 
        u_theis_interference = 0
        
        s_theis, u_theis = drawdown_theis(DAYS * 86400, well_r, S, T_array[year], Well_Q_array[year])
        s_theis_interference, u_theis_interference = drawdown_theis(DAYS * 86400, well_roi_array[year] * 2, S, T_array[year], Well_Q_array[year])
        s_total = s_theis + 4 * s_theis_interference # total drawdown (well + interference)
        
        # check if drawdown constraints are violated by end of 100 day pumping period
        # if constraints violated: (1) first deepen well, (2) then reduce well pumping rate 
        
        max_s_absolute_constraint = 60
        max_s_frac_constraint = 0.4
        
        if s_total > max_s_absolute_constraint or s_total/sat_thickness_array[year] > max_s_frac_constraint:
            
            # 1) first preference deepen well 
            if well_length_array[year] < total_thickness:
                # update well length
                if well_length_array[year] + 50 < total_thickness:
                    well_length_array[year] = 50 + well_length_array[year]
                
                else:
                    remaining_length = total_thickness - well_length_array[year]
                    well_length_array[year] = remaining_length + well_length_array[year]
                    
                # update saturated thickness and T 
                if EXPONENTIAL_K_DECLINE == False: 
                    sat_thickness_array[year] = well_length_array[year] - DTW_array[year]
                    T_array[year] = sat_thickness_array[year] * K
                    
                else: 
                    K_avg = np.mean(K_profile[DTW_array[year]:well_length_array[year]]) # m/s
                    T_array[year] = K_avg * sat_thickness_array[year] # m^2/s
                    
                
            # 2) once well cannot be deepened, reduce well pumping rate 
            else:
                s_array = np.zeros(len(Q_array))
                u_array = np.zeros(len(Q_array))
                time_Q_updated = DAYS * 86400 # DAYS of pumping in 1 year of pumping 
                for i, Q in enumerate(Q_array):
                    s_array[i], u_array[i] = drawdown_theis(time_Q_updated, well_r, S, T_array[year], Q)
            
                Q_viability = np.zeros(len(Q_array))
                
                for i, s in enumerate(s_array):
                    if s/sat_thickness_array[year] < max_s_frac and s < max_s_absolute:
                        Q_viability[i] = 1
                
                # Exit pumping code block if no pumping rates are viable 
                if np.sum(Q_viability) == 0:
                    break 
                
                Q_indx_arr = np.where(Q_viability == 1) 
                Q_indx = np.max(Q_indx_arr[:]) # index of largest viable Q
                new_Q = Q_array[Q_indx] # new Q
                if new_Q > Well_Q_array[0]:
                    new_Q = Q_array[initial_Q_indx]
                    
                Well_Q_array[year] = new_Q # update Q for current YEAR 
                
                # calculate well area and new well roi 
                if INCREASE_WELL_DENSITY == False:
                    well_roi_array[year] = initial_roi
                    well_area_array[year] =  initial_well_area
                    
                else:
                    well_area_array[year] = Well_Q_array[year] * DAYS * 86400 / (IRR_DEPTH)
                    well_roi = (well_area_array[year] / math.pi) ** 0.5
                    well_roi_array[year] = well_roi
        
           
        # if constraints aren't violated, proceed to calculate output for pumping year  
        # simulate annual pumping, with drawdown calculated every 10 days
        s_theis_ts = np.zeros(int(DAYS/10)) 
        u_theis_ts = np.zeros(int(DAYS/10)) 
        s_theis_interference_ts = np.zeros(int(DAYS/10))
        u_theis_interference_ts = np.zeros(int(DAYS/10))
        
        for day in range(int(DAYS/10)):
            s_theis_ts[day], u_theis_ts[day] = drawdown_theis((day+1) * 10 * 86400, well_r, S, T_array[year], Well_Q_array[year])
            s_theis_interference_ts[day], u_theis_interference_ts[day] = drawdown_theis((day+1) * 10 * 86400, well_roi_array[year] * 2, S, T_array[year], Well_Q_array[year])
        
        # average drawdown
        s_theis_avg = np.mean(s_theis_ts) + np.mean(4* s_theis_interference_ts) 
        s_interference_avg = 4 * np.mean(s_theis_interference_ts) 
        
        # convert to Jacob - solve quadratic
        a = -1/(2*sat_thickness_array[year])
        b = 1
        c = -s_theis_avg
        
        root_1 = (-b + (b**2 - 4 * a * c) ** 0.5)/(2*a) 
        root_2 = (-b - (b**2 - 4 * a * c) ** 0.5)/(2*a) 
        
        s_jacob = root_1
        
        ########################### compute outputs ###########################
        
        # save annual pumping values to arrays 
        if year == 0:
            drawdown = np.zeros(NUM_YEARS)
            drawdown_interference = np.zeros(NUM_YEARS)
            total_head = np.zeros(NUM_YEARS)
            volume_per_well = np.zeros(NUM_YEARS)
            num_wells = np.zeros(NUM_YEARS)
            volume_all_wells = np.zeros(NUM_YEARS)
            cumulative_volume_per_well = np.zeros(NUM_YEARS)
            cumulative_volume_all_wells = np.zeros(NUM_YEARS)
            depleted_volume_fraction = np.zeros(NUM_YEARS)
            
            drawdown[year] = s_jacob
            drawdown_interference[year] = s_interference_avg
            total_head[year] = s_jacob + DTW_array[year]
            volume_per_well[year] = Well_Q_array[year] * 86400 * DAYS
            num_wells[year] = grid_cell_area/well_area_array[year]
            volume_all_wells[year] = volume_per_well[year] * num_wells[year]
            cumulative_volume_per_well[year] = volume_per_well[year]
            cumulative_volume_all_wells[year] =  volume_all_wells[year]
            depleted_volume_fraction[year] = cumulative_volume_all_wells[year]/available_volume
            
        else:
            drawdown[year] = s_jacob
            drawdown_interference[year] = s_interference_avg
            total_head[year] = s_jacob + DTW_array[year]
            volume_per_well[year] = Well_Q_array[year] * 86400 * DAYS
            num_wells[year] = grid_cell_area/well_area_array[year]
            volume_all_wells[year] = volume_per_well[year] * num_wells[year]
            cumulative_volume_per_well[year] = volume_per_well[year] + cumulative_volume_per_well[year-1]
            cumulative_volume_all_wells[year] =  volume_all_wells[year] +  cumulative_volume_all_wells[year-1]
            depleted_volume_fraction[year] = cumulative_volume_all_wells[year]/available_volume
    
        
        # update variable arrays for next annual pumping iteration
        if year != NUM_YEARS-1:
            Well_Q_array[year+1] = Well_Q_array[year]
            DTW_array[year+1] = DTW_array[year] + (volume_all_wells[year]/grid_cell_area)/S - RECHARGE/S
            sat_thickness_array[year+1] = well_length_array[year] - DTW_array[year+1]
            
            if EXPONENTIAL_K_DECLINE == False: 
                T_array[year+1] = K * sat_thickness_array[year+1]
                
            else:
                K_avg = np.mean(K_profile[DTW_array[year]:well_length_array[year]]) # m/s
                T_array[year] = K_avg * sat_thickness_array[year] # m^2/s
           
            well_roi_array[year+1] = well_roi_array[year]
            well_area_array[year+1] = well_area_array[year]
            well_length_array[year+1] = well_length_array[year]
            
    ##################### annual costs and unit costs ######################### 
    pumping_years = year + 1
       
    # COST SCENARIO 1: reduced capacity and no added wells
    if INCREASE_WELL_DENSITY == False:
        # initialize cost arrays to track annual non-energy costs 
        capital_cost_array = np.zeros((1, int(NUM_YEARS + WELL_LIFETIME)))
        maintenance_array = np.zeros((1, int(NUM_YEARS + WELL_LIFETIME)))
        
        # A) no deepening, initial_sat_thickness < MAXIMUM_INITIAL_SAT_THICKNESS
        if initial_sat_thickness <= MAXIMUM_INITIAL_SAT_THICKNESS:
            install_cost = well_unit_cost * well_length_array[0] # if no deepening, well install remains fixed 
            for year in range(pumping_years):
                capital_cost_array[0, year] = num_wells[0] * install_cost * ((1 + INTEREST_RATE) ** WELL_LIFETIME) * INTEREST_RATE/((1 + INTEREST_RATE) ** WELL_LIFETIME-1)
                maintenance_array[0, year] = MAINTENANCE_RATE * install_cost * num_wells[0] # maintenance cost [% of initial cost]
    
        # B) deepening, initial_sat_thickness > MAXIMUM_INITIAL_SAT_THICKNESS
    
        else:
            for year in range(pumping_years):
                
                if year == 0: 
                    install_cost = well_unit_cost * well_length_array[0] 
                    capital_cost_array[0, year] = num_wells[0] * install_cost * ((1 + INTEREST_RATE) ** WELL_LIFETIME) * INTEREST_RATE/((1 + INTEREST_RATE) ** WELL_LIFETIME-1)
                    maintenance_array[0, year] = MAINTENANCE_RATE * install_cost * num_wells[0] # maintenance cost [% of initial cost]
                    
                elif (year+1) % WELL_LIFETIME == 0: # Replace well every n years (well lifetime)
                        
                    install_cost = well_unit_cost * well_length_array[year]
                    capital_cost_array[0, year] += num_wells[0] * install_cost * ((1 + INTEREST_RATE) ** WELL_LIFETIME) * INTEREST_RATE/((1 + INTEREST_RATE) ** WELL_LIFETIME-1)
                    maintenance_array[0, year] += MAINTENANCE_RATE * install_cost * num_wells[0] # maintenance cost [% of initial cost]
            
                elif well_length_array[year] - well_length_array[year - 1] > 0:
                    capital_cost_array[0, year] += num_wells[0] * install_cost * ((1 + INTEREST_RATE) ** WELL_LIFETIME) * INTEREST_RATE/((1 + INTEREST_RATE) ** WELL_LIFETIME-1)
                    capital_cost_array[0, year: int(year + WELL_LIFETIME)] += well_unit_cost * (well_length_array[year] - well_length_array[year - 1]) * ((1 + INTEREST_RATE) ** WELL_LIFETIME) * INTEREST_RATE/((1 + INTEREST_RATE) ** WELL_LIFETIME-1) * num_wells[0]
                    install_cost = well_unit_cost * well_length_array[year]
                    maintenance_array[0, year] += MAINTENANCE_RATE * install_cost * num_wells[0]
                    
                else:
                    capital_cost_array[0, year] += num_wells[0] * install_cost * ((1 + INTEREST_RATE) ** WELL_LIFETIME) * INTEREST_RATE/((1 + INTEREST_RATE) ** WELL_LIFETIME-1)
                    maintenance_array[0, year] += MAINTENANCE_RATE * install_cost * num_wells[0] # maintenance cost [% of initial cost]
        
        
    # COST SCENARIO 2: reduced capacity and added wells 
    else:
        # find indexes of years when number of wells increase due to pumping rate reduction 
        # along with pumping rate and corresponding number of wells
        
        well_count = np.unique(num_wells)
        if min(well_count) == 0:
            well_count = np.delete(well_count, 0)
            
        added_well_count = np.zeros(len(well_count))
        for i in range(len(added_well_count)):
            if i == 0:
                added_well_count[i] = well_count[i]
            else:
                added_well_count[i] = well_count[i] - well_count[i-1]
        
        Q_vals = np.unique(Well_Q_array)
        if min(Q_vals) == 0:
            Q_vals = np.delete(Q_vals, 0)
        Q_vals = np.sort(Q_vals)    
        Q_vals = Q_vals[::-1]
            
        Start_indx = np.zeros(len(Q_vals)) # indexes where pumping rate and well num changes 
        if len(Start_indx) == 1:
            pass
        
        else:
            for i in range(pumping_years):
                if i == 0:
                    counter = 1
                    continue 
                if num_wells[i] - num_wells[i-1] > 0:
                    Start_indx[counter] = int(i)
                    counter += 1 
        
        # initialize cost arrays to track annual non-energy costs for each group of added wells 
        capital_cost_array = np.zeros((len(Start_indx), int(NUM_YEARS + WELL_LIFETIME)))
        maintenance_array = np.zeros((len(Start_indx), int(NUM_YEARS + WELL_LIFETIME)))
        
        # A) no deepening, initial_sat_thickness < MAXIMUM_INITIAL_SAT_THICKNESS
        if initial_sat_thickness < MAXIMUM_INITIAL_SAT_THICKNESS:
            install_cost = well_unit_cost * well_length_array[0] # if no deepening, well install remains fixed 
            for added_wells in range(len(added_well_count)):
                offset = int(Start_indx[added_wells]) 
                for year in range(pumping_years):
                    capital_cost_array[added_wells, year + offset] = added_well_count[added_wells] * install_cost * ((1 + INTEREST_RATE) ** WELL_LIFETIME) * INTEREST_RATE/((1 + INTEREST_RATE) ** WELL_LIFETIME-1)
                    maintenance_array[added_wells, year + offset] = MAINTENANCE_RATE * install_cost * added_well_count[added_wells] # maintenance cost [% of initial cost]
    
        # B) deepening, initial_sat_thickness > MAXIMUM_INITIAL_SAT_THICKNESS
        else:
            for added_wells in range(len(added_well_count)):
                offset = int(Start_indx[added_wells]) 
                for year in range(pumping_years):
                    if year + offset == pumping_years:
                        break
                    
                    elif year == 0: 
                        install_cost = well_unit_cost * well_length_array[0] 
                        capital_cost_array[added_wells, year + offset] = added_well_count[added_wells] * install_cost * ((1 + INTEREST_RATE) ** WELL_LIFETIME) * INTEREST_RATE/((1 + INTEREST_RATE) ** WELL_LIFETIME-1)
                        maintenance_array[added_wells, year + offset] = MAINTENANCE_RATE * install_cost * added_well_count[added_wells] # maintenance cost [% of initial cost]
                        
                    elif (year+1) % WELL_LIFETIME == 0: # Replace well every n years (well lifetime), if reduced yeild, cheaper unit cost at 200 gpm and below
                            
                        install_cost = well_unit_cost * well_length_array[year + offset]
                        capital_cost_array[added_wells, year + offset] += added_well_count[added_wells] * install_cost * ((1 + INTEREST_RATE) ** WELL_LIFETIME) * INTEREST_RATE/((1 + INTEREST_RATE) ** WELL_LIFETIME-1)
                        maintenance_array[added_wells, year + offset] += MAINTENANCE_RATE * install_cost * added_well_count[added_wells] # maintenance cost [% of initial cost]
                
                    elif well_length_array[year + offset] - well_length_array[year - 1 + offset] > 0:
                        capital_cost_array[added_wells, year + offset] += added_well_count[added_wells] * install_cost * ((1 + INTEREST_RATE) ** WELL_LIFETIME) * INTEREST_RATE/((1 + INTEREST_RATE) ** WELL_LIFETIME-1)
                        capital_cost_array[added_wells, (year + offset): int((year + offset + WELL_LIFETIME))] += well_unit_cost * (well_length_array[year + offset] - well_length_array[year - 1 + offset]) * ((1 + INTEREST_RATE) ** WELL_LIFETIME) * INTEREST_RATE/((1 + INTEREST_RATE) ** WELL_LIFETIME-1) * added_well_count[added_wells]
                        install_cost = well_unit_cost * well_length_array[year + offset]
                        maintenance_array[added_wells, year + offset] += MAINTENANCE_RATE * install_cost * added_well_count[added_wells]
                        
                    else:
                        capital_cost_array[added_wells, year + offset] += added_well_count[added_wells] * install_cost * ((1 + INTEREST_RATE) ** WELL_LIFETIME) * INTEREST_RATE/((1 + INTEREST_RATE) ** WELL_LIFETIME-1)
                        maintenance_array[added_wells, year + offset] += MAINTENANCE_RATE * install_cost * added_well_count[added_wells] # maintenance cost [% of initial cost]
            
                    
    ####################### annual cost metrics ###############################
    annual_capital_cost = np.zeros(NUM_YEARS)
    maintenance_cost = np.zeros(NUM_YEARS)
    #well_installation_cost = np.zeros(NUM_YEARS)
    nonenergy_cost = np.zeros(NUM_YEARS)
    power = np.zeros(NUM_YEARS)
    energy = np.zeros(NUM_YEARS)
    energy_cost_rate = np.zeros(NUM_YEARS)
    energy_cost = np.zeros(NUM_YEARS)
    total_cost_per_well = np.zeros(NUM_YEARS)
    total_cost_all_wells = np.zeros(NUM_YEARS)
    unit_cost = np.zeros(NUM_YEARS)
    unit_cost_per_km3 = np.zeros(NUM_YEARS)
    unit_cost_per_acreft = np.zeros(NUM_YEARS)
    
    annual_capital_cost = np.sum(capital_cost_array, axis = 0)
    maintenance_cost = np.sum(maintenance_array, axis = 0)
    
    for year in range(pumping_years):
        #well_installation_cost[year] = well_unit_cost * well_length_array[year]
        nonenergy_cost[year] = annual_capital_cost[year] + maintenance_cost[year]
        power[year] = num_wells[year] * (SPECIFIC_WEIGHT * total_head[year] * Well_Q_array[year])/(EFFICIENCY*1000) # kW per hour of operation for all wells
        energy[year] = power[year] * (DAYS * 24) # kWh/year all wells 
        energy_cost_rate[year] = ELECTRICITY_RATE # $ per kWh
        energy_cost[year] = energy[year] * energy_cost_rate[year] # energy cost $/year all wells 
        total_cost_per_well[year] = (nonenergy_cost[year] + energy_cost[year])/ num_wells[year]
        total_cost_all_wells[year] = energy_cost[year] + nonenergy_cost[year]
        
        unit_cost[year] = total_cost_all_wells[year]/volume_all_wells[year] # $/m^3
        unit_cost_per_km3[year] = unit_cost[year] * 10**9 # $/km^3
        unit_cost_per_acreft[year] = unit_cost[year] * 1233.48 # $/acft
    
    pumping_rate_decline = np.zeros(len(Well_Q_array)) # fraction of initial pumping rate for each year, for reduced capacity and no increase well density 
    
    for year in range(pumping_years):
        if year == 0: 
            pumping_rate_decline[year] = 1
        else:
            pumping_rate_decline[year] = Well_Q_array[year]/Well_Q_array[0]
            
    return unit_cost, cumulative_volume_all_wells/(10**9), DTW_array,  pumping_rate_decline, volume_all_wells/(10**9), total_cost_per_well, well_area_array, annual_capital_cost, unit_cost_per_acreft, total_head, Well_Q_array


