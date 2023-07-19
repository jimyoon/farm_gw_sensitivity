# -*- coding: utf-8 -*-

import numpy as np 
import pandas as pd 
import math 
import os 

os.chdir('C:\\Users\\yoon644\\OneDrive - PNNL\\Documents\\PyProjects\\farmer_gw_archetypes')
W_lookup = pd.read_csv('Theis_well_function_table.csv', header = "infer") 

#W_lookup = pd.read_csv('./data_inputs/Theis_well_function_table.csv', header = "infer") 
lookup_idx = pd.Index(W_lookup.W)

#%% Construct aquifer properties parameter space  

# Read in parameters using command line arguments 
S = 0.22       # Storativity = Specific yield 
m = 141.5   # Aquifer thickness/depth [m]
K = 1         # K value [m/d]
WL = 14       # Initial depth to water [m]
R = 0          # Annual recharge in cm/yr
IrrDepth = 12  # Annual irrigation depth [inches]
years = 100     # Years of pumping 
extended = 'reduced_capacity'
cost = .12    # electricity cost $/per kWh

def Analytical(S, m, K, WL, R, IrrDepth, years, extended, cost):
    iteration = 0 # set iteration to zero for first loop
    
    Q = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000] # Well Q [gpm]
    
    # Build dataframe of aquifer scenarios (aquifer_df) for Theis calculations
    m_arr = np.zeros(1*1*1)
    K_arr = np.zeros(1*1*1)
    WL_arr = np.zeros(1*1*1)
    for i in range(1):
        m_arr[1*1*i:1*1*(i+1)] = np.tile(m,1*1)
    
    K_tile = np.zeros(1*1)
    for i in range(1):
        K_tile[i*(1):(i+1)*(1)] = np.ones(1)*K
    
    K_arr = np.tile(K_tile, 1)
    
    WL_tile = np.tile(WL, 1)
    WL_arr = np.tile(WL_tile, 1)
    
    S_arr = np.zeros(1*1)
    for i in range(1):
        S_arr[i*1:(i+1)*1] = S*np.ones(1)
        if i == 1-1:
            m_arr = np.tile(m_arr,1)
            K_arr = np.tile(K_arr,1)
            WL_arr = np.tile(WL_arr,1)
            
    Q_arr = np.zeros(1*len(Q))
    for i in range(len(Q)):
        Q_arr[i*1:(i+1)*1] = Q[i]*np.ones(1)
        if i == len(Q)-1:
            m_arr = np.tile(m_arr,len(Q))
            K_arr = np.tile(K_arr,len(Q))
            WL_arr = np.tile(WL_arr,len(Q))
            S_arr = np.tile(S_arr,len(Q))
            
    R_arr = np.zeros(len(m_arr)*1)
    for i in range(1):
        R_arr[i*len(m_arr):(i+1)*len(m_arr)] = R*np.ones(len(m_arr))
        if i == 1-1:
            m_arr = np.tile(m_arr,1)
            K_arr = np.tile(K_arr,1)
            WL_arr = np.tile(WL_arr,1)
            S_arr = np.tile(S_arr, 1)
            Q_arr = np.tile(Q_arr, 1)
    
    d = {'m': m_arr, 'K': K_arr, 'WL': WL_arr, 'S': S_arr,
         'Trans': np.zeros(len(m_arr)), 'R' : R_arr, 
         'Q': Q_arr, 's': np.zeros(len(m_arr)), 'Irr_D_24hr': np.zeros(len(m_arr)),
         'PWL':np.zeros(len(m_arr)), 'Area':np.zeros(len(m_arr)), 
         'Days': np.zeros(len(m_arr))}
    
    aquifer_df = pd.DataFrame(data = d)
    
    # Calculate Transmisivity (meters^2/year)
    for i in range(aquifer_df.m.size):
        aquifer_df.Trans[i] = 365*aquifer_df.K[i]*(aquifer_df.m[i]-aquifer_df.WL[i]) 
    
    #%% Calculate seasonal irrigation depth based on Q
    
    area = [30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160] # well area [acres]
    #days = [30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100] # days of irrigation per year 
    days = np.linspace(30, 100, int((100-30)/2)+1, endpoint = True)
    area_arr = np.zeros(len(area)*len(days))
    
    for i in range(len(area)):
        area_arr[i*len(days):(i+1)*len(days)] = area[i]*np.ones(len(days))
        
    days_arr = np.zeros(len(area)*len(days))
    days_arr = np.tile(days, len(area))
        
    data = {'area': area_arr, 'days': days_arr, 'Q_100_Annual': np.zeros(len(area)*len(days)), 
            'Q_200_Annual': np.zeros(len(area)*len(days)), 'Q_300_Annual': np.zeros(len(area)*len(days)),
            'Q_400_Annual': np.zeros(len(area)*len(days)), 'Q_500_Annual': np.zeros(len(area)*len(days)), 
            'Q_600_Annual': np.zeros(len(area)*len(days)), 'Q_700_Annual': np.zeros(len(area)*len(days)), 
            'Q_800_Annual': np.zeros(len(area)*len(days)), 'Q_900_Annual': np.zeros(len(area)*len(days)), 
            'Q_1000_Annual': np.zeros(len(area)*len(days)),
            'Q_100_24h': np.zeros(len(area)*len(days)), 'Q_200_24h': np.zeros(len(area)*len(days)),
            'Q_300_24h': np.zeros(len(area)*len(days)), 'Q_400_24h': np.zeros(len(area)*len(days)),
            'Q_500_24h': np.zeros(len(area)*len(days)), 'Q_600_24h': np.zeros(len(area)*len(days)),
            'Q_700_24h': np.zeros(len(area)*len(days)), 'Q_800_24h': np.zeros(len(area)*len(days)),
            'Q_900_24h': np.zeros(len(area)*len(days)), 'Q_1000_24h': np.zeros(len(area)*len(days))}
    
    irr_depth_df = pd.DataFrame(data = data)
        
    for i in range(irr_depth_df.iloc[:,0].size):
        for j in range(len(Q)):
            irr_depth_df.iloc[i,j+2] = (Q[j]*1440*irr_depth_df.days[i]*1/264)/(irr_depth_df.area[i]*4046.86) # irrigation depth in meters, conversion of gpm to m^3 and area in acres to m^2
            irr_depth_df.iloc[i,j+2+len(Q)] = (Q[j]*1440*1/264)/(irr_depth_df.area[i]*4046.86) # 24 hr irrigation depth 
            
    #%% For each aquifer scenario, determine possible (1) Area irrigated AND required (2) Days of irrigation 
      # to meet defined annual Irrigation Depth (IrrDepth)
     
    Irr_depth_m = .0256*IrrDepth # irrigation depth in m, provided irr depth [in] 
    irr_depth_df_Q_ilocs = np.array([[100, 200, 300, 400, 500, 600, 700, 800, 900, 1000],[2,3,4,5,6,7,8,9,10,11]])
    
    for i in range(aquifer_df.m.size):
        irr_Q = int(aquifer_df.Q[i])
        indx = np.where(irr_depth_df_Q_ilocs[0,:] == irr_Q)
        col = irr_depth_df_Q_ilocs[1,indx[0]]
        temp = np.zeros((irr_depth_df.iloc[:,0].size, 2))
        for j in range(irr_depth_df.iloc[:,0].size):
            if abs(irr_depth_df.iloc[j,col[0]]-Irr_depth_m) < 0.05:
                temp[j,:] = irr_depth_df.iloc[j,0:2].values
            if j+1 == irr_depth_df.iloc[:,0].size:
                temp_sorted = np.argsort(temp[:,0])
                aquifer_df.iloc[i,10] = temp[temp_sorted[-1],0] # agents want to maximize area served 
                aquifer_df.iloc[i,11]  = np.min(temp[np.where(temp[:,0] == aquifer_df.iloc[i,10]),1]) # NEW minimizes days or irrigation
                #aquifer_df.iloc[i,11] = temp[temp_sorted[-1],1] # days of irrigation for area served OLD
    
      
    #%% One year of simulated aquifer response at half day resolution (time steps = 0.5 d or 1/730 yr)

    for iteration in range(len(aquifer_df.iloc[:,0])):
        #print(iteration)
        if iteration == 0: 
            
            ###### Code block for first iteration with preferred (selected) initial pumping rate 
            
            T_array = np.zeros((aquifer_df.iloc[:,0].size, years+1))
            D_array = np.zeros((aquifer_df.iloc[:,0].size, years+1))
            Days_array  = np.zeros((aquifer_df.iloc[:,0].size, years+1))
            
            # Starting Transmissivity not to exceed 100 m of saturated thickness
            if aquifer_df.m[0] - aquifer_df.WL[0] > 100:
                T_array[:,0] = aquifer_df.K[0] * 100 * 365 # T [m^2/yr]  based on first 100 m of sat thickness 
                D_array[:,0] = aquifer_df.WL[0] + 100 # initial depth of well = depth to water (WL) + 100 m sat thickness                                           
                
            else: 
                T_array[:,0] = aquifer_df.Trans # T [m^2/yr] # if sat thickness is <100 m then use T from aquifer props as initial 
                D_array[:,0] = aquifer_df.m[0] 
            
            One_yr_ts = np.zeros((aquifer_df.iloc[:,0].size, 730*1))
            WL_array = np.zeros((aquifer_df.iloc[:,0].size, years+1))
            WL_array[:,0] = aquifer_df.WL
            
            # One year of pumping for each aquifer scenario to determine viability of 
            # pumping rate (cannot violate constraints on total drawdown, drawdown as 
            # fraction of saturated thickness, or well screen entrance velocity)
            
            for x in range(aquifer_df.iloc[:,0].size): # number of aquifer scenarios (rows in aquifer_df)
                for t in range(1): # 1 year of simulated pumping to determine drawdown for each pumping rate under consideration
                    if T_array[x,t] > 0: 
                        row = x
                        Q_gpm = aquifer_df.Q[row]
                        Q = Q_gpm*1440*1/264*365 # pumping rate in GPM converted, GPD (mult by 1440), to m^3/day (divide by 264), to m^3/year
                        Qv = Q_gpm*1440*1/264*aquifer_df.Days[row] # pumped volume (m^3) per yr 
                        r = 0.5 # radial distance of drawdown estimate
                        well_area = aquifer_df.Area[row]*4046.86 # WELL AREA SERVED BASED ON WELL CAPACITY
                    
                        # Annual pump pattern 
                        pump_pattern = np.array([1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
                        pump_ts = np.tile(pump_pattern, int(aquifer_df.Days[row]/2))
                        off_period = np.zeros(730 - pump_ts.size)
                        pump_ts = np.concatenate((pump_ts,off_period))
                        duration = len(pump_ts)
                        
                        # Find indexes for start times of pumping and injection well 
                        pump_indexes = np.asarray(np.where(pump_ts > 0))
                        pump_start_times = []
                        for i in range(pump_indexes.size-1):
                            if i == 0 :
                                pump_start_times.append(pump_indexes[0,i])
                            elif pump_indexes[0,i]-pump_indexes[0,i+1] < -1:
                               pump_start_times.append(pump_indexes[0,i+1])
                                    
                        pump_durations = duration-np.asarray(pump_start_times)
                        pump_array = np.zeros((len(pump_start_times),duration+1))
                        
                        inject_indexes = np.asarray(np.where(pump_ts < 1))
                        inject_start_times = []
                        for i in range(inject_indexes.size-1):
                            if i == 0 :
                                inject_start_times.append(inject_indexes[0,i])
                            
                            elif inject_indexes[0,i]-inject_indexes[0,i+1] < -1:
                                inject_start_times.append(inject_indexes[0,i+1])
                        
                        inject_durations = duration-np.asarray(inject_start_times)
                        inject_array = np.zeros((len(inject_start_times),duration+1))
                        
                        for i in range(len(pump_start_times)):
                            pump_time = np.arange((pump_durations[i]))/730+1/730
                            start_time = pump_start_times[i]
                            for j in range((pump_durations[i])):
                                u = ((r**2)*aquifer_df.S[row])/(4*T_array[x,t]*pump_time[j])
                                if u > 0.0099:
                                    lookup_idx = pd.Index(W_lookup.u)
                                    lookup_loc = lookup_idx.get_loc(u, method = 'nearest')
                                    W = W_lookup.W[lookup_loc]
                                    pump_array[i,j+1+start_time] = (Q*W)/(4*math.pi*T_array[x,t]) 
                                else:
                                    W = -0.577-math.log(u)+u-(u**2)/4+u**3/(3*6)
                                    pump_array[i,j+1+start_time] = (Q*W)/(4*math.pi*T_array[x,t]) 
                        
                        for i in range(len(inject_start_times)):
                            inject_time = np.arange((inject_durations[i]))/730+1/730
                            start_time = inject_start_times[i]
                            for j in range((inject_durations[i])):
                                u = ((r**2)*aquifer_df.S[row])/(4*T_array[x,t]*inject_time[j])
                                if u > 0.0099:
                                    lookup_idx = pd.Index(W_lookup.u)
                                    lookup_loc = lookup_idx.get_loc(u, method = 'nearest')
                                    W = W_lookup.W[lookup_loc]
                                    inject_array[i,j+1+start_time] = -1*(Q*W)/(4*math.pi*T_array[x,t]) 
                                else:
                                    W = -0.577-math.log(u)+u-(u**2)/4+u**3/(3*6)
                                    inject_array[i,j+1+start_time] = -1*(Q*W)/(4*math.pi*T_array[x,t]) 
                        
                        superposition_arr = np.concatenate((pump_array, inject_array), axis = 0)
                        drawdown_total = -1*np.sum(superposition_arr, axis = 0) - WL_array[x,t]
                        T_array[x,t+1] = T_array[x,t] - aquifer_df.K[x]*365*((Qv/(well_area))/aquifer_df.S[row]) # update T value due to drop in water table from depletion 
                        #WL_array[x,t+1] = WL_array[x,t] + ((Qv/well_area)-aquifer_df.R[x]/100)/aquifer_df.S[row]
                        One_yr_ts[x,(t*730):((t+1)*730)] = drawdown_total[1:] 
                    else:
                        continue
                    
            # Calculate max drawdown over 1 yr test period  
            for i in range(aquifer_df.iloc[:,0].size):
                 if aquifer_df.Trans[i] == 0:
                     aquifer_df.s[i] = -999
                 else:
                     aquifer_df.PWL[i] = -1*min(One_yr_ts[i,:]) # max depth to water during 1 year of cyclical pumping (m below surface)
                     aquifer_df.s[i] = -1*min(One_yr_ts[i,:])-aquifer_df.WL[i] # max drawdown during 1 year of cyclical of pumping (m)
              
            # Filter out non-viable scenarios
            aquifer_df['viable'] = np.ones(aquifer_df.iloc[:,0].size)
            for i in range(aquifer_df.iloc[:,0].size):
                if aquifer_df.m[0] - aquifer_df.WL[0] > 100: # for cases where initial aquifer thickness is > 100 m 
                    if aquifer_df.s[i]/100 > 0.6 or aquifer_df.Trans[i] == 0 or aquifer_df.s[i] > 50 or aquifer_df.Area[i] == 0:
                        aquifer_df.viable[i] = 0
                else: # for cases where initial aquifer saturated thickness is < 100 m 
                    if aquifer_df.s[i]/(aquifer_df.m[i] - aquifer_df.WL[i]) > 0.6 or aquifer_df.Trans[i] == 0 or aquifer_df.s[i] > 50 or aquifer_df.Area[i] == 0:
                        aquifer_df.viable[i] = 0
                         
            # Choose preferrable viable aquifer scenario based on Area and Days  
            viable_scenarios = aquifer_df.copy()
            viable_scenarios.where(aquifer_df.viable > 0, inplace = True)
            viable_scenarios = viable_scenarios.dropna()
            
            if np.sum(aquifer_df.viable[:].values) == 0: # exit if no viable scenarios 
                print('Non-viable')
                viable = 0 
        
            else: # do the following if viable (at least one viable pumping rate)
                viable = 1
                viable_first_Q = 1
                Area_indx = viable_scenarios.Area.sort_values() 
                Area_indx = Area_indx.where(Area_indx == np.max(Area_indx))
                Area_indx.dropna(inplace = True)
                Days_count = viable_scenarios.Days[Area_indx.index].sort_values()
                #Days_count = viable_scenarios.Days[Area_indx[0][:]].sort_values()
                row = int(Days_count.index[0]) # row of aquifer_df scenario to use for 100 yr simulation 
                initial_row = int(row)
                
            # Build arrays to track pumping rate, well area, well transmissivity, and
            # well depth for each simulated year
            Q_array = np.zeros((aquifer_df.iloc[:,0].size, years+1)) # tracks Q through time
            A_array = np.zeros((aquifer_df.iloc[:,0].size, years+1)) # tracks Area through time
            T_array = np.zeros((aquifer_df.iloc[:,0].size, years+1)) # reset after aquifer viability pumping evalutation, tracks Transmissivity through time
            D_array = np.zeros((aquifer_df.iloc[:,0].size, years+1)) # tracks well depth through time 
            
            # Starting Transmissivity not to exceed 100 m of saturated thickness
            if aquifer_df.m[0] - aquifer_df.WL[0] > 100: # if initial saturated thickness is greater than 100 m 
                T_array[:,0] = aquifer_df.K[0] * 100 * 365 # T [m^2/yr]  based on first 100 m of sat thickness 
                D_array[:,0] = aquifer_df.WL[0] + 100 # initial depth of well = depth to water (WL) + 100 m sat thickness                                           
                
            else: # if initial saturated thickness is less than 100 m 
                T_array[:,0] = aquifer_df.Trans # T [m^2/yr] # if sat thickness is <100 m then use T from aquifer props as initial 
                D_array[:,0] = aquifer_df.m[0] 
            
            # Code block simulating pumping at Half day resolution up to the years specified by 
            # the parameter 'years'. If the aquifer saturated thickness is > 100 m, the loop 
            # will deepen the well in 30 m increments to enable continued pumping at initial 
            # selected pumping rate. The loop terminates when one of the constraints is 
            # violated and the well cannot be deepened to accomodate. This triggers the 
            # reduced pumping loop that follows in the second code block. 
            # The reduced pumping code block (Line 420) reassess pumping rate viability 
            # and determines new pumping rate to continue cost curve simulation. If even the 
            # lowest pumping rate violates the constriant the pumping simulation is terminated. 
            
            DTW_ts = np.zeros((730*years))  # depth to water at half-day resolution. (time steps = 0.5 d or 1/730 yr)
            WL_array = np.zeros((aquifer_df.iloc[:,0].size, years+1)) # tracks starting water level at beggining of each annual period 
            WL_array[:,0] = aquifer_df.WL # initial WL for year 1 
            aquifer_pump_indexes = [] # list of lists where each aquifer scenario has a record of pumping indexes for PWL calcs
            entrance_V = (aquifer_df.Q[row]/(7.48*60))/((D_array[row,0]-aquifer_df.WL[row])*3*2*math.pi*8/12*0.2) # casing entrance velocity, potential constrain, but not currently implemented 
            max_drawdown = aquifer_df.s[row] # max drawdown for first year, based on selected viable pumping rate 
            sat_thickness =  (D_array[row,0]-aquifer_df.WL[row]) # initial saturated thickness from initial 
                                                                 # aquifer depth - starting water depth not to exceed 100 m sat thickness
            previous_row = int(row)
            
            for t in range(years): 
                if T_array[row,t] > 0 and entrance_V < 0.25 and max_drawdown/sat_thickness < 0.6 and viable > 0: # conditionals to determine if pumping is possible 
                    Q_gpm = aquifer_df.Q[row]
                    Q_array[row,t] = aquifer_df.Q[row]
                    A_array[row,t] = aquifer_df.Area[row]
                    Days_array[row,t] = aquifer_df.Days[row]
                    Q = Q_gpm*1440*1/264*365 # pumping rate in GPM converted, GPD (mult by 1440), to m^3/day (divide by 264), to m^3/year
                    Qv = Q_gpm*1440*1/264*aquifer_df.Days[row] # pumped volume (m^3) per yr 
                    r = 0.5 # radial distance of drawdown estimate
                    well_area = aquifer_df.Area[row]*4046.86 # Area served by well capacity (convert from acres to m^2)
                
                    # Annual pump pattern 
                    pump_pattern = np.array([1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
                    pump_ts = np.tile(pump_pattern, int(aquifer_df.Days[row]/2)) # number of weekly cycles based on days of irrigation
                    off_period = np.zeros(730 - pump_ts.size)
                    pump_ts = np.concatenate((pump_ts,off_period))
                    duration = len(pump_ts)
                    
                    # Find indexes for start times of pumping and injection well
                    pump_indexes = np.asarray(np.where(pump_ts > 0))
                    if t == 0:
                        aquifer_pump_indexes.append(pump_indexes) # UDDATE to reference new list location for each x iteration
                    
                    pump_start_times = []
                    for i in range(pump_indexes.size-1):
                        if i == 0 :
                            pump_start_times.append(pump_indexes[0,i])
                        elif pump_indexes[0,i]-pump_indexes[0,i+1] < -1:
                           pump_start_times.append(pump_indexes[0,i+1])
                                
                    pump_durations = duration-np.asarray(pump_start_times)
                    pump_array = np.zeros((len(pump_start_times),duration+1))
                    
                    inject_indexes = np.asarray(np.where(pump_ts < 1))
                    inject_start_times = []
                    for i in range(inject_indexes.size-1):
                        if i == 0 :
                            inject_start_times.append(inject_indexes[0,i])
                        
                        elif inject_indexes[0,i]-inject_indexes[0,i+1] < -1:
                            inject_start_times.append(inject_indexes[0,i+1])
                    
                    inject_durations = duration-np.asarray(inject_start_times)
                    inject_array = np.zeros((len(inject_start_times),duration+1))
                    
                    for i in range(len(pump_start_times)):
                        pump_time = np.arange((pump_durations[i]))/730+1/730 # pump time steps in fractions of a year 
                        start_time = pump_start_times[i]
                        for j in range((pump_durations[i])):
                            u = ((r**2)*aquifer_df.S[row])/(4*T_array[row,t]*pump_time[j])
                            if u > 0.0099:
                                lookup_idx = pd.Index(W_lookup.u)
                                lookup_loc = lookup_idx.get_loc(u, method = 'nearest')
                                W = W_lookup.W[lookup_loc]
                                pump_array[i,j+1+start_time] = (Q*W)/(4*math.pi*T_array[row,t]) 
                            else:
                                W = -0.577-math.log(u)+u-(u**2)/4+u**3/(3*6)
                                pump_array[i,j+1+start_time] = (Q*W)/(4*math.pi*T_array[row,t]) 
                    
                    for i in range(len(inject_start_times)):
                        inject_time = np.arange((inject_durations[i]))/730+1/730
                        start_time = inject_start_times[i]
                        for j in range((inject_durations[i])):
                            u = ((r**2)*aquifer_df.S[row])/(4*T_array[row,t]*inject_time[j])
                            if u > 0.0099:
                                lookup_idx = pd.Index(W_lookup.u)
                                lookup_loc = lookup_idx.get_loc(u, method = 'nearest')
                                W = W_lookup.W[lookup_loc]
                                inject_array[i,j+1+start_time] = -1*(Q*W)/(4*math.pi*T_array[row,t]) 
                            else:
                                W = -0.577-math.log(u)+u-(u**2)/4+u**3/(3*6)
                                inject_array[i,j+1+start_time] = -1*(Q*W)/(4*math.pi*T_array[row,t])
                                
                    superposition_arr = np.concatenate((pump_array, inject_array), axis = 0)
                    drawdown_total = -1*np.sum(superposition_arr, axis = 0) - WL_array[row,t] # 
                    T_array[row,t+1] = T_array[row,t] - aquifer_df.K[row]*365*((Qv/(well_area)-aquifer_df.R[row]/100)/aquifer_df.S[row]) # update T value due to drop in water table from depletion 
                    WL_array[row,t+1] = WL_array[row,t] + ((Qv/well_area)-aquifer_df.R[row]/100)/aquifer_df.S[row] # water depth below ground surface 
                    DTW_ts[(t*730):((t+1)*730)] = drawdown_total[1:]
                    max_drawdown = (-1*min(drawdown_total) - WL_array[row,t]) # max drawdown during current simulated year 
                    sat_thickness = D_array[row, t] - WL_array[row, t] # remaining saturated thickness 
                    entrance_V = (aquifer_df.Q[row]/(7.48*60))/(sat_thickness*3*2*math.pi*8/12*0.2) # casing entrance velocity, potential constrain, but not currently implemented 
                    if max_drawdown/sat_thickness > 0.8 and aquifer_df.m[row] > D_array[row,t]:
                        additional_depth_remaining = aquifer_df.m[row] - D_array[row, t]
                        if additional_depth_remaining > 30:
                            D_array[row, t+1] = D_array[row, t] + 30 # deepen well by 30 m 
                            T_array[row,t+1] = T_array[row,t+1] + 30*365*aquifer_df.K[row] # update Transmissivity 
                            sat_thickness = D_array[row, t+1] - WL_array[row, t+1]
                        else:
                            D_array[row, t+1] =  D_array[row, t] + additional_depth_remaining
                            T_array[row,t+1] = T_array[row,t+1] + additional_depth_remaining*365*aquifer_df.K[row] # update Transmissivity
                            sat_thickness = D_array[row, t+1] - WL_array[row, t+1]
                    else:
                        D_array[row, t+1] = D_array[row, t] # well depth remains the same
                    
                else:
                    D_array[row,t] = 0
                    break                         
        
        
            
        ###### Code block for additional iterations where pumping rate is reduced ######
        
        elif (t+1) >= years or row == 0 or viable == 0:
            break
        
        else:
                 
             years_remaining = years - (t+1)
             T_array[0:row,t] = T_array[row,t]
             WL_array[:,t] = WL_array[row,t]
             T_array[previous_row,t+1], WL_array[previous_row, t+1] = 0, 0
             One_yr_ts = np.zeros((aquifer_df.iloc[:,0].size, 730*1))
             
             # TRY  
             D_array[:,t] = max(D_array[:,t-1])
             #Days_array[:,t] = max(Days_array[:,t-1]) 
             
             for x in range(aquifer_df.iloc[:,0].size-(10-row)): # number of aquifer scenarios (rows in aquifer_df)
                 #for t in range(1): # 1 year of simulated pumping with remaining T to determine drawdown for each pumping rate under consideration
                 if T_array[x,t] > 0: 
                     row = x
                     Q_gpm = aquifer_df.Q[row]
                     Q = Q_gpm*1440*(1/264)*365 # pumping rate in GPM converted, GPD (mult by 1440), to m^3/day (divide by 264), to m^3/year
                     Qv = Q_gpm*1440*(1/264)*aquifer_df.Days[row] # pumped volume (m^3) per yr 
                     r = 0.5 # radial distance of drawdown estimate
                     well_area = aquifer_df.Area[row]*4046.86 # WELL AREA SERVED BASED ON WELL CAPACITY
                 
                     # Annual pump pattern 
                     pump_pattern = np.array([1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
                     pump_ts = np.tile(pump_pattern, int(aquifer_df.Days[row]/2))
                     off_period = np.zeros(730 - pump_ts.size)
                     pump_ts = np.concatenate((pump_ts,off_period))
                     duration = len(pump_ts)
                     
                     # Find indexes for start times of pumping and injection well 
                     pump_indexes = np.asarray(np.where(pump_ts > 0))
                     pump_start_times = []
                     for i in range(pump_indexes.size-1):
                         if i == 0 :
                             pump_start_times.append(pump_indexes[0,i])
                         elif pump_indexes[0,i]-pump_indexes[0,i+1] < -1:
                            pump_start_times.append(pump_indexes[0,i+1])
                                 
                     pump_durations = duration-np.asarray(pump_start_times)
                     pump_array = np.zeros((len(pump_start_times),duration+1))
                     
                     inject_indexes = np.asarray(np.where(pump_ts < 1))
                     inject_start_times = []
                     for i in range(inject_indexes.size-1):
                         if i == 0 :
                             inject_start_times.append(inject_indexes[0,i])
                         
                         elif inject_indexes[0,i]-inject_indexes[0,i+1] < -1:
                             inject_start_times.append(inject_indexes[0,i+1])
                     
                     inject_durations = duration-np.asarray(inject_start_times)
                     inject_array = np.zeros((len(inject_start_times),duration+1))
                     
                     for i in range(len(pump_start_times)):
                         pump_time = np.arange((pump_durations[i]))/730+1/730
                         start_time = pump_start_times[i]
                         for j in range((pump_durations[i])):
                             u = ((r**2)*aquifer_df.S[row])/(4*T_array[x,t]*pump_time[j])
                             if u > 0.0099:
                                 lookup_idx = pd.Index(W_lookup.u)
                                 lookup_loc = lookup_idx.get_loc(u, method = 'nearest')
                                 W = W_lookup.W[lookup_loc]
                                 pump_array[i,j+1+start_time] = (Q*W)/(4*math.pi*T_array[x,t]) 
                             else:
                                 W = -0.577-math.log(u)+u-(u**2)/4+u**3/(3*6)
                                 pump_array[i,j+1+start_time] = (Q*W)/(4*math.pi*T_array[x,t]) 
                     
                     for i in range(len(inject_start_times)):
                         inject_time = np.arange((inject_durations[i]))/730+1/730
                         start_time = inject_start_times[i]
                         for j in range((inject_durations[i])):
                             u = ((r**2)*aquifer_df.S[row])/(4*T_array[x,t]*inject_time[j])
                             if u > 0.0099:
                                 lookup_idx = pd.Index(W_lookup.u)
                                 lookup_loc = lookup_idx.get_loc(u, method = 'nearest')
                                 W = W_lookup.W[lookup_loc]
                                 inject_array[i,j+1+start_time] = -1*(Q*W)/(4*math.pi*T_array[x,t]) 
                             else:
                                 W = -0.577-math.log(u)+u-(u**2)/4+u**3/(3*6)
                                 inject_array[i,j+1+start_time] = -1*(Q*W)/(4*math.pi*T_array[x,t]) 
                     
                     superposition_arr = np.concatenate((pump_array, inject_array), axis = 0)
                     drawdown_total = -1*np.sum(superposition_arr, axis = 0) - WL_array[x,t] # drawdown relative to ground surface 
                     
                     #T_array[x,t+1] = T_array[x,t] - aquifer_df.K[x]*365*((Qv/(well_area))/aquifer_df.S[row]) # update T value due to drop in water table from depletion 
                     #WL_array[x,t+1] = WL_array[x,t] + ((Qv/well_area)-aquifer_df.R[x]/100)/aquifer_df.S[row]
                     One_yr_ts[x,(0*730):((0+1)*730)] = drawdown_total[1:] # half day timer series of water depth from surface 
                 else:
                     continue
                         
            # Calculate max drawdown over 1 yr pumping period 
             for i in range(aquifer_df.iloc[:,0].size):
                  if T_array[i,t] == 0:
                      aquifer_df.s[i] = 999
                  else:
                      aquifer_df.PWL[i] = -1*min(One_yr_ts[i,:]) # max depth to water during 1 year of cyclical pumping (m below surface)
                      aquifer_df.s[i] = -1*min(One_yr_ts[i,:]) - WL_array[previous_row, t] # aquifer_df.WL[i] # max drawdown during 1 year of cyclical of pumping (m)
              
             # filter out non-viable scenarios
             aquifer_df['viable'] = np.ones(aquifer_df.iloc[:,0].size) # reset viability column 
             for i in range(aquifer_df.iloc[:,0].size):
                 # if aquifer_df.m[0] - aquifer_df.WL[0] > 100:
                 #     if aquifer_df.s[i]/100 > 0.6 or aquifer_df.Trans[i] == 0 or aquifer_df.s[i] > 50 or aquifer_df.Area[i] == 0:
                 #         aquifer_df.viable[i] = 0
                 #else:
                 if aquifer_df.s[i]/(aquifer_df.m[0]-WL_array[previous_row,t]) > 0.6 or aquifer_df.s[i] > 50 or aquifer_df.Area[i] == 0 or aquifer_df.s[i] < 0:
                     aquifer_df.viable[i] = 0
                         
             # Choose preferrable viable aquifer scenario based on Area and Days
         
             viable_scenarios = aquifer_df.copy()
             viable_scenarios.where(aquifer_df.viable > 0, inplace = True)
             viable_scenarios = viable_scenarios.dropna()
             
             if np.sum(aquifer_df.viable[:].values) == 0: # exit if no viable scenarios 
                 print('Non-viable')
                 viable = 0
             else:
                 viable = 1 # at least one pumping rate is still viable
                 Area_indx = viable_scenarios.Area.sort_values() 
                 Area_indx = Area_indx.where(Area_indx == np.max(Area_indx))
                 Area_indx.dropna(inplace = True)
                 Days_count = viable_scenarios.Days[Area_indx.index].sort_values()
                 #Days_count = viable_scenarios.Days[Area_indx[0][:]].sort_values()
                 row = int(Days_count.index[0]) # row of aquifer_df scenario to use for 100 yr simulation, pumping rate with least required days of pumping to meet annual demand 
                 # Area_indx = viable_scenarios.Area.sort_values() 
                 # Area_indx = np.where(Area_indx == np.max(Area_indx))
                 # viable_scenarios.reset_index
                 # Days_count = viable_scenarios.Days[Area_indx[0][:]].sort_values()
                 # row = int(Days_count.index[0]) # row of aquifer_df scenario to use for 100 yr simulation 
                 
                 # Assigns new area to the well if increased density to offset 
                 # reduced capacity is selected, otherwise the initial well area 
                 # is used for all years of simulation 
                 if extended == 'reduced_capacity':
                     area_row = int(initial_row) # row index for inital well Area for initial well Q
                 else:
                     area_row = int(row)
                 
                 
             entrance_V = (aquifer_df.Q[row]/(7.48*60))/((D_array[previous_row,t]-WL_array[previous_row,t])*3*2*math.pi*(8/2)/12*0.2) # CHECK D_array val on iteration 1 # casing entrance velocity constraint
             #D_array[row, t] = int(D_array[previous_row,t-1])
             D_array[row, t] = float(D_array[previous_row,t])
             A_array[row, t] = int(A_array[previous_row,t])
             Days_array[row, t] = int(Days_array[previous_row,t])
             
             # FOR FIXED AREA CASE, ASSIGN NEW # of PUMPING DAYS UP TO MAX (100)                        
             initial_Qv = aquifer_df.Days[initial_row]*aquifer_df.Q[initial_row] * 1440 * 1/264 # m^3 pumped 
             reduced_cap_days = initial_Qv/(aquifer_df.Q[row] * 1440 * 1/264) # number of days pumping a reduced rate to meet demand 
             reduced_cap_days_diff = days - reduced_cap_days
             reduced_days = 0 
             if reduced_cap_days > 100:
                 reduced_days = 100
             else: 
                 reduced_capacity_indx = 0
                 for i in range(len(reduced_cap_days_diff)):
                     if reduced_cap_days_diff[i] < 0 :
                         continue
                     else: # first days index (non-negative days_diff) that meets annual pumped volume target 
                         reduced_capacity_indx = i
                         reduced_days = days[reduced_capacity_indx] 
                         break 
                    
             # reduced_days = 100 
             # count = 0
             # for i in range(len(days)):
             #     if reduced_cap_days_diff[i] < 0:
             #         count += 1
             #     else:
             #         if reduced_cap_days_diff[i] < reduced_days:
             #             reduced_days = reduced_cap_days_diff[i]
             # if count == len(days):
             #     reduced_days = 100
                 
             # else:
             #     reduced_days = np.where(reduced_cap_days_diff == reduced_days)
             #     reduced_days = days[reduced_days[0][0]]
             
             D_array[previous_row,t], A_array[previous_row,t], Q_array[previous_row,t], Days_array[previous_row,t] = 0, 0, 0, 0 # check D_array on Iteration 1
             #D_array[row,t] = D_array[previous_row, t-1]
             max_drawdown = aquifer_df.s[row] # max drawdown during annual pumped volume viability test (preceeding pumping code block)
             sat_thickness =  (D_array[row,t]-WL_array[row,t]) # NEEDS UPDATING TO REFERENCE PROPER ROW AND COLUMN INDICIES
             previous_row = int(row) # update value to track 'previous_row' for next iteration of reduced well Q 
             start_yr = int(t)
             
             for year in range(years_remaining): # years remaining of pumping to reach prescribed years of pumping
                 if T_array[row,t] > 0 and entrance_V < 0.25 and max_drawdown/sat_thickness < 0.6 and viable > 0: # conditionals to determine if pumping is possible 
                     t = 0
                     t = start_yr + year 
                     Q_gpm = aquifer_df.Q[row]
                     Q_array[row,t] = aquifer_df.Q[row]
                     A_array[row,t] = aquifer_df.Area[area_row]
                     if extended == 'reduced_capacity':
                         Days_array[row,t] = reduced_days
                     else:   
                         Days_array[row,t] = aquifer_df.Days[row]
                         
                     Q = Q_gpm*1440*1/264*365 # pumping rate in GPM converted, GPD (mult by 1440), to m^3/day (divide by 264), to m^3/year
                     #Qv = Q_gpm*1440*1/264*aquifer_df.Days[row] # pumped volume (m^3) per yr 
                     Qv = Q_gpm*1440*1/264*Days_array[row,t] # pumped volume (m^3) per yr 
                     r = 0.5 # radial distance of drawdown estimate
                     well_area = aquifer_df.Area[area_row]*4046.86 # Area served by well capacity (convert from acres to m^2)
                 
                     # Annual pump pattern 
                     pump_pattern = np.array([1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
                     #pump_ts = np.tile(pump_pattern, int(aquifer_df.Days[row]/2)) # number of weekly cycles based on days of irrigation
                     pump_ts = np.tile(pump_pattern, int(Days_array[row,t]/2)) # number of weekly cycles based on days of irrigation
                     off_period = np.zeros(730 - pump_ts.size)
                     pump_ts = np.concatenate((pump_ts,off_period))
                     duration = len(pump_ts)
                     
                     # Find indexes for start times of pumping and injection well
                     pump_indexes = np.asarray(np.where(pump_ts > 0))
                     if t == 0:
                         aquifer_pump_indexes.append(pump_indexes) 
                     
                     pump_start_times = []
                     for i in range(pump_indexes.size-1):
                         if i == 0 :
                             pump_start_times.append(pump_indexes[0,i])
                         elif pump_indexes[0,i]-pump_indexes[0,i+1] < -1:
                            pump_start_times.append(pump_indexes[0,i+1])
                                 
                     pump_durations = duration-np.asarray(pump_start_times)
                     pump_array = np.zeros((len(pump_start_times),duration+1))
                     
                     inject_indexes = np.asarray(np.where(pump_ts < 1))
                     inject_start_times = []
                     for i in range(inject_indexes.size-1):
                         if i == 0 :
                             inject_start_times.append(inject_indexes[0,i])
                         
                         elif inject_indexes[0,i]-inject_indexes[0,i+1] < -1:
                             inject_start_times.append(inject_indexes[0,i+1])
                     
                     inject_durations = duration-np.asarray(inject_start_times)
                     inject_array = np.zeros((len(inject_start_times),duration+1))
                     
                     for i in range(len(pump_start_times)):
                         pump_time = np.arange((pump_durations[i]))/730+1/730 # pump time steps in fractions of a year 
                         start_time = pump_start_times[i]
                         for j in range((pump_durations[i])):
                             u = ((r**2)*aquifer_df.S[row])/(4*T_array[row,t]*pump_time[j])
                             if u > 0.0099:
                                 lookup_idx = pd.Index(W_lookup.u)
                                 lookup_loc = lookup_idx.get_loc(u, method = 'nearest')
                                 W = W_lookup.W[lookup_loc]
                                 pump_array[i,j+1+start_time] = (Q*W)/(4*math.pi*T_array[row,t]) 
                             else:
                                 W = -0.577-math.log(u)+u-(u**2)/4+u**3/(3*6)
                                 pump_array[i,j+1+start_time] = (Q*W)/(4*math.pi*T_array[row,t]) 
                     
                     for i in range(len(inject_start_times)):
                         inject_time = np.arange((inject_durations[i]))/730+1/730
                         start_time = inject_start_times[i]
                         for j in range((inject_durations[i])):
                             u = ((r**2)*aquifer_df.S[row])/(4*T_array[row,t]*inject_time[j])
                             if u > 0.0099:
                                 lookup_idx = pd.Index(W_lookup.u)
                                 lookup_loc = lookup_idx.get_loc(u, method = 'nearest')
                                 W = W_lookup.W[lookup_loc]
                                 inject_array[i,j+1+start_time] = -1*(Q*W)/(4*math.pi*T_array[row,t]) 
                             else:
                                 W = -0.577-math.log(u)+u-(u**2)/4+u**3/(3*6)
                                 inject_array[i,j+1+start_time] = -1*(Q*W)/(4*math.pi*T_array[row,t])
                                 
                     superposition_arr = np.concatenate((pump_array, inject_array), axis = 0)
                     drawdown_total = -1*np.sum(superposition_arr, axis = 0) - WL_array[row,t] # 
                     T_array[row,t+1] = T_array[row,t] - aquifer_df.K[row]*365*((Qv/(well_area)-aquifer_df.R[row]/100)/aquifer_df.S[row]) # update T value due to drop in water table from depletion 
                     WL_array[row,t+1] = WL_array[row,t] + ((Qv/well_area)-aquifer_df.R[row]/100)/aquifer_df.S[row] # water depth below ground surface 
                     DTW_ts[(t*730):((t+1)*730)] = drawdown_total[1:]

                     max_drawdown = (-1*min(drawdown_total) - WL_array[row,t]) # max drawdown during current simulated year 
                     sat_thickness = D_array[row, t] - WL_array[row, t+1] # remaining saturated thickness 
                     entrance_V = (aquifer_df.Q[row]/(7.48*60))/(sat_thickness*3*2*math.pi*(8/2)/12*0.2) # casing entrance velocity, potential constrain, but not currently implemented 
                     if max_drawdown/sat_thickness > 0.8 and aquifer_df.m[row] > D_array[row,t]:
                         additional_depth_remaining = aquifer_df.m[row] - D_array[row, t]
                         if additional_depth_remaining > 30:
                             D_array[row, t+1] = D_array[row, t] + 30 # deepen well by 30 m 
                             T_array[row,t+1] = T_array[row,t+1] + 30*365*aquifer_df.K[row] # update Transmissivity 
                             sat_thickness = D_array[row, t+1] - WL_array[row, t+1]
                         else:
                             D_array[row, t+1] =  D_array[row, t] + additional_depth_remaining
                             T_array[row,t+1] = T_array[row,t+1] + additional_depth_remaining*365*aquifer_df.K[row] # update Transmissivity
                             sat_thickness = D_array[row, t+1] - WL_array[row, t+1]
                     else:
                         D_array[row, t+1] = D_array[row, t] # well depth remains the same
                     
                 else:
                     D_array[row,t+1] = 0
                     break
            
        
    #%% Calculate cost curve using pumping results 
    
    Ct = np.zeros((years+30)) # array of annual total costs
    Cc = np.zeros((years+30)) # array of annual capital costs
    Cm = np.zeros((years+30)) # array of annual maintenance costs
    Elec_cost = np.zeros((years)) # array of annual electricity costs for pumping
    
    if extended != 'reduced_capacity':
        #Ct = np.zeros((len(A_array[:,0]), years+30)) # array of annual total costs
        Cc = np.zeros((len(A_array[:,0]), years+30)) # array of annual capital costs
        Cm = np.zeros((len(A_array[:,0]), years+30)) # array of annual maintenance costs
        #Elec_cost = np.zeros((len(A_array[:,0]), years)) # array of annual electricity costs for pumping

    Annual_avg_PWL_depth = np.zeros((1, years))
    depletion_ind = np.zeros((1))    
    
    Annual_D = np.ones((D_array[0,:].size,1))
    for i in range(len(Annual_D)):
        Annual_D[i] = max(np.unique(D_array[:,i]))
    
    Q_annual = np.sum(Q_array, axis = 0)
    Days_annual = np.sum(Days_array, axis = 0)
    Area_annual = np.sum(A_array, axis = 0)
    
    # Make sure that there are valid outputs. If not, empty arrays are output
    if DTW_ts[0] < 0 and abs(DTW_ts[0]) < aquifer_df.m[row] and viable_first_Q > 0: 
        
        # Cost parameters
        #cost = 0.12 # cost per kwh in $
        eff = 0.75 # well efficiency
        rate = 0.06 # interest rate on loan 
        opr_P = 30 # meters of additional head for sprinkler system 
        n = 30 # lifetime of well 
        maintenance = 0.07 # annual maintence cost (% of initial capital cost)
        
        # Code block for reduced well yield and fixed well area
        if extended == 'reduced_capacity':
            
            
            if aquifer_df.m[0] - aquifer_df.WL[0] < 100:
                print('block_1')
                if Q_annual[0] < 200:
                    well_unit_cost = 75 #$/per foot
                else:
                    well_unit_cost = 75 #$/per foot
                    
                Ini_cost = well_unit_cost*3.28*Annual_D[0] # $/per foot, converted to $/meter
                
                for i in range(len(Ct[:])-30):
                    if i == 0:
                        CapCost = Ini_cost*(1+rate)**n*(rate/((1+rate)**n-1))
                        Cc[i] = CapCost  
                        Cm[i] = maintenance*Ini_cost # maintenance cost [% of initial cost]
                    
                    elif (i+1)/n - (i+1)//n == 0: # Replace well every n years (well lifetime), if reduced yeild, cheaper unit cost at 200 gpm and below
                        if Q_annual[i] < 200:
                            well_unit_cost = 75 #$/per foot
                        else:
                            well_unit_cost = 75 #$/per foot
                            
                        Ini_cost = well_unit_cost*3.28*Annual_D[i]
                        CapCost = Ini_cost*(1+rate)**n*(rate/((1+rate)**n-1)) # new cap cost per year for 30 yr period
                        Cc[i] = Cc[i]+CapCost
                        Cm[i] = maintenance*Ini_cost
                    
                    else: # D_array[row, i] > 0 and i > 0:
                        Cc[i] = Cc[i]+CapCost 
                        Cm[i] = maintenance*Ini_cost # maintenance cost [% of initial cost]
                        
                
            else:
                print('block_2')
                if Q_annual[0] < 200:
                    well_unit_cost = 75 #$/per foot
                else:
                    well_unit_cost = 75 #$/per foot
                    
                Ini_cost = well_unit_cost*3.28*D_array[row,0] # $75 per foot, converted to $/meter
    
                for i in range(len(Ct[:])-30):
                    if i == 0:
                        CapCost = Ini_cost*(1+rate)**n*(rate/((1+rate)**n-1))
                        Cc[i] = CapCost  
                        Cm[i] = maintenance*Ini_cost # maintenance cost [% of initial cost]
                        
                    elif (i+1)/n - (i+1)//n == 0: # Replace well every n years (well lifetime), if reduced yeild, cheaper unit cost at 200 gpm and below
                        if Q_annual[i] < 200:
                            well_unit_cost = 75 #$/per foot
                        else:
                            well_unit_cost = 75 #$/per foot
                            
                        Ini_cost = well_unit_cost*3.28*Annual_D[i]
                        CapCost = Ini_cost*(1+rate)**n*(rate/((1+rate)**n-1)) # new cap cost per year for 30 yr period
                        Cc[i] = Cc[i]+CapCost
                        Cm[i] = maintenance*Ini_cost
                        
                    elif Annual_D[i+1] - Annual_D[i] > 0: # add additional capital cost if depth is increased
                        Cc[i] = Cc[i-1]
                        Cc[i+1:i+31] += well_unit_cost*3.28*(Annual_D[i+1] - Annual_D[i])*(1+rate)**n*(rate/((1+rate)**n-1)) 
                        Ini_cost = well_unit_cost*3.28*Annual_D[i+1]
                        Cm[i] = maintenance*Ini_cost
                        
                    else: # D_array[row, i] > 0 and i > 0:
                        Cc[i] = Cc[i]+CapCost 
                        Cm[i] = maintenance*Ini_cost # maintenance cost [% of initial cost]
            
            # Calculate total annual cost (capital + pumping)
            print('cost_1')
            for i in range(years):
                Annual_avg_PWL_depth[0,i] = np.mean(-1*DTW_ts[aquifer_pump_indexes[0][:]+730*i]) # REFERENCES PUMP INDEXES FOR EACH SCENARIO
                if Annual_avg_PWL_depth[0,i] > aquifer_df.m[row] or Annual_avg_PWL_depth[0,i] == 0:
                    Ct[i] == 0
                else:
                    Elec_cost[i] = cost*((1000*9.81*(Annual_avg_PWL_depth[0,i]+opr_P)*Q_annual[i]/(264*60))/(1000*eff))*24*Days_annual[i]
                    Ct[i] = Elec_cost[i]+Cc[i]+Cm[i]
                
            if Annual_avg_PWL_depth[0,years-1] < aquifer_df.m[row] and Annual_avg_PWL_depth[0,years-1] > 0:
                depletion_ind = np.array([years])
                
            #elif max(Annual_avg_PWL_depth[0,:]) > aquifer_df.m[row]:
                #depletion_ind = min(np.argwhere(Annual_avg_PWL_depth[0,:] > aquifer_df.m[row]))
                
            else:
                depletion_ind = min(np.argwhere(Annual_avg_PWL_depth[0,:] == 0))-1
                
            unit_cost = np.zeros((1, years+30))
            unit_cost_per_acre = np.zeros((1, years+30))
            unit_cost_per_acre_ft = np.zeros((1, years+30))
            cost_per_one_eighth_cell = np.zeros((1, years+30))
            total_pumped_volume = np.zeros((1, years+30))
            Qv_annual = np.zeros(years+30)
            
            for i in range(depletion_ind[0]):
                Qv_annual[i] = Q_annual[i]*Days_annual[i]*1440*1/264.17 # annual pumped V in m^3 from one well
            
            unit_cost[0,:] = Ct/Qv_annual
            
            

            for j in range(int(depletion_ind[0])): # iterate for each pumping year until depletion 
                well_num = 247.1/Area_annual[0] # number of wells in 1 km2
                cost_per_one_eighth_cell[0,j] = Ct[j]*(192.52*well_num) # total number of wells in 1/8 deg cell (192.52 km x wells/km)
                
                if j == 0:
                    total_pumped_volume[0,j] = (Qv_annual[j] * (192.52*well_num) * (1/10**9)) 
                else:
                    total_pumped_volume[0,j] = (Qv_annual[j] * (192.52*well_num) * (1/10**9)) + total_pumped_volume[0,j-1]
                        
            QV_annual = Qv_annual * (192.52*well_num)
            Cost_annual = Ct * (192.52*well_num)
            Fixed_costs = (Cc + Cm) * (192.52*well_num)
            Elec_cost = Elec_cost * (192.52*well_num)
            Ct = Fixed_costs[0:years] + Elec_cost
            
        # Code block for capital costs for adding wells to offset reduced yield 
        else:
            
            Well_count_annual = np.zeros((1, years+30)) # number of wells per 1/8th deg grid cell 
            for i in range(years):
                if Area_annual[i] > 0:
                    Well_count_annual[0,i] = round((247.1/Area_annual[i]) * 192.52) # total number of wells in 1/8 deg cell (192.52 km x wells/km)
                else: 
                    continue
            Well_counts = np.unique(Well_count_annual)
            Well_counts = np.delete(Well_counts, 0)
            Added_well_count = np.zeros(len(Well_counts))
            for i in range(len(Well_counts)):
                if i == 0:
                    Added_well_count[i] = Well_counts[i]
                else:
                    Added_well_count[i] = Well_counts[i] - Well_counts[i-1]
                    
            Q_vals = np.unique(Q_annual)
            Q_vals = np.delete(Q_vals, 0)
            Start_indx = np.zeros(len(Q_vals))
            for i in range(len(Well_count_annual[0,:])-30):
                if i == 0:
                    counter = 1
                if Well_count_annual[0,i+1] - Well_count_annual[0,i] > 0:
                    Start_indx[counter] = int(i+1)
                    counter += 1 
            
            end_pumping = np.min(np.where(Q_annual == 0)) - 1 # last index of pumping 
            
            # Code block for Capital cost without well deepening (initial sat thickness < 100 m)       
            if aquifer_df.m[0] - aquifer_df.WL[0] < 100:
                print('block_3')
                if Q_annual[0] < 200:
                    well_unit_cost = 75 #$/per foot
                else:
                    well_unit_cost = 75 #$/per foot
                    
                Ini_cost = well_unit_cost*3.28*Annual_D[0] # $/per foot, converted to $/meter
                for j in range(len(Added_well_count)):
                    for i in range(len(Cm[0,:])-30):
                        offset = int(Start_indx[j])
                        if i == 0:
                            CapCost = Ini_cost*(1+rate)**n*(rate/((1+rate)**n-1))
                            Cc[j,i+offset] = CapCost * Added_well_count[j]
                            Cm[j,i+offset] = maintenance * Ini_cost * Added_well_count[j] # maintenance cost [% of initial cost]
                        
                        elif i + offset > end_pumping:
                            break
                        
                        elif (i+1)/n - (i+1)//n == 0: # Replace well every n years (well lifetime), if reduced yeild, cheaper unit cost at 200 gpm and below
                            if Q_annual[i+offset] < 200:
                                well_unit_cost = 75 # $/per foot
                            else:
                                well_unit_cost = 75 # $/per foot
                                
                            Ini_cost = well_unit_cost*3.28*Annual_D[i]
                            CapCost = Ini_cost*(1+rate)**n*(rate/((1+rate)**n-1)) # new cap cost per year for 30 yr period
                            Cc[j,i+offset] = Cc[j,i+offset]+ CapCost * Added_well_count[j]
                            Cm[j,i+offset] = maintenance * Ini_cost * Added_well_count[j]
                            
                        else: 
                            Cc[j,i+offset] = Cc[j,i+offset] + CapCost * Added_well_count[j]
                            Cm[j,i+offset] = maintenance * Ini_cost * Added_well_count[j] # maintenance cost [% of initial cost]
            
            # Code block for Capital cost with well deepening (initial sat thickness > 100 m)       
            else:
                print('block_4')
                if Q_annual[0] < 200:
                    well_unit_cost = 75 #$/per foot
                else:
                    well_unit_cost = 75 #$/per foot
                    
                Ini_cost = well_unit_cost*3.28*D_array[row,0] # $75 per foot, converted to $/meter
                
                for j in range(len(Added_well_count)):
                    for i in range(len(Cm[0,:])-30):
                        offset = int(Start_indx[j])
                        if i == 0:
                            CapCost = Ini_cost*(1+rate)**n*(rate/((1+rate)**n-1))
                            Cc[j,i+offset] = CapCost * Added_well_count[j]
                            Cm[j,i+offset] = maintenance * Ini_cost * Added_well_count[j] # maintenance cost [% of initial cost]
                        
                        elif i + offset > end_pumping:
                            break
                        
                        elif (i+1)/n - (i+1)//n == 0: # Replace well every n years (well lifetime), if reduced yeild, cheaper unit cost at 200 gpm and below
                            if Q_annual[i+offset] < 200:
                                well_unit_cost = 75 # $/per foot
                            else:
                                well_unit_cost = 75 # $/per foot
                                
                            Ini_cost = well_unit_cost*3.28*Annual_D[i]
                            CapCost = Ini_cost*(1+rate)**n*(rate/((1+rate)**n-1)) # new cap cost per year for 30 yr period
                            Cc[j,i+offset] = Cc[j,i+offset] + CapCost * Added_well_count[j]
                            Cm[j,i+offset] = maintenance * Ini_cost * Added_well_count[j]
                        
                        elif (Annual_D[i+1+offset] - Annual_D[i+offset]) > 0: # add additional capital cost if depth is increased
                            Cc[j,i] = Cc[j,i-1+offset]
                            Cc[j,i+1:i+31] += well_unit_cost*3.28*(Annual_D[i+1+offset] - Annual_D[i+offset])*(1+rate)**n*(rate/((1+rate)**n-1)) * Added_well_count[j]
                            Ini_cost = well_unit_cost*3.28*Annual_D[i+1+offset]
                            Cm[j,i+offset] = maintenance * Ini_cost * Added_well_count[j]
                                
                        else: 
                            #Cc[j,i+offset] = Cc[j,i+offset] + CapCost * Added_well_count[j]
                            Cc[j,i+offset] += CapCost * Added_well_count[j]
                            Cm[j,i+offset] = maintenance * Ini_cost * Added_well_count[j] # maintenance cost [% of initial cost]
        
                                
            # Calculate total annual cost (capital + pumping)
            Fixed_costs = np.sum(Cc, axis = 0)  + np.sum(Cm, axis = 0) # annual fixed costs 
            print('cost_2')
            for i in range(years):
                Annual_avg_PWL_depth[0,i] = np.mean(-1*DTW_ts[aquifer_pump_indexes[0][:]+730*i]) # REFERENCES PUMP INDEXES FOR EACH SCENARIO
                if Annual_avg_PWL_depth[0,i] > aquifer_df.m[row] or Annual_avg_PWL_depth[0,i] == 0:
                    Ct[i] == 0
                else:
                    Elec_cost[i] = cost*((1000*9.81*(Annual_avg_PWL_depth[0,i]+opr_P)*Q_annual[i]/(264*60))/(1000*eff))*24*Days_annual[i]*Well_count_annual[0,i]
                    Ct[i] = Elec_cost[i]+Fixed_costs[i]
                
            if Annual_avg_PWL_depth[0,years-1] < aquifer_df.m[row] and Annual_avg_PWL_depth[0,years-1] > 0:
                depletion_ind = np.array([years])
                
            #elif max(Annual_avg_PWL_depth[0,:]) > aquifer_df.m[row]:
                #depletion_ind = min(np.argwhere(Annual_avg_PWL_depth[0,:] > aquifer_df.m[row]))
                
            else:
                depletion_ind = min(np.argwhere(Annual_avg_PWL_depth[0,:] == 0))-1
                
            unit_cost = np.zeros((1, years+30))
            unit_cost_per_acre = np.zeros((1, years+30))
            unit_cost_per_acre_ft = np.zeros((1, years+30))
            cost_per_one_eighth_cell = np.zeros((1, years+30))
            total_pumped_volume = np.zeros((1, years+30))
            Qv_annual = np.zeros(years+30)
            
            for i in range(depletion_ind[0]):
                Qv_annual[i] = Q_annual[i] * Days_annual[i] * Well_count_annual[0,i] * 1440 *1/264.17 # annual pumped V in m^3
            
            unit_cost[0,:] = Ct/Qv_annual
                
            # Cumulative pumped volume in km^3 in a 1/8th grid cell 
            for j in range(int(depletion_ind[0])): # iterate for each pumping year until depletion 
                #well_num = 247.1/Area_annual[0] # number of wells in 1 km2
                #cost_per_one_eighth_cell[0,j] = Ct[0,j]*(192.52*well_num) # total number of wells in 1/8 deg cell (192.52 km x wells/km)
                if j == 0:
                    total_pumped_volume[0,j] = (Qv_annual[j] * (1/10**9)) 
                else:
                    total_pumped_volume[0,j] = (Qv_annual[j] * (1/10**9)) + total_pumped_volume[0,j-1]
            
            QV_annual = Qv_annual 
            Cost_annual = Ct
            #Ct = Ct/Well_count_annual[0,:]
            #Elec_cost = Elec_cost[0:years]/Well_count_annual[0,0:years]
            #Cc = np.sum(Cc, axis = 0)/Well_count_annual[0,:]
            #Cm = np.sum(Cm, axis = 0)/Well_count_annual[0,:]
    
    WL_ts = np.zeros(years+1)
    for i in range(len(WL_array[0,:])):
        WL_ts[i] = np.max(WL_array[:,i])
        
    Qv_ratio = Qv_annual/Qv_annual[0]
    for i in range(len(Qv_ratio)):
        if Qv_ratio[i] > 1:
            Qv_ratio[i] = 1
    
    for i in range(len(unit_cost[0,:])):
        if np.isnan(unit_cost[0,i]) == 1 or unit_cost[0,i] == np.inf:
            unit_cost[0,i] = 0 
    
    
    return unit_cost, total_pumped_volume, WL_ts, Qv_ratio, QV_annual, Cost_annual


