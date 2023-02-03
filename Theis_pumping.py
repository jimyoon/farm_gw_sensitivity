# -*- coding: utf-8 -*-

import numpy as np 
import pandas as pd 
import math 
import os 

# os.chdir('C:/Users/fere556/Documents/Pumping_cost')
W_lookup = pd.read_csv('./data_inputs/Theis_well_function_table.csv', header = "infer")
lookup_idx = pd.Index(W_lookup.W)

#%% Construct aquifer properties parameter space  
def Analytical(S, m, K, WL, R, IrrDepth, years):
    Q = [125, 250, 500, 750, 1000] # Well Q [gpm]
    
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
    days = [30, 40, 50, 60, 70, 80] # days of irrigation per year 
    area_arr = np.zeros(len(area)*len(days))
    
    for i in range(len(area)):
        area_arr[i*len(days):(i+1)*len(days)] = area[i]*np.ones(len(days))
        
    days_arr = np.zeros(len(area)*len(days))
    days_arr = np.tile(days, len(area))
        
    data = {'area': area_arr, 'days': days_arr, 'Q_125_Annual': np.zeros(len(area)*len(days)), 
            'Q_250_Annual': np.zeros(len(area)*len(days)), 'Q_500_Annual': np.zeros(len(area)*len(days)),
            'Q_750_Annual': np.zeros(len(area)*len(days)), 'Q_1000_Annual': np.zeros(len(area)*len(days)), 
            'Q_125_24h': np.zeros(len(area)*len(days)), 'Q_250_24h': np.zeros(len(area)*len(days)),
            'Q_500_24h': np.zeros(len(area)*len(days)), 'Q_750_24h': np.zeros(len(area)*len(days)),
            'Q_1000_24h': np.zeros(len(area)*len(days))}
    
    irr_depth_df = pd.DataFrame(data = data)
        
    for i in range(irr_depth_df.iloc[:,0].size):
        for j in range(len(Q)):
            irr_depth_df.iloc[i,j+2] = (Q[j]*1440*irr_depth_df.days[i]*1/264)/(irr_depth_df.area[i]*4046.86) # irrigation depth in meters, conversion of gpm to m^3 and area in acres to m^2
            irr_depth_df.iloc[i,j+2+len(Q)] = (Q[j]*1440*1/264)/(irr_depth_df.area[i]*4046.86)
            
#%% For each aquifer scenario, determine (1) area irrigated AND (2) days of irrigation 
# for defined annual Irrigation Depth (IrrDepth)
     
    Irr_depth_m = .0256*IrrDepth # irrigation depth in m, provided irr depth [in] 
    irr_depth_df_Q_ilocs = np.array([[125, 250, 500, 750, 1000],[2,3,4,5,6]])
    
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
                aquifer_df.iloc[i,11] = temp[temp_sorted[-1],1] # days of irrigation for area served 
    
      
#%% One year if simulated aquifer response at half day resolution (time steps = 0.5 d or 1/730 yr)
    T_array = np.zeros((aquifer_df.iloc[:,0].size, 51))
    T_array[:,0] = aquifer_df.Trans
    One_yr_ts = np.zeros((aquifer_df.iloc[:,0].size, 730*1))
    WL_array = np.zeros((aquifer_df.iloc[:,0].size, 51))
    WL_array[:,0] = aquifer_df.WL
    
    for x in range(aquifer_df.iloc[:,0].size): # of aquifer scenarios 
        for t in range(1): # years of simulated pumping
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
                WL_array[x,t+1] = WL_array[x,t] + ((Qv/well_area)-aquifer_df.R[x]/100)/aquifer_df.S[row]
                One_yr_ts[x,(t*730):((t+1)*730)] = drawdown_total[1:] 
            else:
                continue
    
#%% Calculate max drawdown over 1 yr pumping period
    
    for i in range(aquifer_df.iloc[:,0].size):
        if aquifer_df.Trans[i] == 0:
            aquifer_df.s[i] = -999
        else:
            aquifer_df.PWL[i] = -1*min(One_yr_ts[i,:]) # max depth to water during 1 year of cyclical pumping (m below surface)
            aquifer_df.s[i] = -1*min(One_yr_ts[i,:])-aquifer_df.WL[i] # max drawdown curing 1 year of cyclical of pumping (m)
    
    # filter out non-viable scenarios
    aquifer_df['viable'] = np.ones(aquifer_df.iloc[:,0].size)
    for i in range(aquifer_df.iloc[:,0].size):
        if aquifer_df.PWL[i] > 0.66 * aquifer_df.m[i] or aquifer_df.Trans[i] == 0 or aquifer_df.s[i] > 75 or aquifer_df.Area[i] == 0:
            aquifer_df.viable[i] = 0
   
#%% Choose preferrable viable aquifer scenario based on Area and Days

    viable_scenarios = aquifer_df.copy()
    viable_scenarios.where(aquifer_df.viable > 0, inplace = True)
    viable_scenarios = viable_scenarios.dropna()
    
    if np.sum(aquifer_df.viable[:].values) == 0: # exit if no viable scenarios 
        print('Non-viable')

    else:
        Area_indx = viable_scenarios.Area.sort_values() 
        Area_indx = np.where(Area_indx == np.max(Area_indx))
        Days_count = viable_scenarios.Days[Area_indx[0][:]].sort_values()
        row = int(Days_count.index[0])
        
        
    #%% Half day resolution (time steps = 0.5 d or 1/730 yr)
    
    T_array = np.zeros((aquifer_df.iloc[:,0].size, years+1))
    T_array[:,0] = aquifer_df.Trans
    DTW_ts = np.zeros((aquifer_df.iloc[:,0].size, 730*years))
    WL_array = np.zeros((aquifer_df.iloc[:,0].size, years+1))
    WL_array[:,0] = aquifer_df.WL
    aquifer_pump_indexes = [] # list of lists where each aquifer scenario has a record of pumping indexes for PWL calcs
    
    for t in range(years): # years of simulated pumping
        if T_array[row,t] > 0: 
            Q_gpm = aquifer_df.Q[row]
            Q = Q_gpm*1440*1/264*365 # pumping rate in GPM converted, GPD (mult by 1440), to m^3/day (divide by 264), to m^3/year
            Qv = Q_gpm*1440*1/264*aquifer_df.Days[row] # pumped volume (m^3) per yr 
            r = 0.5 # radial distance of drawdown estimate
            well_area = aquifer_df.Area[row]*4046.86 # AREA SERVED BASED ON WELL CAPACITY (convert from acres to m^2)
        
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
            drawdown_total = -1*np.sum(superposition_arr, axis = 0) - WL_array[row,t]
            T_array[row,t+1] = T_array[row,t] - aquifer_df.K[row]*365*((Qv/(well_area))/aquifer_df.S[row]) # update T value due to drop in water table from depletion 
            WL_array[row,t+1] = WL_array[row,t] + ((Qv/well_area)-aquifer_df.R[row]/100)/aquifer_df.S[row]
            DTW_ts[row,(t*730):((t+1)*730)] = drawdown_total[1:] 
        else:
            continue
    
    # np.savetxt('fifty_year_ts.csv', Fifty_yr_ts, delimiter = ',')
    # aquifer_df.to_csv('aquifer_df.csv', header= True)
    
#%% Calculate cost curves for each scenario 
    
    Ct = np.zeros((1, years))
    Elec_cost = np.zeros((1, years))
    Annual_avg_PWL_depth = np.zeros((1, years))
    depletion_ind = np.zeros((1))
    aquifer_df['fixed_cost'] = np.zeros(aquifer_df.iloc[:,0].size) 
    unit_cost = np.zeros((1, years))
    unit_cost_per_acre = np.zeros((1, years))
    cost_per_one_eight_cell = np.zeros((1, years))
    total_pumped_volume = np.zeros((1, years))
    

    if DTW_ts[row,0] < 0 and abs(DTW_ts[row,0]) < aquifer_df.m[row] and aquifer_df.viable[row] > 0: 
        cost = 0.12 # cost per kwh in $
        eff = 0.75 # well efficiency 
        Ini_cost = 75*3.28*aquifer_df.m[row] # $75 per foot, converted to $/meter
        rate = 0.06 # interest rate on loan 
        opr_P = 30 # meters of additional head for sprinkler system 
        n = 50 # years of pumping 
        CapCost = Ini_cost*(1+rate)**n*(rate/((1+rate)**n-1))
        Cc = np.concatenate(((np.tile(CapCost,n)),np.zeros(50-n)))
        Cm = 0.07*Ini_cost
        Qv = Q_gpm*1440*1/264*aquifer_df.Days[row] # pumped volume (m^3) per yr 
        aquifer_df.fixed_cost[row] = (Cc[0]+Cm)/(Qv) # fixed cost in terms of $/m^3 pumped per year 
        
        for i in range(years):
            Annual_avg_PWL_depth[0,i] = np.mean(-1*DTW_ts[row, aquifer_pump_indexes[0][:]+730*i]) # REFERENCES PUMP INDEXES FOR EACH SCENARIO
            if Annual_avg_PWL_depth[0,i] > aquifer_df.m[row] or Annual_avg_PWL_depth[0,i] == 0:
                Ct[0,i] == 0
            else:
                Elec_cost[0,i] = cost*((1000*9.81*(Annual_avg_PWL_depth[0,i]+opr_P)*aquifer_df.Q[row]/(264*60))/(1000*eff))*24*aquifer_df.Days[row]
                Ct[0,i] = Elec_cost[0,i]+Cc[i]+Cm
            
        if Annual_avg_PWL_depth[0,years-1] < aquifer_df.m[row] and Annual_avg_PWL_depth[0,years-1] > 0:
            depletion_ind = np.array([years])
            
        elif max(Annual_avg_PWL_depth[0,:]) > aquifer_df.m[row]:
            depletion_ind = min(np.argwhere(Annual_avg_PWL_depth[0,:] > aquifer_df.m[row]))
            
        else:
            depletion_ind = min(np.argwhere(Annual_avg_PWL_depth[0,:] == 0))
            
        # COMMENTED OUT PLOTTING SCRIPT 
        # fig, axs = plt.subplots(nrows = 2, ncols = 1, figsize = (8,6))
        # axs[0].plot(np.arange(depletion_indx[0]*730)/730,DTW_ts[x,0:depletion_indx[0]*730], label = 'Q = ' + str(Q_gpm) +' GPM for 52 days/yr')
        # axs[0].plot([0,50*729/730],[-1*aquifer_df.m[x],-1*aquifer_df.m[x]], label = 'Aquifer thickness') 
        # axs[0].legend(fontsize = 'small')
        # axs[0].set_ylabel('Depth to water [m]')
        # axs[0].set_title('K =' + str(aquifer_df.K[x])+ 'm/d , m = ' + str(aquifer_df.m[x]) + ' meters')
        # axs[1].plot(np.arange(depletion_indx[0])+1,(1/(Qv))*Ct[x,0:depletion_indx[0]], label = 'Total cost', color = 'k')
        # axs[1].plot(np.arange(depletion_indx[0])+1,(1/(Qv))*Elec_cost[x,0:depletion_indx[0]], label = 'Electricity cost', color = 'k', linestyle = '--')
        # axs[1].set_ylabel('Cost ($/m^3)')
        # axs[1].set_ylim(0,max((1/(Qv))*Ct[x,0:depletion_indx[0]]))
        # axs[1].set_xlim(0,50)
        # axs[1].set_xlabel('Years')
        # axs[1].set_title('Unit cost ($/m^3) over 50 years')
        # axs[1].legend(fontsize = 'small')
        # plt.tight_layout(pad = 1.2)
        # plt.savefig('50_year_ts_row' + str(x) + '.png')
    
        # fig, axs = plt.subplots(nrows = 1, ncols = 1, figsize = (8,4))
        # cost_curve[x,0:depletion_indx[0]] = (1/(Qv))*Ct[x,0:depletion_indx[0]]
        # pumped_V = np.array([Qv*(np.arange(depletion_indx[0])+1)*num_well*(10**-9)])
        # axs.plot((1/(Qv))*Ct[x,0:depletion_indx[0]],Qv*(np.arange(depletion_indx[0])+1)*num_well*(10**-9), color = 'k')
        # axs.set_ylabel('Volume extracted (km^3)')
        # axs.set_xlabel('Cost ($/m^3)')
        # axs.set_xlim(0,max((1/(Qv))*Ct[x,0:depletion_indx[0]]))
        # axs.set_title('K =' + str(aquifer_df.K[x])+ 'm/d , m = ' + str(aquifer_df.m[x]) + ' meters, '+ str(Q_gpm) +' GPM for 52 days/yr')
        # plt.savefig('cost_curve_row' + str(x) + '.png')
    
        # irrigation_depth[x] = aquifer_df.Q[row]/(50000**2/aquifer_df.num_wells[row])
        # fig, axs = plt.subplots(nrows = 1, ncols = 1, figsize = (8,4))
        # axs.scatter(1,Qv/(750**2), color = 'k', label = 'Q = ' + str(Q_gpm) +' GPM for 52 days/yr')
        # axs.scatter(1, aquifer_df.Q[row]/(50000**2/aquifer_df.num_wells[row]), color = 'b', 
        #             label = 'Superwell Q =' + str(np.around(aquifer_df.Q[row]*264/(1440*365))) \
        #                 + ' GPM for 365 days/yr' )
        # axs.set_ylabel('Annual irrigation depth (m)')
        # axs.legend()
        # plt.savefig('annual_irrigation_depth_row' + str(x) + '.png')
    
    # fig, axs = plt.subplots(nrows = 1, ncols = 1, figsize = (2,3))
    # plt.scatter(np.ones(12),irrigation_depth, color = 'b')
    # plt.scatter(1, Qv/(750**2), color = 'k')
    # plt.xticks([1],[''])
    # plt.xlim(0,2)
    # plt.ylabel('Annual irrigation depth (m)')
    # plt.ylim(0, Qv/(750**2)+.01)
    
    # Unit cost per m^3 
    aquifer_df['Qv'] = np.zeros(aquifer_df.m.size)
    for i in range(aquifer_df.m.size):
        aquifer_df.Qv[i] = aquifer_df.Q[i]*aquifer_df.Days[i]*1440*1/264.17 # annual pumped V in m^3
    
    #for i in range(aquifer_df.m.size):
    unit_cost[0,:] = Ct/aquifer_df.Qv[row]
        
    # Unit cost per acre irrigated 
    #for i in range(aquifer_df.m.size):
    unit_cost_per_acre[:] = Ct/aquifer_df.Area[row]
    
        
    # Cost per 1/8 deg cell, volume per 1/8 cell
    #for i in range(aquifer_df.m.size):
    for j in range(int(depletion_ind[0])): # iterate for each pumping year until depletion 
        well_num = 247.1/aquifer_df.Area[row] # number of wells in 1 km2
        cost_per_one_eight_cell[0,j] = Ct[0,j]*(192.52*well_num) # total number of wells in 1/8 deg cell (192.52 km x wells/km)
        if j == 0:
            total_pumped_volume[0,j] = (aquifer_df.Qv[row] * (192.52*well_num) * (1/10**9)) 
        else:
            total_pumped_volume[0,j] = (aquifer_df.Qv[row] * (192.52*well_num) * (1/10**9)) + total_pumped_volume[0,j-1]
    
    return unit_cost, total_pumped_volume



