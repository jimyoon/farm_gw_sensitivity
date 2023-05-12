# -*- coding: utf-8 -*-
"""
Created on Thu May  4 16:55:40 2023

@author: fere556
"""

import numpy as np 
import pandas as pd 
import math 
import os 
import matplotlib.pyplot as plt

os.chdir('C:/Users/fere556/Desktop/Superwell')

# import Superwell_python_with_deepening_for_ABM_on_the_fly.cost_curve as cost_curve

import Superwell_for_ABM_on_the_fly

#S=0.25, m=110, WL=10, R=10, IrrDepth=12

years = 200
IRR_DEPTH = [0.25, 0.5, 0.75]
unit_cost = np.zeros((3, years))
cum_vol = np.zeros((3, years))
cost_per_well = np.zeros((3,years))
water_depth = np.zeros((3,years))
total_head = np.zeros((3,years))
well_area = np.zeros((3, years))
capacity_reduction = np.zeros((3, years))
well_Q = np.zeros((3, years))

for i in range(3):
    output = Superwell_for_ABM_on_the_fly.cost_curve(S = 0.25, m = 120, K = 10, WL = 10, R = 0.0, IRR_DEPTH = IRR_DEPTH[i], NUM_YEARS = 200, ELECTRICITY_RATE = 0.12)
    unit_cost[i,:] = output[8]
    cum_vol[i,:] = output[1]
    cost_per_well[i,:] = output[5]
    water_depth[i,:] = output[2]
    total_head[i,:] = output[9]    
    well_area[i,:] = output[6] * 0.000247105
    capacity_reduction[i,:] = output[3]
    well_Q[i,:] = output[10]
    
plt.figure()
for i in range(3):
   plt.plot(cum_vol[i,:], unit_cost[i,:], label = str(IRR_DEPTH[i])) 
plt.ylabel('$/acrefoot')
plt.xlabel('cumulative pumped volume (km^3)')
plt.legend(title = 'Irrigation Depth')

plt.figure()
for i in range(3):
   plt.plot(cost_per_well[i,:], label = str(IRR_DEPTH[i])) 
plt.ylabel('cost per well ($/yr)')
plt.xlabel('Year')
plt.legend(title = 'Irrigation Depth')

plt.figure()
for i in range(3):
   plt.plot(water_depth[i,:], label = str(IRR_DEPTH[i])) 
plt.ylabel('Water Depth')
plt.xlabel('Year')
plt.legend(title = 'Irrigation Depth')

plt.figure()
for i in range(3):
   plt.plot(total_head[i,:], label = str(IRR_DEPTH[i])) 
plt.ylabel('Total head')
plt.xlabel('Year')
plt.legend(title = 'Irrigation Depth')

plt.figure()
for i in range(3):
   plt.plot(cum_vol[i,:],total_head[i,:], label = str(IRR_DEPTH[i])) 
plt.ylabel('Total head')
plt.xlabel('cumulative pumped volume (km^3)')
plt.legend(title = 'Irrigation Depth')

plt.figure()
for i in range(3):
   plt.plot(well_area[i,:], label = str(IRR_DEPTH[i])) 
plt.ylabel('Well area')
plt.xlabel('Year')
plt.legend(title = 'Irrigation Depth')

plt.figure()
for i in range(3):
    end_indx = np.max(cum_vol)
    plt.plot(capacity_reduction[i,:], label = str(IRR_DEPTH[i]))
plt.ylabel('Capacity fraction')
plt.xlabel('Year')
plt.legend(title = 'Irrigation Depth')

plt.figure()
for i in range(3):
    #end_indx = np.max(cum_vol)
    end_indx = np.where(cum_vol[i,:] == np.max(cum_vol[i,:])) 
    plt.plot(cum_vol[i,0:int(end_indx[0])],capacity_reduction[i,0:int(end_indx[0])], label = str(IRR_DEPTH[i])) 
plt.ylabel('Capacity fraction')
plt.xlabel('cumulative pumped volume (km^3)')
plt.legend(title = 'Irrigation Depth')

plt.figure()
for i in range(3):
    end_indx = np.max(cum_vol)
    plt.plot(cum_vol[i,:], well_Q[i,:] * 264 * 60, label = str(IRR_DEPTH[i])) 
plt.ylabel('Well Q (gpm)')
plt.xlabel('cumulative pumped volume (km^3)')
plt.legend(title = 'Irrigation Depth')






plt.figure()
plt.plot(output[1], output[8])
plt.ylabel('$/acrefoot')
plt.xlabel('cumulative pumped volume (km^3)')

plt.figure()
plt.plot(output[6] * 0.000247105)
plt.ylabel('Well area ')
plt.xlabel('Year')

plt.figure()
plt.plot(-1*output[2])
plt.ylabel('Water Table Depth')
plt.xlabel('year')

plt.figure()
plt.plot(output[3])
plt.ylabel('Pumping rate decline')
plt.xlabel('year')

plt.figure()
plt.plot(output[8])
plt.ylabel('$/acrefoot')
plt.xlabel('year')