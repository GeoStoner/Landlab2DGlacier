# -*- coding: utf-8 -*-
"""
Created on Sun Mar 13 10:09:31 2016

@ Created by Ryan Stoner
@ expanded from: gtucker
"""

"""
Initialize
"""

from landlab import RasterModelGrid
import numpy as np
import pylab as plt
from landlab.plot.imshow import imshow_node_grid
plt.clf() 

# Definte parameters
gamma = 0.01        # mass balance coefficient, 1/yr
ELA = 2500.0        # equilibrium-line altitude, m
max_elev = 2550.0   # maximum elevation of valley floor, m
dx = 500.0          # node spacing, m
dt = 0.1            # time step, yr
valley_slope = 0.1
side_slope = 0.2
num_rows = 50
num_cols = 50 
tmax = 100          # yr
t = np.arange(0,tmax,dt)
nsteps = int(len(t)/dt)

tplot = 100         # time between plots
nplots = tmax/tplot # number of plots

usl = 0.1          # sliding velocity, m/yr
A = 2.1*10**-16     # yr^-1*Pa^-3
icedens = 850       # kg/m^3
g =9.81             # m/s^2

# Create grid
mg = RasterModelGrid(num_rows, num_cols, dx)

# Create data fields
z = mg.add_empty('node', 'elevation')
H = mg.add_zeros('node', 'ice_thickness')
z_ice = mg.add_empty('node', 'ice_elevation')
dZicedx = mg.add_zeros('link', 'ice_surface_slope')

# Initialize elevation ("open book" shape)
# (Note: using "z[:] =" instead of "z =" means you are using the same block of
# memory instead of making a new copy)
z[:] = max_elev - mg.node_x * valley_slope
z +=  side_slope * np.abs(mg.node_y - mg.dx * ((num_rows - 1) / 2.0))
z_ice[:] = z + H

mg.set_closed_boundaries_at_grid_edges(True, True, True, True)
core_nodes = mg.core_nodes


pltime = 0          # Initial time, yr 

"""
Loop
"""
for i in range(nsteps):

    # Specific accumulation / melt
    b = gamma * (z_ice - ELA)
    
    # Ice-surface slope
    dZicedx[mg.active_links] = mg.calculate_gradients_at_active_links(z_ice)
    
    # Thickness at links
    H_edge = mg.map_mean_of_link_nodes_to_link('ice_thickness')
    
    # Calculate ice flux at links
    Q = usl * H_edge + A * (icedens*g)*(np.abs(dZicedx)**3)*(H_edge**5)/5;
    # Calculate dHdt
    dHdt = b - mg.calculate_flux_divergence_at_nodes(Q[mg.active_links])     
    
   
    H += np.maximum(dHdt* dt, 0)

    # Update ice-surface elevation
    z_ice[:] = z + H
    
    
    pltime += dt
    if (pltime%tplot)<=0.1:
        plt.clf() 
        plt.figure(1)
        
        plt.ion()
        gp1 = plt.subplot(211)
        
        imshow_node_grid(mg,'ice_thickness')
        gp2 = plt.subplot(212)
        
        imshow_node_grid(mg,'elevation')
        plt.pause(0.0001)
"""
Finalizing
"""
plt.subplots_adjust(hspace=0.4)     # Make sure x-labels do not overlap
gp1 = plt.ylabel('N-S distance (km)')
gp1 = plt.xlabel('E-W distance (km)')
gp1 = plt.title('Thickness of ice over total area')
gp2 = plt.ylabel('N-S distance (km)')
gp2 = plt.xlabel('E-W distance (km)')
gp2 = plt.title('Initial topography')