#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#
#    Copyright (C) 2015  Konstantin Ladutenko <kostyfisik@gmail.com>
#
#    This file is part of python-scattnlay
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

# 1D FDTD vanish E field
# Based on Understanding the Finite-Difference Time-Domain Method, John
# B. Schneider, www.eecs.wsu.edu/~schneidj/ufdtd, 2010.

import numpy as np
import math as m
import matplotlib.pyplot as plt
from time import sleep

imp0=377.0  # Free space impedance

size=1800  # Domain size
#Dielectric distribution
epsilon = 5
eps= np.ones(size)
eps[:] = epsilon

# setting  ABC constants _AFTER_ epsilon (we need speed of ligth in media)
# Taflove, eq. 6.35
c = 1/np.sqrt(epsilon)
a = (c-1)/(c+1)
b = 2/(c + 1)
# Left boundary
wl_nm1,wl_n,wl_np1 = 0,0,0 # Field at x=0 at time steps n-1, n, n+1
wlp1_nm1,wlp1_n,wlp1_np1 = 0,0,0 # Field at x=1 at time steps n-1, n, n+1
# Right boundary
wr_nm1,wr_n,wr_np1 = 0,0,0 # Field at x=size at time steps n-1, n, n+1
wrm1_nm1,wrm1_n,wrm1_np1 = 0,0,0 # Field at x=size-1 at time steps n-1, n, n+1

#Source 
source_width = 30.0*np.sqrt(epsilon)
delay = 10*source_width
source_x = int(1.0*size/2.0)  #Source position
def source(current_time, delay, source_width):
    return m.exp(-(current_time-delay)**2/(2.0 * source_width**2))


#Model
total_steps = int(size*3.0+delay)  # Time stepping
frame_interval = int(total_steps/15.0)
all_steps = np.linspace(0, size-1, size)

#Inital field E_z and H_y is equal to zero
ez = np.zeros(size)
hy = np.zeros(size)
x = np.arange(0,size-1,1)
#print(x)
for time in xrange(total_steps):
    ######################
    #Magnetic field
    ######################
    hy[x] = hy[x] + (ez[x+1] - ez[x])/imp0
    #Evaluate Mur ABC value (eq. 6.35 Taflove)
    wrm1_np1 = hy[-2]
    wr_np1 = -wrm1_nm1 + a*(wrm1_np1+wr_nm1) + b*(wr_n+wrm1_n)
    hy[-1] = wr_np1
    #Cycle field values at boundary
    wr_nm1, wrm1_nm1 = wr_n, wrm1_n
    wr_n, wrm1_n = wr_np1, wrm1_np1
    ######################
    #Electric field
    ######################
    ez[x+1] = ez[x+1] + (hy[x+1]-hy[x])*imp0/eps[x+1]
    ez[source_x] += source(time, delay, source_width)
    #Evaluate Mur ABC value (eq. 6.35 Taflove)
    wlp1_np1 = ez[1]
    wl_np1 = -wlp1_nm1 + a*(wlp1_np1+wl_nm1) + b*(wl_n+wlp1_n)
    ez[0] = wl_np1
    #Cycle field values at boundary
    wl_nm1, wlp1_nm1 = wl_n, wlp1_n
    wl_n, wlp1_n = wl_np1, wlp1_np1
    ######################
    # Output
    ######################
    if time % frame_interval == 0:
        plt.clf()
        plt.title("Ez after t=%i"%time)
        plt.plot(all_steps, ez, all_steps, hy*imp0)
        plt.show()




