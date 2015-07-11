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

# 1D FDTD CPML
# Based on Understanding the Finite-Difference Time-Domain Method, John
# B. Schneider, www.eecs.wsu.edu/~schneidj/ufdtd, 2010.
# and Numerical electromagnetics : the FDTD method / Umran S. Inan,
# Robert A. Marshall. 2011

import numpy as np
import matplotlib.pyplot as plt
from time import sleep

imp0=377.0  # Free space impedance = sqrt(u0/e0)

size=300  # Domain size
#Source 
source_width = int(size/20.0)
delay = 10*source_width

source_x = int(size/2.0)
def source(current_time, delay, source_width):
    return np.exp(-(current_time-delay)**2/(2.0 * source_width**2))

#Model
total_steps = int(size+delay)  # Time stepping
frame_interval = int(total_steps/12.0)
all_steps = np.linspace(0, size-1, size)

ey = np.zeros(size)
hz = np.zeros(size)


#CPML (Inan pp.228-230)
dx = 1.0
R0 = 1e-6
m = 4  # Order of polynomial grading
pml_width = 10.0
sxmax = -(m+1)*np.log(R0)/2/imp0/(pml_width*dx)
sx = np.zeros(size)
sxm = np.zeros(size)
Phx = np.zeros(size)
Pex = np.zeros(size)

for mm in xrange(int(pml_width)):
    sx[mm] = sxmax*((pml_width-mm)/pml_width)**m
    sxm[mm] = sxmax*((pml_width-mm)/pml_width)**m  # Shifted to the right
    sx[size-mm-1] = sxmax*((pml_width-mm-0.5)/pml_width)**m  
    sxm[size-mm-1] = sxmax*((pml_width-mm)/pml_width)**m
aex = np.exp(-sx*imp0)-1
bex = np.exp(-sx*imp0)
ahx = np.exp(-sxm*imp0)-1
bhx = np.exp(-sxm*imp0)

x = np.arange(1,size-1,1)

for time in xrange(total_steps+1):
    Phx[x] = bhx[x]*Phx[x] + ahx[x]*(ey[x+1] - ey[x])
    hz[x] = hz[x] + (ey[x+1] - ey[x])/imp0 + Phx[x]/imp0
    Pex[x+1] = bex[x+1]*Pex[x+1] + aex[x+1]*(hz[x+1]-hz[x])
    ey[x+1] = ey[x+1] + (hz[x+1]-hz[x])*imp0 +Pex[x+1]*imp0
    ey[source_x] += source(time, delay, source_width)
    if time % frame_interval == 0:
    # tshow = size/2- source_width +10
    # if time > tshow and time < tshow+20:
        fig, axs = plt.subplots(2,1)#, sharey=True, sharex=True)
        fig.tight_layout()
        plt.title("After t=%i"%time)
        axs[0].plot(all_steps, ey, label="Ey")
        axs[0].plot(all_steps, hz*imp0,  label="Hz*imp0")
        axs[0].legend(loc='upper left')
        crop = np.arange(-5,0,1)
        axs[1].plot(all_steps[crop], ey[crop],all_steps[crop], hz[crop]*imp0,  label="Ey")
        plt.savefig("step0-at-time-%i.png"%time,pad_inches=0.02, bbox_inches='tight')
        plt.draw()
        #    plt.show()
        plt.clf()
        plt.close()





