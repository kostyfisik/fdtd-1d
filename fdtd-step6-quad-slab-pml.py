#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#
#    Copyright (C) 2015  Konstantin Ladutenko <kostyfisik@gmail.com>
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

# 1D FDTD half wavelength slab
# Based on Understanding the Finite-Difference Time-Domain Method, John
# B. Schneider, www.eecs.wsu.edu/~schneidj/ufdtd, 2010.

import numpy as np
import matplotlib.pyplot as plt
from time import sleep
from matplotlib import patches

imp0=377.0  # Free space impedance

size=400  # Domain size
wavelength = int(size/5.0) #in host media
factor = 4.0 # for slab lambda/4 

#Dielectric distribution
epsilon1 = 1.2 # host
epsilon2 = 4 # slab
n1 = np.sqrt(epsilon1)
n2 = np.sqrt(epsilon2)
print("Epsilon in slab was set to be %f"%epsilon2)
n2 = 1.0/int(wavelength*n1/n2/factor)*wavelength*n1/factor
epsilon2=n2**2
print("Corrected epsilon in slab is %f (to fit exactly to space step)"%epsilon2)
n3 = (n2**2)/n1
epsilon3=n3**2
eps= np.ones(size)
eps[:] = epsilon1
eps[int(size/2.0):int(size/2.0)+int(wavelength*n1/n2/factor)] = epsilon2
eps[int(size/2.0)+int(wavelength*n1/n2/factor):] = epsilon3
print("After slab epsilon is %f"%epsilon3)

refeps= np.ones(size)
refeps[:] = epsilon1

c = 1/np.sqrt(eps[0])
cref = c
#Source 
source_width = wavelength*3.0*np.sqrt(max(epsilon1,epsilon2))
delay = 10*source_width
source_x = int(1.0*size/10.0)  #Source position
def source(current_time, delay, source_width):
    amp =  np.exp(-(current_time-delay)**2/(2.0 * source_width**2))
    if current_time > delay:
        amp = 1.0
    return amp/np.sqrt(epsilon1)*np.sin(2*np.pi*current_time*cref/wavelength)

#Model
total_steps = int(size*4.0+delay)  # Time stepping
frame_interval = int(total_steps/35.0)
all_steps = np.linspace(0, size-1, size)

#CPML (Inan pp.228-230)
dx = 1.0
R0 = 1e-5
m = 4.0  # Order of polynomial grading
pml_width = 20.0
sxmax = -(m+1)*np.log(R0)/2/imp0/(pml_width*dx)
sx = np.zeros(size)
sxm = np.zeros(size)

for mm in xrange(int(pml_width)):
    sx[mm+1] = sxmax*((pml_width-mm-0.5)/pml_width)**m
    sxm[mm] = sxmax*((pml_width-mm)/pml_width)**m  # Shifted to the right
    sx[size-mm-1] = sxmax*((pml_width-mm-0.5)/pml_width)**m  
    sxm[size-mm-1] = sxmax*((pml_width-mm)/pml_width)**m
# print(sx, sxm)
aex = np.exp(-sx*imp0)-1
bex = np.exp(-sx*imp0)
ahx = np.exp(-sxm*imp0)-1
bhx = np.exp(-sxm*imp0)

x = np.arange(0,size-1,1)
x1 = np.arange(0,size-1,1)

#Inital field E_z and H_y is equal to zero
ez = np.zeros(size)
hy = np.zeros(size)
refez = np.zeros(size)
refhy = np.zeros(size)
Phx = np.zeros(size)
Pex = np.zeros(size)
refPhx = np.zeros(size)
refPex = np.zeros(size)

for time in xrange(total_steps):
    ######################
    #Magnetic field
    ######################
    Phx[x] = bhx[x]*Phx[x] + ahx[x]*(ez[x+1] - ez[x])
    hy[x] = hy[x] + (ez[x+1] - ez[x])/imp0 + Phx[x]/imp0
    ######################
    #Magnetic field reference
    ######################
    refPhx[x] = bhx[x]*refPhx[x] + ahx[x]*(refez[x+1] - refez[x])
    refhy[x] = refhy[x] + (refez[x+1] - refez[x])/imp0 + refPhx[x]/imp0
    ######################
    #Electric field
    ######################
    Pex[x1+1] = bex[x1+1]*Pex[x1+1] + aex[x1+1]*(hy[x1+1]-hy[x1])
    ez[x1+1] = ez[x1+1] + (hy[x1+1]-hy[x1])*imp0/eps[x1+1] +Pex[x1+1]*imp0/eps[x1+1]
    ez[source_x] += source(time, delay, source_width)
    ######################
    #Electric field reference
    ######################
    refPex[x1+1] = bex[x1+1]*refPex[x1+1] + aex[x1+1]*(refhy[x1+1]-refhy[x1])
    refez[x1+1] = refez[x1+1] + (refhy[x1+1]-refhy[x1])*imp0/refeps[x1+1] +refPex[x1+1]*imp0/refeps[x1+1]
    refez[source_x] += source(time, delay, source_width)

    ######################
    # Output
    ######################
    if time % frame_interval == 0:
    #if time  == 6825:
        #plt.clf()
        fig, axs = plt.subplots(3,1)#, sharey=True, sharex=True)
        fig.tight_layout()
        #axs[0].title("Ez after t=%i"%time)
        axs[0].plot(all_steps, ez,all_steps, refez)        
        s1 = patches.Rectangle((int(size/2.0), -100), int(wavelength*n1/n2/factor), 200.0, zorder=0,
                               color='blue',alpha = 0.2)
        axs[0].add_patch(s1)

        #axs[1].title("Diff after t=%i"%time)
        axs[1].plot(all_steps, ez-refez)        
        s2 = patches.Rectangle((int(size/2.0), -100), int(wavelength*n1/n2/factor), 200.0, zorder=0,
                               color='blue',alpha = 0.2)
        axs[1].add_patch(s2)
        axs[2].plot(all_steps[:int(size/2.0)], ez[:int(size/2.0)]-refez[:int(size/2.0)])        
        s3 = patches.Rectangle((int(size/2.0), -100), int(wavelength*n1/n2/factor), 200.0, zorder=0,
                               color='blue',alpha = 0.2)
        axs[2].add_patch(s3)
        plt.savefig("step6-at-time-%i-pml.png"%time,pad_inches=0.02, bbox_inches='tight')
        plt.draw()
        #    plt.show()
        plt.clf()
        plt.close()
        




