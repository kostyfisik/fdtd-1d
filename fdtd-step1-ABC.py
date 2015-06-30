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

size=500  # Domain size
#Source 
source_width = 30.0
delay = 5*source_width


source_x = int(size/2.0)
def source(current_time, delay, source_width):
    return m.exp(-(current_time-delay)**2/(2.0 * source_width**2))

#Model
total_steps = int(size+delay)  # Time stepping
frame_interval = int(total_steps/20.0)
all_steps = np.linspace(0, size-1, size)


ez = np.zeros(size)
hy = np.zeros(size)
for time in xrange(total_steps+1):
    e_left = ez[1]
    e_right = ez[-2]
    hy[:-1] = hy[:-1] + (ez[1:] - ez[:-1])/imp0
    ez[1:] = ez[1:] + (hy[1:]-hy[:-1])*imp0
    ez[source_x] += source(time, delay, source_width)
    # ABC
    ez[0] = e_left
    ez[-1] = e_right
    # Output
    if time % frame_interval == 0 or time+5 > total_steps:
        plt.clf()
        plt.title("Ez after t=%i"%time)
        plt.plot(all_steps, ez)
        plt.show()
        plt.clf()
        plt.title("Hy after t=%i"%time)
        plt.plot(all_steps, hy*imp0)
        plt.show()




