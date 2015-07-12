# fdtd-1d
Example task for lectures "Computations in Physics" in ITMO University, spring 2015

## Step 0

Simulate a system with magnetic mirror boundary condition (H=0) on one side and electric mirror (E= 0) on the other side. The source is a gaussian profile propagating to boundaries and back, source located exactly in the center of the simulated domain. The successful presentation should provide a sequence of images as a time evolution (or animation) for electric and magnetic fields. The simulation should finish at the moment where the electric field is vanished (all energy is in the magnetic field).

## Step 1 and 2

Provide absorbing (Mur ABC) and PML (simplified CPML) boundary
condition. Compare ABC with  5,10, and 20 cell PML.

## Step 4

Compare against Fresnel equations http://en.wikipedia.org/wiki/Fresnel_equations , find the limits of FDTD applicability. 

## Step 5 and 6

Compare against single dielectric slab (e.g http://www.ece.rutgers.edu/~orfanidi/ewa/ch05.pdf), you should privide simulation of reflectionless cases   of a quarter-wavelength and half-wavelength slab width cases.


