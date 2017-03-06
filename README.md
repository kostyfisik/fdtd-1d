# fdtd-1d

Example task for lectures "Computations in Physics" in ITMO
University, spring (2015-2017)

[Final report in *.pdf format](https://github.com/kostyfisik/fdtd-1d/blob/master/final-report.pdf)

Recommended books:

(easy) Understanding the Finite-Difference Time-Domain Method, John B. Schneider, www.eecs.wsu.edu/~schneidj/ufdtd, 2010. 
(it is also available at GitHub
https://github.com/john-b-schneider/uFDTD )

(medium) Numerical electromagnetics : the FDTD method / Umran S. Inan, Robert A. Marshall. 2011


(full) A. Taflove and S. C. Hagness, Computational Electrodynamics: The Finite-Difference Time-Domain Method, 3rd  ed.  Norwood, MA: Artech House, 2005. 



## Step 0

Simulate a system with magnetic mirror boundary condition (H=0) on one
side and electric mirror (E= 0) on the other side. The source is a
Gaussian profile propagating to boundaries and back, source located
exactly in the center of the simulated domain. The successful
presentation should provide a sequence of images as a time evolution
(or animation) for electric and magnetic fields. The simulation should
finish at the moment where the electric field is vanished (all energy
is in the magnetic field).

## Step 1,2, and 3

Provide absorbing (Mur ABC) and PML (simplified CPML) boundary
condition. Compare ABC with  5,10, and 20 cell PML.

## Step 4

Compare against Fresnel equations
http://en.wikipedia.org/wiki/Fresnel_equations , find the  limits of
FDTD applicability.

## Step 5 and 6

Compare against single dielectric slab (e.g
http://www.ece.rutgers.edu/~orfanidi/ewa/ch05.pdf),  you should
provide simulation of reflection-less cases   of a quarter-wavelength
and half-wavelength slab width cases.

## Bonus script

You can put all the graphical results to "fig" directory,
generate-report.py will list all the files in it and generate LaTeX
presentation (using Beamer package). After that it should be easy to
provide a single file report with all the figures (do not forget to
provide reasonable captions to the figures)

