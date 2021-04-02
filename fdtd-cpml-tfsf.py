import numpy as np

# space and time discretization
dx = 10e-9
dt = 5e-18

# number of x points in X
N = 400

# constants
epsilon0 = 8.8541878128e-12
mu0 = 1.256637062e-6
c0 = 2.99792458e8

# index of background
index_background = 1.0

# set CPML
thick_pml = 24
P_m = 3.0 # polynomail
kappa_max = 16.0 # stretch coordiantes

# pulse
tau_F = 2e-15 # pulse duration, FWHM (in second)

sigma_max_Ez = 0.8 * (P_m+1)/(mu0*c0*dx*index_background)
alpha_x_max = c0 * epsilon0 / 1e-7

# init field 
Ez = np.zeros(N)
Hy = np.zeros(N)
# for CPML
Psi_Ezx = np.zeros(N)
Psi_Hyx = np.zeros(N)

# coordiantes
xx = np.arange(-N*dx/2, N*dx/2, dx)

# ---------- E-field coefficients
Ca_z = np.ones(N)
Cb_z = np.ones(N)*dt/dx/index_background/index_background/epsilon0

# ---------- H-field coefficients
Da_y = np.ones(N)
Db_y = np.ones(N)*dt/dx/mu0

# set material
index_material = 1.0
size_material = 20
for i in range(np.int32(N/2),np.int32(N/2)+size_material):
    Cb_z[i] = dt/dx/index_material/index_material/epsilon0

# --------------------------------
inv_kappa_Ez = np.ones(N)
b_pml_Ez = np.zeros(N)
c_pml_Ez = np.zeros(N)

inv_kappa_Hy = np.ones(N)
b_pml_Hy = np.zeros(N)
c_pml_Hy = np.zeros(N)


# now set CPML updating parameters
for i0 in range(1,thick_pml):
    poly_i = (1.0 - (i0 - 1.0)/ thick_pml)**P_m
    sigma_Ez_i = sigma_max_Ez * poly_i
    kappa_Ez_i = 1.0 + (kappa_max-1.0)*poly_i
    alpha_Ez_i = alpha_x_max * ((i0-1.0)/thick_pml)**P_m
    
    print (i0, "--",sigma_Ez_i, kappa_Ez_i, alpha_Ez_i)
    
    inv_kappa_Ez[i0] = 1.0 / kappa_Ez_i
    b_pml_Ez[i0] = np.exp(-(sigma_Ez_i/epsilon0/kappa_Ez_i+alpha_Ez_i/epsilon0)*dt)
    c_pml_Ez[i0] = (b_pml_Ez[i0] - 1.0)* (sigma_Ez_i/(sigma_Ez_i*kappa_Ez_i + kappa_Ez_i*kappa_Ez_i*alpha_Ez_i))
    
    print ((sigma_Ez_i/epsilon0/kappa_Ez_i+alpha_Ez_i/epsilon0)*dt, b_pml_Ez[i0], c_pml_Ez[i0] )
    
    inv_kappa_Hy[i0] = 1.0 / kappa_Ez_i
    b_pml_Hy[i0] = np.exp(-(sigma_Ez_i/epsilon0/kappa_Ez_i+alpha_Ez_i/epsilon0)*dt)
    c_pml_Hy[i0] = (b_pml_Hy[i0] - 1.0)* (sigma_Ez_i/(sigma_Ez_i*kappa_Ez_i + kappa_Ez_i*kappa_Ez_i*alpha_Ez_i))

    # right hand side
    inv_kappa_Ez[N-i0] = inv_kappa_Ez[i0]
    b_pml_Ez[N-i0] = b_pml_Ez[i0]
    c_pml_Ez[N-i0] = c_pml_Ez[i0]
    
    inv_kappa_Hy[N-i0] = inv_kappa_Hy[i0]
    b_pml_Hy[N-i0] = b_pml_Hy[i0]
    c_pml_Hy[N-i0] = c_pml_Hy[i0]

    #print ("set CPML done.")
    
def updateE_Field(Ez):
    
    for i in range(1, N-1):
        Ez[i] = Ca_z[i] * Ez[i] + Cb_z[i] * (Hy[i] - Hy[i-1]) * inv_kappa_Ez[i]
            
    #print ("updating field")
    
def updateH_Field(Hy):
    
    for i in range(1, N-1):
        Hy[i] = Da_y[i] * Hy[i] + Db_y[i] * (Ez[i+1] - Ez[i]) * inv_kappa_Hy[i]
            
    #print ("updating field")
    
def updateCPML():
    # left side
    for i0 in range(1, thick_pml+1):
        Psi_Hyx[i0] = b_pml_Hy[i0] * Psi_Hyx[i0] + c_pml_Hy[i0] * (Ez[i0+1]-Ez[i0])
        Psi_Ezx[i0] = b_pml_Ez[i0] * Psi_Ezx[i0] + c_pml_Ez[i0] * (Hy[i0]-Hy[i0-1])
        Ez[i0] += Cb_z[i0] * Psi_Ezx[i0]
        Hy[i0] += Db_y[i0] * Psi_Hyx[i0]
        
    # right side
    for i0 in range(N-thick_pml-1, N-1):
        Psi_Hyx[i0] = b_pml_Hy[i0] * Psi_Hyx[i0] + c_pml_Hy[i0] * (Ez[i0+1]-Ez[i0])
        Psi_Ezx[i0] = b_pml_Ez[i0] * Psi_Ezx[i0] + c_pml_Ez[i0] * (Hy[i0]-Hy[i0-1])
        Ez[i0] += Cb_z[i0] * Psi_Ezx[i0]
        Hy[i0] += Db_y[i0] * Psi_Hyx[i0]
        
    #print ("updating CPML")
    
    
def updateEContinousWave(time, omega_s, time_shift, i_source, sign):
    #Ez[i_source] -= (sign) * Cb_z[i_source] * np.exp(-(time - time_shift)**2.0/2.0/(8e-16)**2.0) / (c0*mu0/index_background)* np.sin(omega_s * (time - time_shift))
    # if (time - time_shift)>=0:
    #     Ez[i_source] -= (sign) * Cb_z[i_source] / (c0*mu0/index_background) * np.sin(omega_s * (time - time_shift))* np.exp(-2*np.log(2)*(time - time_shift)**2.0/(tau_F)**2.0)
    Ez[i_source] -= (sign) * Cb_z[i_source] / (c0*mu0/index_background) * np.sin(omega_s * (time - time_shift))* np.exp(-2*np.log(2)*(time - time_shift)**2.0/(tau_F)**2.0)

def updateHContinousWave(time, omega_s, time_shift, i_source, sign):
    #Hy[i_source] += (sign) * Db_y[i_source] * np.exp(-(time - time_shift)**2.0/2.0/(8e-16)**2.0)* np.sin(omega_s * (time - time_shift))
    # if (time - time_shift)>=0:
    #     Hy[i_source] += (sign) * Db_y[i_source] * np.sin(omega_s * (time - time_shift))* np.exp(-2*np.log(2)*(time - time_shift)**2.0/(tau_F)**2.0)
    Hy[i_source] += (sign) * Db_y[i_source] * np.sin(omega_s * (time - time_shift))* np.exp(-2*np.log(2)*(time - time_shift)**2.0/(tau_F)**2.0)
# ---------------------------------
# ------the main program ----------
import matplotlib.pyplot as plt

plt.ion()

fig0, ax0 = plt.subplots()
line0, = plt.plot(np.zeros(N),'g.-')
ax0.set_ylim([-1.99,1.99])
line_pml_left,  = plt.plot(thick_pml*np.ones(100), np.linspace(-2,2,100),'r--')
line_pml_right, = plt.plot((N-thick_pml-1)*np.ones(100), np.linspace(-2,2,100),'r--')

line_material_left,  = plt.plot((np.int32(N/2)-0)*np.ones(100), np.linspace(-1,1,100),'b--')
line_material_right, = plt.plot((np.int32(N/2)+size_material)*np.ones(100), np.linspace(-1,1,100),'b--')

def update_data(fig0,ax0,line0,data1d):
   
    line0.set_ydata(data1d)  
    #ax0.relim()
    #ax0.autoscale_view()
    fig0.canvas.draw()
    fig0.canvas.flush_events()

loc_source_left = np.int32(N/2) - 120
size_tfsf = 240
loc_source_right = loc_source_left + size_tfsf

for i00 in range(10000):
    
    print (str(i00)+"/"+str(1000))
    
    time0 = dt * i00
    wavelength = 500e-9
    omega0 = c0 / wavelength * 2.0 * np.pi
    # right side of TFSF source
    updateHContinousWave(time0, omega0, 5e-15+(size_tfsf)*dx*index_background/c0, loc_source_right, +1)
    # left side of TFSF source
    updateHContinousWave(time0, omega0, 5e-15, loc_source_left, -1)
    updateH_Field(Hy)
    # right side of TFSF source
    updateEContinousWave(time0, omega0, 5e-15+(size_tfsf)*dx*index_background/c0+dt/2+dx/2/c0*index_background, loc_source_right, +1)
    # left side of TFSF source
    updateEContinousWave(time0, omega0, 5e-15+dt/2+dx/2/c0*index_background, loc_source_left, -1)
    updateE_Field(Ez)
    
    updateCPML()
    if (i00%20==0):
        # display intensity
        #update_data(fig0, ax0, line0, np.abs(Ez))
        # uncomment line below to display Ez
        update_data(fig0, ax0, line0, Ez)
    