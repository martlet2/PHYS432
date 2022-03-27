"""
This produce an animation of the evolution of a adiabtic strong shock in 1D.
The animation is composed of 2 panels a top one which which is the density
in the 1D array over time. The bottom is the mach number in the 1D array
over time.

@author: Nicolas Desjardins
@collab: Mattias Ladza, Jules Faucher
March 26th 2022
"""
from shutil import which
import numpy as np
import matplotlib.pyplot as pl

# Set up the grid, time and grid spacing, and the sound speed squared
Ngrid = 100 #number of points in x 
Nsteps = 5000 #number of time steps 
dt = 0.005#temporal division 
dx = 2.5 #spacial division 

x = np.arange(Ngrid) * dx # grid
f1 = np.ones(Ngrid) # rho
f2 = np.zeros(Ngrid) # rho x u
f3 = np.ones(Ngrid) #rho x etot
gamma  = 5/3
u = np.zeros(Ngrid+1) # advective velocity (keep the 1st and last element zero)

def advection(f, u, dt, dx):
    # calculating flux terms
    J = np.zeros(len(f)+1) # keeping the first and the last term zero
    J[1:-1] = np.where(u[1:-1] > 0, f[:-1] * u[1:-1], f[1:] * u[1:-1])
    f = f - (dt / dx) * (J[1:] - J[:-1]) #update

    return f

# Apply initial Gaussian perturbation
Amp, sigma = 100000 , Ngrid/10 
f3 = f3 + Amp * np.exp(-(x - x.max()/2) ** 2 / sigma ** 2)

#setting pressure, cs2 and mach number 
P = (gamma-1)/gamma * (f3-0.5*f2**2/f1) #(gamma-1) * (f3-0.5*f2**2/f1)
cs2 = gamma*P/f1
mach = (f2/f1)/np.sqrt(cs2)

# plotting
pl.ion()
fig, ax = pl.subplots(2,1)


#This block is for the top panel wich plots f1 (rho) 
x1, = ax[0].plot(x, f1)
analytic, = ax[0].plot(x,np.full(Ngrid,4))
ax[0].set_xlim([0, dx*Ngrid+1])
ax[0].set_ylim([0, 5])
ax[0].set_xlabel('x')
ax[0].set_ylabel('Density')
#ax[0].grid(b=True, which="major")
#ax[0].grid(b=True, which="minor")
#ax[0].minorticks_on()

#This block is for the botto panel wich plots the Mach number 
x2, = ax[1].plot(x, mach)
analytic, = ax[1].plot(x,np.full(Ngrid,4))
ax[1].set_xlim([0, dx*Ngrid+1])
ax[1].set_ylim([-3, 3])
ax[1].set_xlabel('x')
ax[1].set_ylabel('Mach number')



fig.canvas.draw()

for ct in range(Nsteps):

    # 1) advection velocity at the cell interface and at the simulation boundaries
    u[1:-1] = 0.5 * ((f2[:-1] / f1[:-1]) + (f2[1:] / f1[1:]))


    # 2) update density and momentum
    f1 = advection(f1, u, dt, dx)
    f2 = advection(f2, u, dt, dx)

    # 3.1) compute presssure and sound speed square
    P = (gamma-1)/gamma * (f3-0.5*f2**2/f1)
    cs2 = gamma*P/f1

    # 3.2) apply the pressure gradient force to momentum equation 
    f2[1:-1] +=  - 0.5*(dt/dx)*(P[2:]-P[:-2])

    # 3.3) apply the correct boundary condition for the source term
    f2[0] += - 0.5*(dt/dx)*(P[1]-P[0])
    f2[-1] += - 0.5*(dt/dx)*(P[-1]-P[-2])

    # 4 ) advection velocity at the cell interface
    u[1:-1] = 0.5 * ((f2[:-1] / f1[:-1]) + (f2[1:] / f1[1:]))


    #5)Advect energy 
    f3 = advection(f3, u, dt, dx)

    #6.1) Recompute pressure
    P = (gamma-1)/gamma* (f3-0.5*f2**2/f1) #(gamma-1) * (f3-0.5*f2**2/f1)

    #6.2) Apply corresponding source term to the energy equation
    Pu = P*f2/f1
    f3[1:-1] += - 0.5*(dt/dx)*(Pu[2:]-Pu[:-2])

    # 6.3) apply the correct boundary condition for the source term
    f3[0] += - 0.5*(dt/dx)*(Pu[1]-Pu[0])
    f3[-1] += - 0.5*(dt/dx)*(Pu[-1]-Pu[-2])

    P = (gamma-1)/gamma * (f3-0.5*f2**2/f1) #(gamma-1) * (f3-0.5*f2**2/f1)
    cs2 = gamma*P/f1
    mach = (f2/f1)/np.sqrt(cs2)



    # update the plot
    x1.set_ydata(f1)
    x2.set_ydata(mach)
    fig.canvas.draw()
    pl.pause(0.0001)