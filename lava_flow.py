"""
Lava Flowing Down
This code solve the 1D equation df/dt = d^2f/dx^2 + g*sin(alpha), which describe the velocity field 
of lava flowing down an incline plane. The output is a simulation of the time evolution of the 
velocity field. In the box is also plotted the analytical solution. 

@author: Nicolas Desjardins 
@collab: This script was heavily influenced by Prof. Lee's script : phys432_w2022_diff.py 
10/02/2022 (DD/MM/YY)
"""

#Basic idea of the simulation. The x-axis is perpendicular to the fluid flow, the plane is on the
#left edge of the simulation box and the air-lava interface is on the right edge of the simulation box.



import numpy as np
import matplotlib.pyplot as plt

Ngrid = 60 #grid lenght 
Nsteps = 12000 #number of time streps 
dt = 0.062 #time steps
dx = 1 #grid spacing
g = 0.5 # g*sin(alpha)

D = 16 # Diffusion constatn
beta = D*dt/dx**2

x = np.arange(0, Ngrid*1., dx)  # x-axis has lenght Ngrid and spacing dx
                                # multiplying by 1. to make sure this is an array of floats not integers

f = np.zeros(Ngrid) #Initializing the the velocity grid as 0 everywhere

# Set up plot
plt.ion()
fig = plt.figure()
ax = fig.add_subplot()
plt.xlim([0,Ngrid]) # x-limit corrspond to the height of the lava
plt.ylim([0,1/2*g*Ngrid**2/D + 0.05*1/2*g*Ngrid**2/D]) #set up the the y limit as 1.05% the analytical
                                                        #soln at H = Ngrid for nice visual

line_sim, = ax.plot(x, f) #plotting the intial velocity grid
line_analytic, = ax.plot(x, -g/D*(1/2*x**2-Ngrid*x)) #plotting the analytical soln 

for ct in range(Nsteps):

    ## Calculate diffusion first
    # Setting up matrices for diffusion operator
    A = np.eye(Ngrid) * (1.0 + 2.0 * beta) + np.eye(Ngrid, k=1) * -beta + np.eye(Ngrid, k=-1) * -beta
       
    # No slip boundary condition
    A[0][0] = 1.0 
    A[0][1] = 0 

    #Stress-free boundary conditon 
    A[Ngrid - 1][Ngrid - 1] = 1.0 + beta

    #Solve for the next velocity grid
    f = np.linalg.solve(A, f)

    f += dt*g # Add the graviational contribution
    f[0] = 0 #Inforce 0 at H=0 

    # update the plot
    line_sim.set_ydata(f)

    plt.show()    
    plt.pause(0.00001)