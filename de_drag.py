import math	
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

#---------Constants in S. I. units---------
g = 9.81 

#---------List of classes------------------

class Fluid:

    def __init__(self, density, viscosity):

        self.density = density
        self.viscosity = viscosity

class Solid:

    def __init__(self, density):
        self.density = density
        
    def volume(self, diameter):
        return 4/3*math.pi*(diameter/2)**3

    def area(self, diameter):     
        return math.pi*(diameter/2)**2

    def mass(self, diameter):
        return self.density * self.volume(diameter)

#-------List of fluids----------------------

water = Fluid(density=1000, viscosity=0.001)
air = Fluid(density=1.25, viscosity=6e-5)
benzene = Fluid(density=876, viscosity=6e-4)
acetone = Fluid(density=787, viscosity=3.16e-4)
ethyl_alcohol = Fluid(density=787, viscosity=0.001)
methyl_alcohol = Fluid(density=789, viscosity=5.6e-4)
glycerine = Fluid(density=1258, viscosity=0.96)
ethilenglycol = Fluid(density=1100, viscosity=0.0162)
heavy_oil = Fluid(density=906, viscosity=0.107)

#------List of solids-------------------------

iron = Solid(density=8700)
glass = Solid(density=2380)
granite = Solid(density=2900)
sand = Solid(density=2500)
cement = Solid(density=1900)
aluminium = Solid(density=2730)

chosen_fluid = glycerine
chosen_solid = glass

def reynolds(d, v, ro, mu):

    return d*abs(v + 0.00000001)*ro/mu

def drag_coefficient(reynolds):

    return 0.4 + 24/reynolds + 6/(1 + math.sqrt(reynolds))

def height(s, t, fluid, solid, diameter):

    y, vel = s
    
    re = reynolds(diameter, vel, fluid.density, fluid.viscosity)
    cd = drag_coefficient(re)
    dens = fluid.density
    area = solid.area(diameter)
    vol = solid.volume(diameter)
    m = solid.density*vol

    return [
        vel,
        -g - 0.5*dens*(cd/m)*area*vel*abs(vel) + g*(dens/m)*vol
    ]

# set the initial condition

initial_height = 0
initial_velocity = 100

y0 = [initial_height, initial_velocity]

# define the discretization points

duration = 10
discretization_points = 500
time = np.linspace(0, duration, discretization_points)

# call odeint() function

solutionOde = odeint(height, y0, time, args=(chosen_fluid, chosen_solid, 1))

def max_height(h,v):
    for i in range(len(v)-1):
        if (v[i]*v[i+1]<0):
            return 'Max height: ' + str(round(h[i], 2))
    return 'Max height: ' + str(round(max(h[0], h[len(h)-1]),2))

def floor(t, h, v):
    for i in range(len(h) - 1):
        if (h[i]*h[i+1]<0):
            return 'Reach the floor at ' + str(round(t[i], 2)) + 's.' + '\n' + 'with a final speed of: ' + str(abs(round(v[i], 2))) +'m/s.'
    return ' Zero height out of time span' + '\n' + 'Final speed of: ' + str(abs(round(v[len(v)-1], 4))) +'m/s.'

def plot_trajectory():
    fig, ax1 = plt.subplots()

    color = 'tab:red'
    ax1.set_xlabel('time (s)')
    ax1.set_ylabel('height (m)', color = color)
    ax1.plot(time, solutionOde[:, 0], color = color, label=max_height(solutionOde[:,0],solutionOde[:,1]) + '\n' + floor(time, solutionOde[:,0],solutionOde[:,1]))
    ax1.tick_params(axis='y', labelcolor=color)
    ax2 = ax1.twinx() # creates a twin axis that shares the time values.

    color = 'tab:green'
    ax2.set_ylabel('velocity (m/s)', color = color)
    ax2.plot(time, solutionOde[:, 1], color = color)
    ax2.tick_params(axis='y', labelcolor=color)
    plt.grid()
    #fig.tight_layout()
    fig.legend()
    plt.show()






