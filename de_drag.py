import math	
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def reynolds(d, v, ro, mu):

    return d*abs(v + 0.00000001)*ro/mu

def drag_coefficient(reynolds):

    return 0.4 + 24/reynolds + 6/(1 + math.sqrt(reynolds))

def height(s, t, g, m, diam, dens, visc):

    y, vel = s

    re = reynolds(diam, vel, dens, visc)
    cd = drag_coefficient(re)
    area = 0.25*math.pi*diam**2
    vol = 4/3*math.pi*(diam/2)**3

    return [
        vel,
        -g - 0.5*dens*(cd/m)*area*vel*abs(vel) + g*(dens/m)*vol
    ]

# constants in S. I. units

g = 9.81 

diameter = 0.1
density = 1.25
solid_density = 8700
viscosity = 1.5e-5
mass = solid_density*4/3*math.pi*(diameter/2)**3

final_velocity = (2*(diameter/2)**2/(9*viscosity))*(solid_density - density)*g

print('Final velocity = ' + str(round(final_velocity, 4)))



# set the initial condition

initial_height = 52
initial_velocity = 0

y0 = [initial_height, initial_velocity]

# define the discretization points

initial_time = 0
final_time = 3.5
discretization_points = 5000

timePoints = np.linspace(initial_time, final_time, discretization_points)

# call odeint() function and plt the results

solutionOde = odeint(height, y0, timePoints, args=(g, mass, diameter, density, viscosity))

def max_height(h,v):
    for i in range(len(v)-1):
        if (v[i]*v[i+1]<0):
            return 'Max height: ' + str(round(h[i], 2))
    return 'Max height: ' + str(round(max(h[0], h[len(h)-1]),2))

def floor(t, h, v):
    for i in range(len(h) - 1):
        if (h[i]*h[i+1]<0):
            return 'Reach the floor at ' + str(round(t[i], 2)) + 's.' + '\n' + 'with a final speed of: ' + str(abs(round(v[i], 2))) +'m/s. Re =' + str(round(diameter*density*-v[i]/viscosity))
    return ' Zero height out of time span' + '\n' + 'Final speed of: ' + str(abs(round(v[len(v)-1], 4))) +'m/s. Re =' + str(round(diameter*density*-v[len(v)-1]/viscosity,2))

fig, ax1 = plt.subplots()

color = 'tab:red'
ax1.set_xlabel('time (s)')
ax1.set_ylabel('height (m)', color = color)
ax1.plot(timePoints, solutionOde[:, 0], color = color, label=max_height(solutionOde[:,0],solutionOde[:,1]) + '\n' + floor(timePoints, solutionOde[:,0],solutionOde[:,1]))
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx() # creates a twin axis that shares the time values.

color = 'tab:green'
ax2.set_ylabel('velocity (m/s)', color = color)
ax2.plot(timePoints, solutionOde[:, 1], color = color)
ax2.tick_params(axis='y', labelcolor=color)
plt.grid()
#fig.tight_layout()
fig.legend()
plt.show()







