import math	
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def falling_to_moon(s, t, r, m, g):

 
   
    return [
        s[1],
        -g*m/(r + s[0])**2
    ]

#---initial values---

initial_height = 600000
initial_velocity = 0

y0 = [initial_height, initial_velocity]

#---constants-----

moon_mass = 7.342e22
moon_radius = 1.737e6
gravity_const = 6.67e-11

#----time---------

initial_time = 0
final_time = 1200
discretization_points = 10000

time = np.linspace(initial_time, final_time, discretization_points)

#-----call odeint----

sol = odeint(falling_to_moon, y0, time, args=(moon_radius, moon_mass, gravity_const))

#----unpacking-----

h = sol[:, 0]
vel = sol[:, 1]

#----plot!-------

fig, ax1 = plt.subplots()

color = 'tab:red'
ax1.set_xlabel('time (s)')
ax1.set_ylabel('height (m)', color = color)
ax1.plot(time, h, color = color)
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx() # creates a twin axis that shares the time values.

color = 'tab:green'
ax2.set_ylabel('velocity (m/s)', color = color)
ax2.plot(time, -vel, color = color)
ax2.tick_params(axis='y', labelcolor=color)
plt.grid()
fig.tight_layout()
fig.legend()

plt.show()

