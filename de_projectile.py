import math	
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt


#-----Define the differential equation-----------------------

def projectile(s, t, g, wind, m, mu, p):

    x, y, vx, vy = s

    return [
        vx,
        vy,
        -(mu/m)*(vx - wind)*pow((vx - wind)**2 + vy**2, (p - 1)/2),
        -g - (mu/m)*(vy)*pow((vx - wind)**2 + vy**2, (p - 1)/2)
    ]


#-----Set the initial conditions-----------------------------

x0 = 0
y0 = 0
v0 = 80
theta = 30

v0x = v0*math.cos(math.radians(theta))
v0y = v0*math.sin(math.radians(theta))

s0 = [x0, y0, v0x, v0y]

#------Set the constants of the model-------------------------

g = 9.81
wind = 0
m = 200
mu = 0.5
p = 1.5

def wind_direction(w):

    if (w > 0):
        return str(w) + 'm/s tailwind'
    return str(-w) + 'm/s headwind'

#-----Set the timespan and the discretization points----------

initial_time = 0
final_time = 7.9
discretization_points = 5000

time = np.linspace(initial_time, final_time, discretization_points)


#------Call odeint in order to find the solution of the ivp----

ideal_solution = odeint(projectile, s0, time, args=(g, 0, m, 0, p))
real_solution = odeint(projectile, s0, time, args=(g, wind, m, mu, p))

#------Unpacking the solutions----------------------------------

ideal_x = ideal_solution[:,0]
ideal_y = ideal_solution[:,1]
ideal_vx = ideal_solution[:,2]
ideal_vy = ideal_solution[:,3]

real_x = real_solution[:,0]
real_y = real_solution[:,1]
real_vx = real_solution[:,2]
real_vy = real_solution[:,3]

#------Define the max range function---------------------------

def max_range(t, x, y):

    for i in range(len(y)-1):
        if (y[i]*y[i + 1] < 0):
            return 'Max range: ' + str(round(x[i], 2)) + 'm in t = ' + str(round(t[i], 2)) + 's'
    
    return 'Max range out of timespan'


#------Define the max height function--------------------------

def max_height(x, y, vy):

    for i in range(len(vy)-1):
        if (vy[i]*vy[i + 1] < 0):
            return 'Max height: ' + str(round(y[i], 1)) + 'm at x = ' + str(round(x[i], 2)) + 'm'

#------Define the kinetic energy function----------------------

def kinetic_energy(vx, vy):

    arr = []
    for i in range(len(vx)):
        arr.append(0.5*m*(vx[i]**2 + vy[i]**2))
    return arr


#------Define the potential energy function-------------------

def potential_energy(h):

    e_ref = min(m*g*h[0], m*g*h[-1]) #reference potential energy

    arr = []
    for i in range(len(h)):
        arr.append(m*g*h[i] - e_ref)
    return arr

def mechanical_energy(h, vx, vy):

    e_ref = min(m*g*h[0], m*g*h[-1]) #reference potential energy
    arr = []
    for i in range(len(vx)):
        arr.append(0.5*m*(vx[i]**2 + vy[i]**2) + m*g*h[i] - e_ref)
    return arr


#------Plot the trajectory------------------------------------

def plot_trajectory():

    axis = plt.subplot()
    axis.set_ylabel('Height (m)')
    axis.set_xlabel('Range (m)')
    axis.plot(real_x, real_y, label='Real trajectory ' + max_range(time, real_x, real_y) + '\n' + max_height(real_x, real_y, real_vy))
    axis.plot(ideal_x, ideal_y, label='Ideal trajectory '+ max_range(time, ideal_x, ideal_y) + '\n' + max_height(ideal_x, ideal_y, ideal_vy))
    plt.text(50 , 10, r'$ \theta = $' + str(theta) + '°' + '\n' + '$V_0 = $ ' + str(v0) + 'm/s' + '\n' + wind_direction(wind))
    plt.legend()
    plt.grid()
    plt.show()

    return




#------Plot the energy---------------------------------------

def plot_energy():

    axis = plt.subplot()
    axis.set_ylabel('Energy (J)')
    axis.set_xlabel('Time (s)')
    axis.plot(time, kinetic_energy(real_vx, real_vy), label='Kinetic energy')
    axis.plot(time, potential_energy(real_y), label='Potential energy')
    axis.plot(time, mechanical_energy(real_y, real_vx, real_vy), label='Mechanical energy')
    plt.legend()
    plt.grid()
    plt.show()

    return

plot_energy()