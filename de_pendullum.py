import math	
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def pendullum(y, t, g, l, mu):

    theta, omega = y

    return  [omega, 
            -(g/l) * math.sin(theta) - mu*l*omega
            ]
    

y0 = [0.3*math.pi, 0]

m = 500
g = 9.81
l = 0.4
mu = 0.2

initial_time = 0
final_time = 4
discretization_points = 5000

timePoints = np.linspace(initial_time, final_time, discretization_points)

solution = odeint(pendullum, y0, timePoints, args=(g, l, mu)) 

def potential_energy(array):

    arr = []
    for i in range(len(array)):
        arr.append(m*g*l*(1 - math.cos(array[i])))
    return np.array(arr)

def kinetic_energy(array):

    arr = []
    for i in range(len(array)):
        arr.append(0.5*m*(l*array[i])**2)
    return np.array(arr)



axis = plt.subplot()

axis.plot(timePoints, potential_energy(solution[:,0]) + kinetic_energy(solution[:,1]))
axis.plot(timePoints, kinetic_energy(solution[:,1]))

plt.show()

