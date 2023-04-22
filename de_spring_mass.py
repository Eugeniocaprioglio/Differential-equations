import math	
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def spring_mass(s, t, l0_1, l0_2, k_1, k_2, m_1, m_2, b_1, b_2, kc):
    x_1, x_2, v_1, v_2 = s
    return[
        v_1,
        v_2,
        -(k_2/m_1)*(x_1 - x_2 + l0_2) - (k_1/m_1)*(x_1 - l0_1) - kc/(x_2 - x_1)**6 - (b_1/m_1)*v_1,
        -(k_2/m_2)*(x_2 - x_1 - l0_2) + kc/(x_2 - x_1)**6 - (b_2/m_2)*v_2
    ]

#-------Constants------------
k_contact = 1e-6
k1 = 500
k2 = 500
m1 = 20
m2 = 20
l01 = 0.2
l02 = 0.2
b1 = 0.5*m1
b2 = 0.5*m2

#------Initial conditions---------

x01 = 0.25
x02 = 0.5
v01 = 0
v02 = 0

s0 = [x01, x02, v01, v02]

#------Time----------------------

initial_time = 0
final_time = 40
discretization_points = 1000
time = np.linspace(initial_time, final_time, discretization_points)

#-------Call odeint---------------

solution = odeint(spring_mass, s0, time, args=(l01, l02, k1, k2, m1, m2, 0, 0, k_contact))

#---------Unpack solutions------------

x1_t = solution[:,0]
x2_t = solution[:,1]
v1_t = solution[:,2]
v2_t = solution[:,3]

#-------Plotting----------------------

def plot_trajectory():

    axis = plt.subplot()
    axis.plot(time, x1_t, color='blue')
    axis.plot(time, x2_t, color='green')

    plt.show()

def plot_phase(x, p):

    axis = plt.subplot()
    axis.plot(x, p, color='green')

    plt.show()

plot_phase(x1_t, m1*v1_t)
plot_phase(x2_t, m2*v2_t)


