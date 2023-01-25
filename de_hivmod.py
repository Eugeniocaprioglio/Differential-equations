import math	
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

max_dict = {
    "r": [0, 100], 
    "n": [1, 300], 
    "alpha": [2, 1.5], 
    "gamma": [3, 0.001], 
    "s_max": [4, 1500], 
    "mus": [5, 0.1],
    "mui": [6, 0.5],
    "muv": [7, 10]
    }

def gen(a,k):
    return abs(math.sin(a**2*+3*k**3-a*k**4+36*a+k+12))

def values(a):

    # the global declaration means that the variable can be used outside the function by another function.
    global alpha
    global mus
    global r
    global s_max
    global gamma
    global mui
    global muv
    global n

    #max_dict["alpha"][1] is the maximim value of alpha and max_dict["alpha"][0] is the index k
    alpha = max_dict["alpha"][1]*gen(a, max_dict["alpha"][0])
    mus = max_dict["mus"][1]*gen(a, max_dict["mus"][0])
    r = max_dict["r"][1]*gen(a, max_dict["r"][0])
    s_max = max_dict["s_max"][1]*gen(a, max_dict["s_max"][0])
    gamma = max_dict["gamma"][1]*gen(a, max_dict["gamma"][0])
    mui = max_dict["mui"][1]*gen(a ,max_dict["mui"][0])
    muv = max_dict["muv"][1]*gen(a ,max_dict["muv"][0])
    n = max_dict["n"][1]*gen(a, max_dict["n"][0])
    return #this function doesn't return any specific values, but set the values of several global variables.


def uninfected_steady_state(a):

    values(a)

    return [
        (r - mus + math.sqrt((r - mus)**2 + (4*alpha*r/s_max)))*s_max/2/r, 
        0, 
        0
    ]

def infected_steady_state(a):

    values(a)

    return [
        muv/(gamma*n), 
        alpha/mui - (mus*muv)/(gamma*mui*n) + ((muv*r)/(gamma*mui*n))*(1-(muv/(gamma*n*s_max))), 
        alpha*n/muv - mus/gamma + ((r/gamma)*(1-(muv/(gamma*n*s_max))))
    ]

def reproduction_ratio(a):
    
    values(a)

    return gamma*n*uninfected_steady_state(a)[0]/muv

def n_parameters(a):

    values(a)
    
    return[
        (gamma**2)*mui*s_max - 4*alpha*(mui + muv),
        - 2*gamma*mui*muv*((mui**2 + 3*mui*muv + muv**2)*s_max - 4*alpha*(mui + muv)),
        muv*(mui**4 + 2*(mui**2)*muv*(3*mui - 2*mus) + 6*mui*(muv**3) + mui*(muv**2)*(11*mui - 4*mus) + muv**4)
    ]

def r_parameters(a):

    values(a)

    return [
        (muv**4)*(mui + muv),
        gamma*(muv**2)*n*s_max*(-gamma*mui*muv*n*s_max + 2*alpha*gamma*mui*n + 2*alpha*gamma*muv*n + (mui**2)*muv + muv**3),
        (n**2)*(gamma**2)*(s_max**2)*(alpha*gamma*(muv**3)*n + alpha*gamma*mui*muv*n*(mui + muv) + (alpha**2)*(gamma**2)*(n**2)*(mui + muv) + mui*mus*muv**3),
    ]

def state_is_in_P(a):

    values(a)

    [an, bn, cn] = n_parameters(a)
    [ar, br, cr] = r_parameters(a)
    root_1 = (-bn + (bn**2 -4*an*cn)**0.5)/(2*an)
    root_2 = (-bn - (bn**2 -4*an*cn)**0.5)/(2*an)

    return n > max(root_1, root_2) and ar*r**2+br*r+cr < 0 #the last statement is equivalent to have the value of r between r1 and r2

def theorem_1(a):

    if (reproduction_ratio(a) < 1): 
        return "The uninfected steady state is the unique physically relevant steady state and it is asymptotically stable." 
    elif (state_is_in_P(a)):
        return "The uninfected and the infected steady states are both physically relevant" + '\n' + " and the infected steady state is unstable"

    return "The uninfected and the infected steady states are both physically relevant "+ '\n' " and the the infected steady state is asymptotically stable"

# ------------------- Defining the function that returns dy/dt -------------------

def virus_propagation(y, t, re, nv, al, g, sm, mu_s, mu_i, mu_v):
    
    s, i, v = y #unpack
    return [
        al - mu_s*s + re*s*(1 - s/sm) - g*v*s,
        g*v*s - mu_i*i,
        nv*mu_i*i - mu_v*v
    ]

# ----------------  Call the function that set the parameters----------------------

parameter = 6
values(parameter)

# ---------------- Setting the initial conditions-------------------------------

inititial_healthy_tcells = 230
inititial_infected_tcells = 200
initial_virus = 1000

y0 = [inititial_healthy_tcells, inititial_infected_tcells, initial_virus]

# ----------------Setting the time and discretization points ---------------------

initial_time = 0
final_time = 200
discretization_points = 5000 

timePoints = np.linspace(initial_time, final_time, discretization_points)



solution = odeint(virus_propagation, y0, timePoints, args=(r, n, alpha, gamma, s_max, mus, mui, muv))

# unpacking solutions
s_t = solution[:,0]
i_t = solution[:,1]
v_t = solution[:,2]

def dynamic_plot():

    axis = plt.subplot()

    axis.plot(timePoints, s_t, color='blue', label='Healthy CD4+ cells')
    axis.plot(timePoints, i_t, color='green', label='Infected CD4+ cells')
    axis.plot(timePoints, v_t, color='magenta', label='Virus')
    axis.set_title('Dinamic behaviour')
    
    plt.text(0 , 0, theorem_1(parameter))

    plt.legend()
    plt.tight_layout()
    plt.grid()
    plt.show()


def phase_plot(a, b, lab):

    axis = plt.subplot()

    axis.plot(a, b, label=lab)
   

    plt.legend()
    plt.grid()
    plt.show()

dynamic_plot()
phase_plot(i_t, v_t, 'Infected T cells vs Virus')