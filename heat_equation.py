import math	
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

delta_x = 0.01 #meters
delta_t = 0.1 #seconds
alpha = 1.27e-4 #m^2/s

def calculate_fourier():

    if (alpha*delta_t/delta_x**2 > 0.25):
        return 'The system is unstable. Choose a smaler value of delta_T or bigger value of delta_X'
    return alpha*delta_t/delta_x**2

def right_v(t):

    v = [[0 for i in range(len(t[0]))] for j in range(len(t))]
    
    for i in range(len(t)):
        for j in range(len(t[0])-1):
            v[i][j+1] = t[i][j]

    return v

def left_v(t):

    v = [[0 for i in range(len(t[0]))] for j in range(len(t))]
    
    for i in range(len(t)):
        for j in range(len(t[0])-1):
            v[i][j] = t[i][j+1]

    return v


def down_v(t):

    v = [[0 for i in range(len(t[0]))] for j in range(len(t))]
    
    for i in range(len(t)-1):
        for j in range(len(t[0])):
            v[i+1][j] = t[i][j]

    return v

def up_v(t):

    v = [np.zeros(len(t[0])) for j in range(len(t))]
    
    for i in range(len(t)-1):
        for j in range(len(t[0])):
            v[i][j] = t[i+1][j]

    return v

def surroundigs_matrix(t):
    return np.array(right_v(t)) + np.array(left_v(t)) + np.array(up_v(t)) + np.array (down_v(t))

def algorithm(fo, t):

    ts = surroundigs_matrix(t)

    return (1 - 4*fo)*np.array(t) + fo*np.array(ts)

def create_matrix(columns, raws): 

    mat = [[0 for i in range(columns)] for j in range(raws - 1)]
    mat.append([1 for i in range(columns)])
    

    return mat

def initial_conditions(t):
 
    for i in range(len(t)):
        t[i][0] = 0
        t[i][-1] = 0
    for i in range(len(t[0])):
        t[0][i] = 0
        t[-1][i] = 1
    
    return t

def heat_equation(u, fourier, k):

    for i in range(k):
        u = algorithm(fourier, u)
        u = initial_conditions(u)
    return u




def plot_grid(x_length, y_length, iterations):
    mat = create_matrix(x_length, y_length)
    fourier = calculate_fourier()
    if isinstance(fourier, float):
        plt.pcolormesh(heat_equation(mat, fourier, iterations))
        plt.show()
    else:
        print(fourier)
    return
def plot_profile_x(y):
    axis = plt.subplot()
    x = []
    x.append(i for i in range(len))

plot_grid(10, 20, 5)
print(type(calculate_fourier()))




