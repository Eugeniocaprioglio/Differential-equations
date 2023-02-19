import math	
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation

#------------Define the constants and calculate Fourier number-------------

delta_x = 0.01 #meters
delta_t = 0.01 #seconds
alpha = 1.27e-4 #m^2/s

def calculate_fourier():
    if (alpha*delta_t/delta_x**2 > 0.25):
        return 'The system is unstable. Choose a smaller value of delta_T or bigger value of delta_X'
    return alpha*delta_t/delta_x**2

def laplacian(t):
    
    up_v, down_v, right_v, left_v = [np.zeros((len(t), len(t[0]))) for i in range(4)]
    t = np.array(t)
    for i in range(len(t) - 1):
        up_v[i,:] = t[i+1,:]
        down_v[i+1,:] = t[i,:]
    for i in range(len(t[0]) - 1):
        right_v[:,i] = t[:,i+1]
        left_v[:,i+1] = t[:,i]
    return np.array(left_v) + np.array(right_v) + np.array(up_v) + np.array (down_v)

def create_matrix(columns, raws): 

    mat = np.full((raws, columns),0)
    for i in range(int(columns/3),columns - int(columns/3)):
        for j in range(int(raws/3),raws - int(raws/3)):
            mat[j][i] = 200
    return np.array(mat)

def initial_conditions(t):
 
    t[:,0] = t[:,1]
    t[:,-1] = 100
    t[0,:] = t[1,:]
    t[-1,:] = 200

    return t

def algorithm(fo, t):

    return (1 - 4*fo)*np.array(t) + fo*np.array(laplacian(t))

def heat_equation(u, fourier, k):

    for i in range(k):
        u = algorithm(fourier, u)
        u = initial_conditions(u)
    return u

#-------Plotting grid and temperature profiles------------

def plot_grid(x_length, y_length, iterations):
    mat = create_matrix(x_length, y_length)
    fourier = calculate_fourier()
    if isinstance(fourier, float):
        plt.pcolormesh(heat_equation(mat, fourier, iterations), cmap='inferno', vmax=200, vmin=0)
        plt.show()
    else:
        print(fourier)
    return

def animate(k):
    plot_grid(20,20,k)



def plot_yprofile(x_length, y_length, iterations, y_coord):
    mat = create_matrix(x_length, y_length)
    fourier = calculate_fourier()
    axis = plt.subplot()
    if isinstance(fourier, float):
        solution = heat_equation(mat, fourier, iterations)
        axis.plot([i for i in range(len(solution[:,y_coord]))],solution[:,y_coord])
        plt.show()
    else:
        print(fourier)
    
    return

def plot_xprofile(x_length, y_length, iterations, x_coord):
    mat = create_matrix(x_length, y_length)
    fourier = calculate_fourier()
    axis = plt.subplot()
    if isinstance(fourier, float):
        solution = heat_equation(mat, fourier, iterations)
        axis.plot([i for i in range(len(solution[x_coord,:]))],solution[x_coord,:])
        plt.show()
    else:
        print(fourier)
    
    return



