import numpy as np


#-------In this part we build a function that returns the diffusion matrix-----------


def right_v(t):

    v = [[0 for i in range(len(t))] for j in range(len(t))]
    
    for i in range(len(t)):
        for j in range(len(t)-1):
            v[i][j+1] = t[i][j]

    return v

def left_v(t):

    v = [[0 for i in range(len(t))] for j in range(len(t))]
    
    for i in range(len(t)):
        for j in range(len(t)-1):
            v[i][j] = t[i][j+1]

    return v


def down_v(t):

    v = [[0 for i in range(len(t))] for j in range(len(t))]
    
    for i in range(len(t)-1):
        for j in range(len(t)):
            v[i+1][j] = t[i][j]

    return v

def up_v(t):

    v = [[0 for i in range(len(t))] for j in range(len(t))]
    
    for i in range(len(t)-1):
        for j in range(len(t)):
            v[i][j] = t[i+1][j]

    return v

def diffusion_matrix(t):
    return np.array(right_v(t)) + np.array(left_v(t)) + np.array(up_v(t)) + np.array (down_v(t))


#------- In this part we build a function that multiplies a matrix of scalars by an array of matrices------------------


def matrix_vector(mat, t):
    vec_of_mat = np.array(t)
    arr = [0 for i in range(len(t))]
    
    
    for i in range(len(t)):
        for j in range(len(t)):
            arr[i] += vec_of_mat[j]*mat[i][j]
    return arr

    


