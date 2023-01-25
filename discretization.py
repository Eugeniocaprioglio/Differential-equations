import numpy as np
vmn = [
     [1, 2, 3, 11],
     [4, 5, 6, -2],
     [7, 8, 9, 4],
     [1, 1, 5, 3]]

print(vmn[1][2])



# V_m,n-1

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

def matrix_vector(mat, t):
    vec_of_mat = np.array(t)
    arr = [0 for i in range(len(t))]
    
    
    for i in range(len(t)):
        for j in range(len(t)):
            arr[i] += vec_of_mat[j]*mat[i][j]
    return arr

a = [[
     [1, 2, 3, 11],
     [4, 5, 6, -2],
     [7, 8, 9, 4],
     [1, 1, 5, 3]],
     [
     [1, 2, 3, 11],
     [4, 5, 6, -2],
     [7, 8, 9, 4],
     [1, 1, 5, 3]],
     [
     [1, 2, 3, 11],
     [4, 5, 6, -2],
     [7, 8, 9, 4],
     [1, 1, 5, 3]]]

b = [[1, 2, -1],
     [-3, 2, 1],
     [0, 2, 0]]
    



print(right_v(vmn))
print(left_v(vmn))
print(up_v(vmn))
print(down_v(vmn))
print(diffusion_matrix(vmn))

print(matrix_vector(b, a))