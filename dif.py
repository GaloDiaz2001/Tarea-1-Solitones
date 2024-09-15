import numpy as np
import math as math
def right1(f,h):
    f_1 = np.zeros(len(f))
    for i in range(len(f)-2):
        f_1[i] = 1/(2*h) * (-f[i+2]+4*f[i+1]-3*f[i])
    f_1[len(f)-1] = f_1[len(f)-2] =f_1[len(f)-3]
    return f_1
def right2(f,h):
    f_1 = np.zeros(len(f))
    for i in range(len(f)-3):
        f_1[i] = 1/(h**2)*(-f[i+3] + 4*f[i+2] - 5*f[i+1] + 2*f[i])
    f_1[len(f)-1] = f_1[len(f)-2] =f_1[len(f)-3] = f_1[len(f)-4]
    return f_1
def central1(f,h):
    f_1 = np.zeros(len(f))
    for i in range (len(f)):
        if (i != 0 and i != len(f)-1):
            f_1[i] = 1/(2*h) * (f[i+1]-f[i-1])
    f_1[0] = f_1[1]
    f_1[len(f)-1] = f_1[len(f)-2]
    return f_1
def central2(f,h):
    f_1 = np.zeros(len(f))
    for i in range(len(f)):
        if (i != 0 and i != len(f)-1):
            f_1[i]=1/(h**2)*(f[i+1]-2*f[i] + f[i-1])
    f_1[0] = f_1[1]
    f_1[len(f)-1] = f_1[len(f)-2]
    return f_1
