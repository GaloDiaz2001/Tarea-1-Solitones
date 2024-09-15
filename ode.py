import numpy as np
import math as math
import matplotlib.pyplot as plt
import dif as dif

def RK4(f,y0,x0,xf,steps):
    # h val
    h = (xf-x0)/(steps-1)
    xval = x0 + np.arange(steps)*h
    yval = np.zeros(steps)
    
    y = y0
    for j, xi in enumerate(xval):
        yval[j] = y
        k0 = h*f(xi, y)
        k1 = h*f(xi+h/2, y+k0/2)
        k2 = h*f(xi+h/2, y+k1/2)
        k3 = h*f(xi+h, y+k2)
        y += (k0 + 2*k1 + 2*k2 + k3)/6
def int_definida(f,xx):
    int = 0
    for i in range(len(f)-1):
        int += (f[i] + f[i+1])/2 * (xx[i+1]-xx[i])
    return int
def RK_2DO_ORDEN(f0,f1,y00,y01,x0,xf,steps,omega):
    # h val
    h = (xf-x0)/(steps-1)
    xval = x0 + np.arange(steps)*h
    y0val = np.zeros(steps)
    y1val = np.zeros(steps)
    y0 = y00
    y1 = y01
    for j, xi in enumerate(xval):
        y0val[j] = y0
        k00 = h*f0(xi, y0,y1)
        k10 = h*f0(xi+h/2, y0+k00/2, y1+k00/2)
        k20 = h*f0(xi+h/2, y0+k10/2, y1 + k10/2)
        k30 = h*f0(xi+h, y0+k20,y1+k20)
        y0 += (k00 + 2*k10 + 2*k20 + k30)/6
        y1val[j] = y1
        k01 = h*f1(xi, y0,y1,omega)
        k11 = h*f1(xi+h/2, y0+k01/2, y1+k01/2,omega)
        k21 = h*f1(xi+h/2, y0+k11/2, y1 + k11/2,omega)
        k31 = h*f1(xi+h, y0+k21,y1+k21,omega)
        y1 +=(k01 + 2*k11 + 2*k21 + k31)/6
           
    return xval, y0val,y1val

def Int_evaluada(f,t):
    I = 0
    for i in range(len(f)-2):
        I += (t[i+1]-t[i])*(f[i+1])
    return I
def condicionesIniciales (l,ht,u,x,u_0t):
    N = len(u)
    w = np.zeros(N)
    for i in range (1 ,N -1):
            w[i] = (1 - l **2) * u[i] + 0.5 * l **2 * ( u [i+1] + u[i-1]) + ht * u_0t(x[i])
    return w
#SOLVER A ECUACIÓN DIFERENCIAL
def solver(u , w , N , x , Nt , l,ht,V,y0,yf,L ):
    s = np . zeros (N+2)
    for n in range (1 , Nt ):
        for i in range (1 , N+1):
            s [ i ] = 2 * (1 - l **2) * w [ i ] + l **2 * ( w [ i +1] + w [i -1]) - u[ i ] -l**2*V(w[i])
        s[0] = y0 #valor en +infinito y - infinito
        s[-1] = yf
        u = w.copy ()
        w = s.copy ()
        plt.plot(x ,s, 'ro--')
        plt.xlim([-L,L])
        plt.ylim([y0-1,yf+1])
        plt.grid()
        plt.savefig(f'frame{n}')
        plt.show()
    return s
def solver_akink(u , w , N , x , Nt , l,ht,V,y0,yf,L ):
    s = np . zeros (N+2)
    for n in range (1 , Nt ):
        for i in range (1 , N+1):
            s [ i ] = 2 * (1 - l **2) * w [ i ] + l **2 * ( w [ i +1] + w [i -1]) - u[ i ] -l**2*V(w[i])
        s[0] = y0 #valor en +infinito y - infinito
        s[-1] = yf
        u = w.copy ()
        w = s.copy ()
        plt.plot(x ,s, 'ro--')
        plt.xlim([-L,L])
        plt.ylim([yf-1,y0+1])
        plt.grid()
        plt.savefig(f'frame{n}')
        plt.show()
    return s
#Solver para energía de sombrero mexicano
def solEnergia_Kink(u , w , N , x , Nt , l,ht,V,y0,yf,L,h):
    s = np . zeros (N+2)
    for n in range (1 , Nt ):
        for i in range (1 , N+1):
            s [ i ] = 2 * (1 - l **2) * w [ i ] + l **2 * ( w [ i +1] + w [i -1]) - u[ i ] -l**2*V(w[i])
        s[0] = y0 #valor en +infinito y - infinito
        s[-1] = yf
        u = w.copy ()
        w = s.copy ()
        kinke = 1/2*(dif.central1(s,h))**2
        plt.plot(x ,kinke, 'ro--')
        plt.xlim([-L,L])
        plt.ylim([0,4])
        plt.grid()
        plt.savefig(f'frame{n}')
        plt.show()
    return s 