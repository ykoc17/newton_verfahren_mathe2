import numpy as np
import matplotlib as plt

def newton(x0, f, Df, tol, itmax, x):

    return x_star


#Testbeispiele

def f1(x): #in-place ?
    x[0] = pow(x[0],2) + pow(x[1], 2) -1
    x[1] = x[0]*x[1]-0.25

def f2(x): #in-place
    x[0] = pow(x[0], 2)-x[1]-2
    x[1] = x[0]*x[1]+1

def f3(x): #in-place
    x[0] = x[0]*0.5*np.sin(np.pi*0.5*x[0])-x[1]
    x[1] = pow(x[1], 2)-x[0]+1

def Df1(x): #in-place
    x[0] =
    x[1] =

def Df2(x): #in-place
    x[0] =
    x[1] =

def Df3(x): #in-place
    x[0] =
    x[1] =



y1 = f3(np.array((1,1)))
print(y1)

