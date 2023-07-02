import numpy as np
import matplotlib as plt

def newton(x0, f, Df, tol, itmax, x=None):
    x_k = x0
    if x != None: x.append(x_k)

    it = 0
    while True:
        it += 1
        if x0.size == 1: delta_k = -(f(x_k)/Df(x_k))
        else: delta_k = np.linalg.solve(Df(x_k), -f(x_k))
        x_k_prev = x_k
        x_k = x_k + delta_k
        if x != None: x.append(x_k)
        if np.linalg.norm(x_k - x_k_prev) <= tol: break
        if it == itmax: break
    return x_k

#a)
#Testbeispiele

def f1(x): #in-place ?
    return np.array([pow(x[0],2) + pow(x[1], 2) -1, x[0]*x[1]-0.25])

def f2(x): #in-place ?
    return np.array([pow(x[0], 2)-x[1]-2, x[0]*x[1]+1])

def f3(x): #in-place
    return np.array([x[0]*0.5*np.sin(np.pi*0.5*x[0])-x[1], pow(x[1], 2)-x[0]+1])

def Df1(x): #in-place ?
    return np.array([[2*x[0], 2*x[1]], [x[1], x[0]]])

def Df2(x): #in-place ?
    return np.array([[2*x[0], -1],[x[1], x[0]]])

def Df3(x): #in-place ?
    return np.array([[0.25*x[0]*np.pi*np.cos(np.pi*0.5*x[0])+0.5*np.sin(np.pi*0.5*x[0]), -1], [-1, 2*x[1]]])

"""
x_ki = []
x_star = newton(np.array([0.1,0.1]), f2, Df2, 0.0001, 1000, x_ki)
print(x_star)
print("")
for i in x_ki: print(i)
"""

#b)

def f4(x):
    return (pow(x,2)-4)*(pow(x,2)-1)

def Df4(x):
    return 4*pow(x,3)-10*x

"""
x_ki = []
x_star = newton(np.array([5]), f4, Df4, 0.0001, 1000, x_ki)
print(x_star)
print("")
for i in x_ki: print(i)
"""

x_werte = np.linspace(-2.5, 2.5, num=1000)
print(x_werte)

x_star_werte = np.array([])
for x0 in x_werte:
    x_star = newton(x0, f4, Df4, 0.0001, 1000)
    x_star_werte = np.append(x_star_werte, x_star)

print(x_star_werte)

