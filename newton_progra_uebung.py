import numpy as np
import matplotlib.pyplot as plt

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
    return np.array([[2*x[0], 2*x[1]],
                     [x[1], x[0]]])

def Df2(x): #in-place ?
    return np.array([[2*x[0], -1],
                     [x[1], x[0]]])

def Df3(x): #in-place ?
    return np.array([[0.25*x[0]*np.pi*np.cos(np.pi*0.5*x[0])+0.5*np.sin(np.pi*0.5*x[0]), -1],
                     [-1, 2*x[1]]])

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


x_startwerte = np.linspace(-2.5, 2.5, num=10000)

x_star_werte = np.array([])
for x0 in x_startwerte:
    x_star = newton(x0, f4, Df4, 0.0001, 1000)
    x_star_werte = np.append(x_star_werte, x_star)

fig, (axa, axb) = plt.subplots(2,1)
fig.suptitle("Grenzwerte für jeden Startwert x0 aus [-2.5, 2.5]")
axa.plot(x_startwerte, x_star_werte, "_")
axb.plot(x_startwerte, f4(x_startwerte))
axb.plot(x_startwerte, np.zeros(x_startwerte.size))
"""
#c)

def f5(x, j=3):
    r = np.sqrt(pow(x[0], 2)+pow(x[1], 2))
    rj = pow(r, j)
    print(r)
    phi=0
    if x[0]>0 and x[1]>0: phi = np.arctan(x[1]/x[0])
    if x[0]<0 and x[1]>0: phi = np.pi - np.arctan(x[1]/x[0])
    if x[0]<0 and x[1]<0: phi = np.pi + np.arctan(x[1]/x[0])
    if x[0]>0 and x[1]<0: phi = np.pi*2 - np.arctan(x[1]/x[0])
    if x[0]==0:
        if x[1]>0: phi = np.pi/2
        if x[1]<0: phi = 3*np.pi/2
    if x[1]==0:
        if x[0]>0: phi = 0
        if x[0]<0: phi = np.pi

    print(phi)

    return np.array([rj*np.cos(j*phi),
                     rj*np.sin(j*phi)])

def Df5(x, j=3):
    r = np.sqrt(pow(x[0], 2)+pow(x[1], 2))
    phi=0
    if x[0]>0 and x[1]>0: phi = np.arctan(x[1]/x[0])
    if x[0]<0 and x[1]>0: phi = np.pi - np.arctan(x[1]/x[0])
    if x[0]<0 and x[1]<0: phi = np.pi + np.arctan(x[1]/x[0])
    if x[0]>0 and x[1]<0: phi = np.pi*2 - np.arctan(x[1]/x[0])
    if x[0]==0:
        if x[1]>0: phi = np.pi/2
        if x[1]<0: phi = 3*np.pi/2
    if x[1]==0:
        if x[0]>0: phi = 0
        if x[0]<0: phi = np.pi

    return np.array([[j*pow(r, j-1)*np.cos(j*phi),-1*pow(r, j)*j*np.sin(j*phi)],
                    [j*pow(r, j-1)*np.sin(j*phi), pow(r, j)*j*np.cos(j*phi)]])


x_ki = []
x_star = newton(np.array([0.1,0.1]), f5, Df5, 0.0001, 900, x_ki)
print(x_star)
print("")
for i in x_ki: print(i)

plt.show()
