import numpy as np
import matplotlib.pyplot as plt

def newton(x0, f, Df, tol, itmax, x=None):
    x_k = x0
    if x != None: x.append(x_k)

    it = 0
    while True:
        it += 1
        if x0.size == 1:
            delta_k = -(f(x_k)/Df(x_k))
        else:
            Df_xk = Df(x_k)
            if np.linalg.det(Df_xk) == 0:
                print("Jacobi-Matrix ist singulär. Kein Newton-Verfahren möglich!")
                break
            else:
                delta_k = np.linalg.solve(Df_xk, -f(x_k))

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


fig1, (axa1, axa2, axa3) = plt.subplots(3,1)
fig1.suptitle("a) (x_k, y_k), k=0,1,... -> gegen x*")
plt.xlabel("x_k")
plt.ylabel("y_k")


print("f1:")
x_ki = []
x_star = newton(np.array([3,1]), f1, Df1, 0.0001, 1000, x_ki)
print("x*=", x_star)
print("\nIterationsschritte:")

axa1.plot(x_star[0], x_star[1], '*', color="green")
it=0
for i in x_ki:
    print(it, ":", i)
    axa1.plot(i[0], i[1], '.', color="black")
    axa1.text(i[0], i[1], "x" + str(it), color="red", fontsize="small")
    it+=1

print("\nf2:")
x_ki = []
x_star = newton(np.array([7.345,0.1245]), f2, Df2, 0.0001, 1000, x_ki)
print("x*=", x_star)
print("\nIterationsschritte:")

axa2.plot(x_star[0], x_star[1], '*', color="green")
it=0
for i in x_ki:
    print(it, ":", i)
    axa2.plot(i[0], i[1], '.', color="black")
    axa2.text(i[0], i[1], "x" + str(it), color="red", fontsize="small")
    it+=1

print("\nf3:")
x_ki = []
x_star = newton(np.array([1,0.1]), f3, Df3, 0.0001, 1000, x_ki)
print("x*=", x_star)
print("\nIterationsschritte:")

axa3.plot(x_star[0], x_star[1], '*', color="green")
it=0
for i in x_ki:
    print(it, ":", i)
    axa3.plot(i[0], i[1], '.', color="black")
    axa3.text(i[0], i[1], "x" + str(it), color="red", fontsize="small")
    it+=1


#b)

def f4(x):
    return (pow(x,2)-4)*(pow(x,2)-1)

def Df4(x):
    return 4*pow(x,3)-10*x

print("\nf4:")
x_ki = []
x_star = newton(np.array([5]), f4, Df4, 0.0001, 1000, x_ki)
print("x*=", x_star)
print("\nIterationsschritte:")
it=0
for i in x_ki:
    print(it, ":", i)
    it+=1


x_startwerte = np.linspace(-2.5, 2.5, num=10000)
x_star_tol = np.float_power(10,-2)

fig, (axa, axb) = plt.subplots(2,1)
fig.suptitle("b) Grenzwerte für jeden Startwert x0 aus [-2.5, 2.5]")
plt.xlabel("x")
plt.ylabel("f4(x)")


for x in x_startwerte:
    x0 = np.array([x])
    x_star = newton(x0, f4, Df4, 0.0001, 1000)
    if np.abs(x_star+2) <=  x_star_tol:
        axa.plot(x0, x_star, '.', color="red")
        axb.plot(x0, f4(x0), '.', color="red")
    elif np.abs(x_star+1) <=  x_star_tol:
        axa.plot(x0, x_star, '.', color="limegreen")
        axb.plot(x0, f4(x0), '.', color="limegreen")
    elif np.abs(x_star-1) <=  x_star_tol:
        axa.plot(x0, x_star, '.', color="blue")
        axb.plot(x0, f4(x0), '.', color="blue")
    elif np.abs(x_star-2) <=  x_star_tol:
        axa.plot(x0, x_star, '.', color="yellow")
        axb.plot(x0, f4(x0), '.', color="yellow")

plt.figtext(0.27, 0.4, "red-> -2   green-> -1   blue-> 1   yellow-> 2")


#c)

k = 3        #k = j (j is reseved for imaginary number)
tol = pow(10, -6)

def f5(z):
    a = z[0]
    b = z[1]
    res = np.power((a+b*1j),k)-1
    return np.array([np.real(res), np.imag(res)])

def Df5(z):
    a = z[0]
    b = z[1]
    return np.matrix([[np.real(k*np.power((a+b*1j),k-1)), np.real(1j*k*np.power((a+b*1j),k-1))],    #re(d/da), re(d/db)
                      [np.imag(k*np.power((a+b*1j),k-1)), np.imag(1j*k*np.power((a+b*1j),k-1))]])   #im(d/da), im(d/db)

print("\nf5:")
x_ki = []
x_star = newton(np.array([-1,-1]), f5, Df5, tol, 1000, x_ki)
print("x*=", x_star)
print("\nIterationsschritte:")
it=0
for i in x_ki:
    print(it, ":", i)
    it+=1

a_werte = np.linspace(-1, 1, num=100)
b_werte = np.linspace(-1, 1, num=100)

k=3
fig3, (axa3) = plt.subplots(1,1)
fig3.suptitle("c) j=3, r: (1,0) g: (-1/2, √3/2) b: (-1/2, -√3/2)", fontsize="small")
plt.xlabel("a")
plt.ylabel("b")

for ia, a in enumerate(a_werte):
    for ib, b in enumerate(b_werte):
        x_star = newton(np.array([a,b]), f5, Df5, tol, 1000, x_ki)
        if (np.linalg.norm(np.array([1,0])-x_star)<1.e-8): axa3.plot(a, b, 'o',color="red")
        if (np.linalg.norm(np.array([-1/2,np.sqrt(3)/2])-x_star)<1.e-8): axa3.plot(a, b, 'o', color="green")
        if (np.linalg.norm(np.array([-1/2, np.sqrt(3)/-2]) - x_star) < 1.e-8): axa3.plot(a, b, 'o', color="blue")

axa3.plot(1, 0, 'x', color="yellow")
axa3.plot(-1/2, np.sqrt(3)/2, 'x', color="yellow")
axa3.plot(-1/2, np.sqrt(3)/-2, 'x', color="yellow")

k=4
fig4, (axa4) = plt.subplots(1,1)
fig4.suptitle("c) j=4, r: (1,0) g: (-1,0) b: (0,-1) y: (0,1)", fontsize="small")
plt.xlabel("a")
plt.ylabel("b")

for ia, a in enumerate(a_werte):
    for ib, b in enumerate(b_werte):
        x_star = newton(np.array([a,b]), f5, Df5, tol, 1000, x_ki)
        if (np.linalg.norm(np.array([1,0])-x_star)<1.e-8): axa4.plot(a, b, 'o',color="red")
        if (np.linalg.norm(np.array([-1,0])-x_star)<1.e-8): axa4.plot(a, b, 'o', color="green")
        if (np.linalg.norm(np.array([0,-1]) - x_star) < 1.e-8): axa4.plot(a, b, 'o', color="blue")
        if (np.linalg.norm(np.array([0, 1]) - x_star) < 1.e-8): axa4.plot(a, b, 'o', color="yellow")


axa4.plot(1, 0, 'x', color="black")
axa4.plot(-1, 0, 'x', color="black")
axa4.plot(0, -1, 'x', color="black")
axa4.plot(0, 1, 'x', color="black")

k=5
fig5, (axa5) = plt.subplots(1,1)
fig5.suptitle("c) j=5, r: (1,0) g: (0.25*(√5-1), (√0.5*(√5+5))/2) blu: (0.25*(-√5-1), (√0.5*(-√5+5))/2) y: (0.25*(-√5-1),(-√0.5*(-√5+5))/2) bla: (0.25*(√5-1),(-√0.5*(√5+5))/2)", fontsize="small")
plt.xlabel("a")
plt.ylabel("b")


for ia, a in enumerate(a_werte):
    for ib, b in enumerate(b_werte):
        x_star = newton(np.array([a,b]), f5, Df5, tol, 1000, x_ki)
        if (np.linalg.norm(np.array([1,0])-x_star)<1.e-8): axa5.plot(a, b, 'o',color="red")
        if (np.linalg.norm(np.array([0.25*(np.sqrt(5)-1), np.sqrt(0.5*(np.sqrt(5)+5))/2]) - x_star) < 1.e-8): axa5.plot(a, b, 'o', color="green")
        if (np.linalg.norm(np.array([0.25*(-np.sqrt(5)-1), np.sqrt(0.5*(-np.sqrt(5)+5))/2]) - x_star) < 1.e-8): axa5.plot(a, b, 'o', color="blue")
        if (np.linalg.norm(np.array([0.25*(-np.sqrt(5)-1), -np.sqrt(0.5*(-np.sqrt(5)+5))/2]) - x_star) < 1.e-8): axa5.plot(a, b, 'o', color="yellow")
        if (np.linalg.norm(np.array([0.25*(np.sqrt(5)-1), -np.sqrt(0.5*(np.sqrt(5)+5))/2]) - x_star) < 1.e-8): axa5.plot(a, b, 'o', color="black")

axa5.plot(1, 0, 'x', color="purple")
axa5.plot(0.25*(np.sqrt(5)-1), np.sqrt(0.5*(np.sqrt(5)+5))/2, 'x', color="purple")
axa5.plot(0.25*(-np.sqrt(5)-1), np.sqrt(0.5*(-np.sqrt(5)+5))/2, 'x', color="purple")
axa5.plot(0.25*(-np.sqrt(5)-1), -np.sqrt(0.5*(-np.sqrt(5)+5))/2, 'x', color="purple")
axa5.plot(0.25*(np.sqrt(5)-1), -np.sqrt(0.5*(np.sqrt(5)+5))/2, 'x', color="purple")

k=6
fig6, (axa6) = plt.subplots(1,1)
fig6.suptitle("c) j=6, r: (1,0) g: (1/2, √3/2) b: (-1/2, √3/2) y: (-1,0) bla: (-1/2, -√3/2) o: (1/2, -√3/2)", fontsize="small")
plt.xlabel("a")
plt.ylabel("b")

for ia, a in enumerate(a_werte):
    for ib, b in enumerate(b_werte):
        x_star = newton(np.array([a,b]), f5, Df5, tol, 1000, x_ki)
        if (np.linalg.norm(np.array([1,0])-x_star)<1.e-8): axa6.plot(a, b, 'o',color="red")
        if (np.linalg.norm(np.array([1/2,np.sqrt(3)/2])-x_star)<1.e-8): axa6.plot(a, b, 'o', color="green")
        if (np.linalg.norm(np.array([-1/2, np.sqrt(3)/2]) - x_star) < 1.e-8): axa6.plot(a, b, 'o', color="blue")
        if (np.linalg.norm(np.array([-1,0])-x_star)<1.e-8): axa6.plot(a, b, 'o',color="yellow")
        if (np.linalg.norm(np.array([-1/2,np.sqrt(3)/-2])-x_star)<1.e-8): axa6.plot(a, b, 'o', color="black")
        if (np.linalg.norm(np.array([1/2, np.sqrt(3)/-2]) - x_star) < 1.e-8): axa6.plot(a, b, 'o', color="orange")

axa6.plot(1, 0, 'x', color="purple")
axa6.plot(1/2, np.sqrt(3)/2, 'x', color="purple")
axa6.plot(-1/2, np.sqrt(3)/2, 'x', color="purple")
axa6.plot(-1, 0, 'x', color="purple")
axa6.plot(-1/2, np.sqrt(3)/-2, 'x', color="purple")
axa6.plot(1/2, np.sqrt(3)/-2, 'x', color="purple")

plt.show()
