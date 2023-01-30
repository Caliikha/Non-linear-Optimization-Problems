from __future__ import annotations
import numpy as np
import math
import matplotlib
matplotlib.use('TKAgg')
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from sympy import *

class Vector3:
    def __init__(self, X, Y, Z):
        self.x = X
        self.y = Y
        self.z = Z

fig = plt.figure()
ax3 = plt.axes(projection='3d') # matplotlib projecction selection
ax3.set_proj_type('ortho')

V1 = Vector3(0, 0, 1) # z axis
V1xline = np.linspace(0,V1.x,10)
V1yline = np.linspace(0,V1.y,10)
V1zline = np.linspace(0,V1.z,10)

V2 = Vector3(0, 1, 0) # y axis
V2xline = np.linspace(0, V2.x, 10)
V2yline = np.linspace(0, V2.y, 10)
V2zline = np.linspace(0, V2.z, 10)

V3 = Vector3(1, 0, 0) # x axis
V3xline = np.linspace(0, V3.x, 10)
V3yline = np.linspace(0, V3.y, 10)
V3zline = np.linspace(0, V3.z, 10)

ax3.plot3D(V1xline, V1yline, V1zline, 'blue')
ax3.plot3D(V2xline, V2yline, V2zline, 'red')
ax3.plot3D(V3xline, V3yline, V3zline, 'orange')

domain_resolution = 100 # number of sampling points in the domain

def f(x): # single variable f : R -> R
    assert type(x) == float or type(x) == int or type(x) == np.float64
    return 8*math.exp(1 - x) + 7*math.log(x)

def domain(a = 0, b = 1, dom_res = 100):
    iteration_val = (b + a)/dom_res
    while a < b:
        yield a
        a += iteration_val

ax3.set_xlim3d(0,2) # setting the viewable domain area
ax3.set_ylim3d(0,10)
ax3.set_zlim3d(0,0)

a = 1
b = 2
v = Vector3(a, f(a), 0)
for i in domain(a, b, domain_resolution):
    ax3.plot3D(
            np.linspace(v.x,i,5), 
            np.linspace(v.y,f(i),5), 
            np.linspace(v.z,0,5),
            'black')
    v = Vector3(i, f(i), 0)


def GoldenSearch(f, dom, uncertainty):
    assert callable(f) and type(dom) == list and type(uncertainty) == float
    rho = (3 - np.sqrt(5))/2
    ax3.plot3D(np.linspace(dom[1], dom[1], 5), np.linspace(0, f(dom[1]), 5), np.linspace(0,0,5), 'blue')
    ax3.plot3D(np.linspace(dom[0], dom[0], 5), np.linspace(0, f(dom[0]), 5), np.linspace(0,0,5), 'purple')
    #print(f'dom: {dom},\n f({dom[0]}) = {f(dom[0])}, f({dom[1]}) = {f(dom[1])}')
    while (abs(dom[1] - dom[0]) > uncertainty):
        if f(dom[0]) < f(dom[1]):
            dom[1] = dom[1] - rho*(dom[1] - dom[0])
            ax3.plot3D(np.linspace(dom[1], dom[1], 5), np.linspace(0, f(dom[1]), 5), np.linspace(0,0,5), 'blue')
        else:
            dom[0] = dom[0] + rho*(dom[1] - dom[0])
            ax3.plot3D(np.linspace(dom[0], dom[0], 5), np.linspace(0, f(dom[0]), 5), np.linspace(0,0,5), 'purple')
        #print(f'dom: {dom},\n f({dom[0]}) = {f(dom[0])}, f({dom[1]}) = {f(dom[1])}')

    #print(f'dom: {dom}')
    return dom

def Fibonacci(n):
    assert type(n) == int
    return int((1/np.sqrt(5))*(((1 + np.sqrt(5))/2)**(n) - ((1 - np.sqrt(5))/2)**(n)))

def FibSearch(f, dom, uncertainty, N, epsilon=0.05):
    rho = 1 - Fibonacci(N)/Fibonacci(N+1)
    ax3.plot3D(np.linspace(dom[1], dom[1], 5), np.linspace(0, f(dom[1]), 5), np.linspace(0,0,5), 'blue')
    ax3.plot3D(np.linspace(dom[0], dom[0], 5), np.linspace(0, f(dom[0]), 5), np.linspace(0,0,5), 'purple')
    #print(f'dom: {dom},\n f({dom[0]}) = {f(dom[0])}, f({dom[1]}) = {f(dom[1])}, rho: {rho}')
    while (abs(dom[1] - dom[0]) > uncertainty):
        if f(dom[0]) < f(dom[1]):
            dom[1] = dom[1] - rho*(dom[1] - dom[0])
            ax3.plot3D(np.linspace(dom[1], dom[1], 5), np.linspace(0, f(dom[1]), 5), np.linspace(0,0,5), 'blue')
        else:
            dom[0] = dom[0] + rho*(dom[1] - dom[0])
            ax3.plot3D(np.linspace(dom[0], dom[0], 5), np.linspace(0, f(dom[0]), 5), np.linspace(0,0,5), 'purple')
        rho = 1 - rho/(1- rho)
        if (1/2 - epsilon) < rho <= (1/2 + epsilon): 
            rho = rho - epsilon
        #print(f'dom: {dom},\n f({dom[0]}) = {f(dom[0])}, f({dom[1]}) = {f(dom[1])}, rho: {rho}')
    #print(f'final range -> dom: {dom}')
    return dom

def NewtonSearch(f, uncertainty, xi):
    x = symbols('x')
    h = f.subs(x,xi).evalf()/diff(f).subs(x,xi).evalf()
    while (abs(h) > uncertainty):
        h = f.subs(x,xi).evalf()/diff(f).subs(x,xi).evalf()
        xi = xi - h
        ax3.plot3D(np.linspace(0, float(xi), 5), np.linspace(0, float(func.subs(x,xi).evalf()), 5), np.linspace(0,0,5), 'purple')
    return xi

def SecantSearch(f, uncertainty, x1, x0):
    x = symbols('x')
    while (f.subs(x,x1) != f.subs(x,x0)):
        xprime = ((x1 - ((diff(f).subs(x,x1).evalf()*(x1-x0)))))/(diff(f).subs(x,x1).evalf() - diff(f).subs(x,x0).evalf())
        x0 = x1
        x1 = xprime
    return x1

def SecantMethod(f, uncertainty, x1, x0):
    it = 0
    while (abs(f(x1)) > uncertainty):
        it += 1
        buffr = x1
        x1 = x1 - (x1 - x0)*f(x1) / (f(x1) - f(x0))
        x0 = buffr
    return x1

print(10*'-' + "Golden Search with Uncertainty = 0.3" + 10*'-')
print(GoldenSearch(f, [1,2], 0.3))
print(10*'-' + "Fibonacci Search with Uncertainty = 0.3, N = 4" + 10*'-')
print(FibSearch(f, [1,2], 0.3, 4))
print(10*'-' + "Golden Search with Uncertainty = 0.2" + 10*'-')
print(GoldenSearch(f, [1,2], 0.2))

print(10*'-' + "Newton-Raphson Search with Uncertainty=0.002, x0=0.07" + 10*'-')
x = symbols('x')
func = 48*x*(1+x)**4 - (1+x)**4 - 7
evaluation = NewtonSearch(func, 0.002, 0.07)
print(f'x: {evaluation}, f(x): {func.subs(x,evaluation)}')

print(10*'-' + "Secant Search with Uncertainty=10e-5, x0=0, x1=1" + 10*'-')
func = (x**2)*exp(-x/2) - 1
evaluation = SecantSearch(func, 10e-5, 1, 0)
print(f'x: {evaluation}, f(x): {func.subs(x,evaluation).evalf()}')

print(10*'-' + "Secant Method with Uncertainty=10e-5, x0=0, x1=1" + 10*'-')
def f(x):
    return (x**2)*math.exp(-x/2) - 1
evaluation = SecantMethod(f, 10e-5, 1, 0)
print(f'x: {evaluation}, f(x): {func.subs(x,evaluation).evalf()}')


plt.show()

