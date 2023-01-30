from __future__ import annotations
import numpy as np
import math
import matplotlib
matplotlib.use('TKAgg')
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from sympy import *
from sympy.vector import gradient
from scipy import linalg

# define the variable's symbols
x = symbols('x')
y = symbols('y')
z = symbols('z')
w = symbols('w')

class Vector3:
    def __init__(self, X, Y, Z):
        self.x = X
        self.y = Y
        self.z = Z
    def __str__(self):
        return f"({self.x}, {self.y}, {self.z})"

class Vector4:
    def __init__(self, X, Y, Z, W):
        self.x, self.y, self.z, self.w = X, Y, Z, W

    def __str__(self):
        return f"({self.x}, {self.y}, {self.z}, {self.w})"
    
    def __add__(self, rhs):
        return Vector4(self.x + rhs.x,
                       self.y + rhs.y,
                       self.z + rhs.z,
                       self.w + rhs.w)

    def __sub__(self, rhs):
        return Vector4(self.x - rhs.x,
                       self.y - rhs.y,
                       self.z - rhs.z,
                       self.w - rhs.w)

    def magnitude(self):
        return np.sqrt(self.x**2 + self.y**2 + self.z**2 + self.w**2)

    def scale(self, k):
        return Vector4(self.x*k, self.y*k, self.z*k, self.w*k)

class Vectorfunc:
    def __init__(self, f, h = 0, g = 0, z = 0):
        if (h == g == h == 0):
            self.dimcod = 1
            self.expression = f
        else:
            self.dimcod = 4
            self.expression = [f, h, g, z]

    def gradient(self):
        return Vectorfunc(diff(self.expression, x), diff(self.expression, y), diff(self.expression, z), diff(self.expression, w))

    def eval(self, v):
        if self.dimcod == 1:
            return (self.expression.subs(x,v.x).subs(y,v.y).subs(z,v.z).subs(w,v.w)).evalf()
        elif self.dimcod == 4:
            return Vector4(self.expression[0].subs(x,v.x).subs(y,v.y).subs(z,v.z).subs(w,v.w).evalf(),
                           self.expression[1].subs(x,v.x).subs(y,v.y).subs(z,v.z).subs(w,v.w).evalf(),
                           self.expression[2].subs(x,v.x).subs(y,v.y).subs(z,v.z).subs(w,v.w).evalf(),
                           self.expression[3].subs(x,v.x).subs(y,v.y).subs(z,v.z).subs(w,v.w).evalf())
    
    def scale(self, k):
        return Vectorfunc(self.expression[0]*k,
                          self.expression[1]*k,
                          self.expression[2]*k,
                          self.expression[3]*k)

def gradient(f, coords, v = None):
    if x == None:
        return Matrix([ diff(f, c) for c in coords ])
    else:
        return Matrix([ diff(f, c).subs(x,v.x).subs(y,v.y).subs(z,v.z).subs(w,v.w) for c in coords ])


def GradientDescent(f, x0, uncertainty = 0.001):
    vecf = Vectorfunc(f)
    x1 = x0 - vecf.gradient().eval(x0).scale(0.05)
    print(abs(vecf.eval(x1) - vecf.eval(x0)) > uncertainty)
    while (abs(vecf.eval(x1) - vecf.eval(x0)) > uncertainty):
        x0 = Vector4(x1.x,x1.y,x1.z,x1.w)
        x1 = x0 - vecf.gradient().eval(x0).scale(0.05)
        print(x1)

    return x1

def GradDesc(f, x0, uncertainty = 0.001, alpha = 0.05):
    grad = Matrix([f]).jacobian(Matrix([x,y,z,w])).subs(x,x0.x).subs(y,x0.y).subs(z,x0.z).subs(w,x0.w)
    #x1 = Vector4(x0.x - grad[0]*0.05, x0.y - grad[1]*0.05, x0.z - grad[2]*0.05, x0.w - grad[3]*0.05)
    flagin = False
    while (flagin == False or abs(f.subs(x,x1.x).subs(y,x1.y).subs(z,x1.z).subs(w,x1.w).evalf() - f.subs(x,x0.x).subs(y,x0.y).subs(z,x0.z).subs(w,x0.w).evalf()) > uncertainty):
        if flagin == False:
            x1 = Vector4(x0.x - grad[0]*alpha, x0.y - grad[1]*alpha, x0.z - grad[2]*alpha, x0.w - grad[3]*alpha)
            flagin = True
        x0 = Vector4(x1.x,x1.y,x1.z,x1.w)
        x1 = Vector4(x0.x - grad[0]*alpha,
                     x0.y - grad[1]*alpha,
                     x0.z - grad[2]*alpha,
                     x0.w - grad[3]*alpha)
        grad = Matrix([f]).jacobian(Matrix([x,y,z,w])).subs(x,x1.x).subs(y,x1.y).subs(z,x1.z).subs(w,x1.w)
        print((Vector4(float(x1.x),float(x1.y),float(x1.z),float(x1.w)) - Vector4(float(x0.x),float(x0.y),float(x0.z),float(x0.w))).magnitude())
    return x1

def NewtonSearch(f, x0, uncertainty = 0.001):
    grad = Matrix([f]).jacobian(Matrix([x,y,z,w])).subs(x,x0.x).subs(y,x0.y).subs(z,x0.z).subs(w,x0.w)
    hess = simplify(hessian(f,[x,y,z,w])).subs(x,x0.x).subs(y,x0.y).subs(z,x0.z).subs(w,x0.w)
    decr_factor = hess.inv() * grad.transpose()
    x1 = Vector4(x0.x - decr_factor[0], x0.y - decr_factor[1], x0.z - decr_factor[2], x0.w - decr_factor[3])
    while(abs(f.subs(x,x1.x).subs(y,x1.y).subs(z,x1.z).subs(w,x1.w).evalf() - f.subs(x,x0.x).subs(y,x0.y).subs(z,x0.z).subs(w,x0.w).evalf()) > uncertainty):
        x0 = Vector4(x1.x, x1.y, x1.z, x1.w)
        grad = Matrix([f]).jacobian(Matrix([x,y,z,w])).subs(x,x0.x).subs(y,x0.y).subs(z,x0.z).subs(w,x0.w)
        hess = simplify(hessian(f,[x,y,z,w])).subs(x,x0.x).subs(y,x0.y).subs(z,x0.z).subs(w,x0.w)
        decr_factor = hess.inv() * grad.transpose()
        x1 = Vector4(x0.x - decr_factor[0], x0.y - decr_factor[1], x0.z - decr_factor[2], x0.w - decr_factor[3])
    return Vector4(x1.x.evalf(),x1.y.evalf(),x1.z.evalf(),x1.w.evalf())

def ConjGrad(f, x0, uncertainty = 0.001):
    hess = simplify(hessian(f,[x,y,z])).subs(x,x0.x).subs(y,x0.y).subs(z,x0.z)
    grad = Matrix([f]).jacobian(Matrix([x,y,z])).subs(x,x0.x).subs(y,x0.y).subs(z,x0.z)
    dk = -grad
    alpha_k = (-grad * dk.transpose())[0,0] / ((dk * hess) * dk.transpose())[0,0]
    x0 = Vector3(x0.x + alpha_k*dk[0], x0.y + alpha_k*dk[1], x0.z + alpha_k*dk[2])
    while(Vector4(float(grad[0]), float(grad[1]), float(grad[2]),0).magnitude() > uncertainty):
        grad = Matrix([f]).jacobian(Matrix([x,y,z])).subs(x,x0.x).subs(y,x0.y).subs(z,x0.z)
        hess = simplify(hessian(f,[x,y,z])).subs(x,x0.x).subs(y,x0.y).subs(z,x0.z)
        beta_k = (grad * hess * dk.transpose())[0,0] / (dk * hess * dk.transpose())[0,0]
        dk = -grad + beta_k*dk
        alpha_k = (-grad * dk.transpose())[0,0] / ((dk * hess) * dk.transpose())[0,0]
        x0 = Vector3(x0.x + alpha_k*dk[0], x0.y + alpha_k*dk[1], x0.z + alpha_k*dk[2])
    return x0


f = (x + 10*y)**2 + 5*(z - w)**2 + (y - 2*z)**4 + 10*(x - w)**4 - 2*x
v = Vector4(3,-1,0,1)

#print(GradDesc(f,Vector4(0.82,-0.07,0.27,0.46)))
print("Newton Search")
print(NewtonSearch(f, v))

f = 2*x**2 + 2*y**2 + 1.5*z**2 + x*z + 2*y*z - 3*x - 2*z
hess = simplify(hessian(f,[x,y,z,w])).subs(x,v.x).subs(y,v.y).subs(z,v.z).subs(w,v.w)
v = Vector3(0,0,0)

print("Conjugate Grad Search")
print(ConjGrad(f,v))


