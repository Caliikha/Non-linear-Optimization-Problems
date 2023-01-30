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
