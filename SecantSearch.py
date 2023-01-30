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

