def NewtonSearch(f, uncertainty, xi):
    x = symbols('x')
    h = f.subs(x,xi).evalf()/diff(f).subs(x,xi).evalf()
    while (abs(h) > uncertainty):
        h = f.subs(x,xi).evalf()/diff(f).subs(x,xi).evalf()
        xi = xi - h
        ax3.plot3D(np.linspace(0, float(xi), 5), np.linspace(0, float(func.subs(x,xi).evalf()), 5), np.linspace(0,0,5), 'purple')
    return xi

