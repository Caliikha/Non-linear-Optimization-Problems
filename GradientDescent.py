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

