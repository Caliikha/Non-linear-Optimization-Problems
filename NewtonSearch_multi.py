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

