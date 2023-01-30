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

