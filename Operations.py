def gradient(f, coords, v = None):
    if x == None:
        return Matrix([ diff(f, c) for c in coords ])
    else:
        return Matrix([ diff(f, c).subs(x,v.x).subs(y,v.y).subs(z,v.z).subs(w,v.w) for c in coords ])

