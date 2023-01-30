def GoldenSearch(f, dom, uncertainty):
    assert callable(f) and type(dom) == list and type(uncertainty) == float
    rho = (3 - np.sqrt(5))/2
    ax3.plot3D(np.linspace(dom[1], dom[1], 5), np.linspace(0, f(dom[1]), 5), np.linspace(0,0,5), 'blue')
    ax3.plot3D(np.linspace(dom[0], dom[0], 5), np.linspace(0, f(dom[0]), 5), np.linspace(0,0,5), 'purple')
    while (abs(dom[1] - dom[0]) > uncertainty):
        if f(dom[0]) < f(dom[1]):
            dom[1] = dom[1] - rho*(dom[1] - dom[0])
            ax3.plot3D(np.linspace(dom[1], dom[1], 5), np.linspace(0, f(dom[1]), 5), np.linspace(0,0,5), 'blue')
        else:
            dom[0] = dom[0] + rho*(dom[1] - dom[0])
            ax3.plot3D(np.linspace(dom[0], dom[0], 5), np.linspace(0, f(dom[0]), 5), np.linspace(0,0,5), 'purple')

    return dom
