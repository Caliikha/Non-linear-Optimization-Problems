def Fibonacci(n):
    assert type(n) == int
    return int((1/np.sqrt(5))*(((1 + np.sqrt(5))/2)**(n) - ((1 - np.sqrt(5))/2)**(n)))

def FibSearch(f, dom, uncertainty, N, epsilon=0.05):
    rho = 1 - Fibonacci(N)/Fibonacci(N+1)
    ax3.plot3D(np.linspace(dom[1], dom[1], 5), np.linspace(0, f(dom[1]), 5), np.linspace(0,0,5), 'blue')
    ax3.plot3D(np.linspace(dom[0], dom[0], 5), np.linspace(0, f(dom[0]), 5), np.linspace(0,0,5), 'purple')
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
    return dom

