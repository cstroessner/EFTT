import numpy as np
from dolfin import *

def evalPDEfun(p,m):
    # input: vector of length 16 which specifies the values of the coefficients in the cookies
    # output: integral over the solution of the cookie pde 

    # compute the cookie positions
    #m = 2
    dd = pow(m,2)
    rho = 1/(4*m+2)
    centers = np.zeros((dd,2))
    for s in range(0,m):
        for t in range(0,m):
            centers[(s+1)+m*(t)-1,0] = 4*(s+1)-1
            centers[(s+1)+m*(t)-1,1] = 4*(t+1)-1
    centers = rho*centers

    # unit square mesh 
    d = 20
    mesh = UnitSquareMesh(d, d)
    #plot(mesh)

    # definition of function spaces
    V = FunctionSpace(mesh, 'Lagrange', 1)
    u = TrialFunction(V)
    v = TestFunction(V)
    usol = Function(V)

    # boundary conditions
    u0 = Expression('0',degree = 1)
    def u0_boundary(x, on_boundary):
        return on_boundary
    bc = DirichletBC(V, u0, u0_boundary)

    # custom expression evaluating coefficient for the cookie problem
    class MyExpression(UserExpression):
        def eval(self, val, x):
            for i in range(0,dd):
                y = centers[i][:]
                if np.linalg.norm(y-x) <= rho:
                    val[0] = p[i]
                    return
            val[0] = 1 #TODO set this value to the default value for the coefficient
        def __init__(self, **kwargs):
            super().__init__(**kwargs)

    # create my expression object
    my_expr = MyExpression(element=V._ufl_element, degree = 1)

    # variational formulation
    a = inner(my_expr*grad(u), grad(v))*dx
    f = Expression('1', degree = 1)
    L = f*v*dx

    # solving the problem
    solve(a == L, usol, bc)

    # compute the integral of the solution
    return assemble(usol*dx) 