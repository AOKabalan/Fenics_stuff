import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
import dolfin
from dolfin import *
import wurlitzer
from wurlitzer import pipes



mesh = Mesh()
mvc = MeshValueCollection("size_t", mesh, mesh.topology().dim())
with XDMFFile("mesh.xdmf") as infile:
    infile.read(mesh)
    infile.read(mvc, "name_to_read")
cf = cpp.mesh.MeshFunctionSizet(mesh, mvc)

mvc = MeshValueCollection("size_t", mesh, mesh.topology().dim()-1)
with XDMFFile("mf.xdmf") as infile:
    infile.read(mvc, "name_to_read")
mf = cpp.mesh.MeshFunctionSizet(mesh, mvc)


cmap = colors.ListedColormap(['blue', 'yellow', 'green', 'red'])
idx = np.unique(mf.array())
#print(idx)
labels = dict()


# Define function spaces
# V_element = VectorFunctionSpace(mesh, 'P', 2)
# Q_element = FunctionSpace(mesh, 'P', 1)
# W_element = MixedElement([V_element, Q_element]) # Taylor-Hood
# W = FunctionSpace(mesh, W_element)

V_element = VectorElement("Lagrange", mesh.ufl_cell(), 2)
Q_element = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
W_element = MixedElement(V_element, Q_element) # Taylor-Hood
W = FunctionSpace(mesh, W_element)


# Set parameter values
dt = 0.01
T = 3
theta = 0.5 # Crank Nicolson


Re = 1.
u_bar = 1.
# u_in = Expression(("1.5*u_bar*4/(0.41*0.41)*x[1]*(0.41 - x[1])", "0."), u_bar=u_bar, degree=2)

u_in = Expression(("-u_bar/(4.5*4.5)*x[1]*x[1] + u_bar", "0."), u_bar=u_bar, degree=2)
nu = Constant(u_bar*0.1/Re) # obtained from the definition of Re = u_bar * diam / nu. In our case diam = 0.1.
f = Constant((0., 0.))
# Define inflow profile


# Define boundary conditions
bcu_inflow = DirichletBC(W.sub(0), u_in, mf, 1)
bcu_walls = DirichletBC(W.sub(0), Constant((0, 0)), mf, 3)
bcu_cylinder = DirichletBC(W.sub(0), Constant((0, 0)), mf,4)
bcp_outflow = DirichletBC(W.sub(1), Constant(0), mf, 2)
bcs = [bcu_inflow, bcu_walls, bcu_cylinder]
bcp = [bcp_outflow]


# vq     = TestFunction(W) # Test function in the mixed space
# up     = TrialFunction(W) # Trial function in the mixed space
# (u, p) = split(up) # Function in each subspace to write the functional
# (v, q) = split(vq) # Test function in each subspace to write the functional

# # print(f'The mesh data are: {info(mesh.data(), True)}')

# lhs = (   nu*inner(grad(u), grad(v))*dx
#         - div(v)*p*dx
#         + div(u)*q*dx
#       )
# rhs =     inner(f, v)*dx


# up = Function(W) # Function in the mixed space to store the solution
# solve(lhs == rhs, up, bcu)
# (u, p) = up.split()

# u_plot = plot(u, title="Velocity")
# plt.colorbar(u_plot)
# plt.show()





def solve_stokes(W, nu, bcs):
    """Solve steady Stokes and return the solution"""

    # Define variational forms
    u, p = TrialFunctions(W)
    v, q = TestFunctions(W)
    a = nu*inner(grad(u), grad(v))*dx - p*div(v)*dx - q*div(u)*dx
    L = inner(Constant((0, 0)), v)*dx

    # Solve the problem
    w = Function(W)
    solve(a == L, w, bcs)




    return w



def solve_navier_stokes(W, nu, bcs):
    """Solve steady Navier-Stokes and return the solution"""

    # Define variational forms
    v, q = TestFunctions(W)
    w = Function(W)
    u, p = split(w)
    F = nu*inner(grad(u), grad(v))*dx + inner(grad(u)*u, v)*dx \
        - p*div(v)*dx - q*div(u)*dx

    # Solve the problem
    solve(F == 0, w, bcs)

    return w




def solve_navier_stokes2(W, nu, bcs):
    """Solve steady Navier-Stokes and return the solution"""
    
    # Define variational forms
    v, q = TestFunctions(W)
    delta_up = TrialFunction(W) # Trial function in the mixed space 
    w = Function(W)
    u, p = split(w)
    # Residual
    F = (   nu*inner(grad(u), grad(v))*dx
        + inner(grad(u)*u, v)*dx
        - div(v)*p*dx
        + div(u)*q*dx
        - inner(f, v)*dx
        )
    # Jacobian
    J = derivative(F, w, delta_up)

    snes_solver_parameters = {"nonlinear_solver": "snes",
                          "snes_solver": {"linear_solver": "mumps",
                                          "maximum_iterations": 20,
                                          "report": True,
                                          "error_on_nonconvergence": True}}
    problem = NonlinearVariationalProblem(F, w, bcs, J)
    solver  = NonlinearVariationalSolver(problem)
    solver.parameters.update(snes_solver_parameters)
    with pipes() as (out, err):
        solver.solve()
    print(out.read())

    return w

def solve_unsteady_navier_stokes(W, nu, bcs, T, dt, theta):
    """Solver unsteady Navier-Stokes and write results
    to file"""

    # Current and old solution
    w = Function(W)
    u, p = split(w)

    w_old = Function(W)
    u_old, p_old = split(w_old)

    # Define variational forms
    v, q = TestFunctions(W)
    F = ( Constant(1/dt)*dot(u - u_old, v)
          + Constant(theta)*nu*inner(grad(u), grad(v))
          + Constant(theta)*dot(dot(grad(u), u), v)
          + Constant(1-theta)*nu*inner(grad(u), grad(v))
          + Constant(1-theta)*dot(dot(grad(u_old), u_old), v)
          - p*div(v)
          - q*div(u)
        )*dx
    J = derivative(F, w)

    # Create solver
    problem = NonlinearVariationalProblem(F, w, bcs, J)
    solver = NonlinearVariationalSolver(problem)
    solver.parameters['newton_solver']['linear_solver'] = 'mumps'

    f = XDMFFile('velocity_unteady_navier_stokes.xdmf')
    u, p = w.split()

    # Perform time-stepping
    t = 0
    while t < T:
        w_old.vector()[:] = w.vector()
        solver.solve()
        t += dt
        f.write(u, t)


def save_and_plot(w, name):
    """Saves and plots provided solution using the given
    name"""

    u, p = w.split()

    # Store to file
    with XDMFFile("results_{}/u.xdmf".format(name)) as f:
        f.write(u)
    with XDMFFile("results_{}/p.xdmf".format(name)) as f:
        f.write(p)

    # Plot
    plt.figure()
    pl = dolfin.plot(u, title='velocity {}'.format(name))
    plt.colorbar(pl)
    plt.figure()
    p_plot = dolfin.plot(p, title="Pressure")
    plt.colorbar(p_plot)
   
   



# # Solve Stokes
# w = solve_stokes(W, nu, bcs)
# save_and_plot(w, 'stokes')

w = solve_navier_stokes2(W, nu, bcs)
save_and_plot(w, 'navier-stokes')
plt.show()
