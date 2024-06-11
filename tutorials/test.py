import fenics as fe
import matplotlib.pyplot as plt 

if __name__ == '__main__':
    n_elements = 32
    mesh = fe.UnitIntervalMesh(n_elements)

    # Function space
    lagrange_first_order = fe.FunctionSpace(mesh, "Lagrange", 1)

    # Boundary Value
    u_dirichlet = fe.Constant(0.0)

    # A function to return whether we are on the boundary 
    def boundary_boolean_function(x, on_boundary):
        return on_boundary

    # Homogeneous Dirichlet
    boundary_condition = fe.DirichletBC(lagrange_first_order, u_dirichlet, boundary_boolean_function)

    # Initial Condition
    initial_condition = fe.Expression("sin(3.141*x[0])", degree=1)

    # Discretize Initial Condition
    u_old = fe.interpolate(initial_condition, lagrange_first_order)
    plt.figure()
    fe.plot(u_old, label='t=0.0')

    # Time-stepping (IMPLICIT EULER)
    time_step = 0.1
    heat_source = fe.Constant(0.0)

    # FEM Problem
    u_trial = fe.TrialFunction(lagrange_first_order)
    v_test = fe.TestFunction(lagrange_first_order)

    # Define bilinear form and linear form
    a = u_trial * v_test * fe.dx + time_step * fe.dot(fe.grad(u_trial), fe.grad(v_test)) * fe.dx
    L = (u_old * v_test + time_step * heat_source * v_test) * fe.dx

    n_time_step = 5
    time_current = 0
     
    u_solution = fe.Function(lagrange_first_order)

    for i in range(n_time_step):
        time_current += time_step

        fe.solve(a == L, u_solution, boundary_condition)

        u_old.assign(u_solution)
        fe.plot(u_solution, label=f't = {time_current:.1f}')
    
    plt.legend()
    plt.title('HEAT')
    plt.xlabel('x-pos')
    plt.ylabel('Temp')
    plt.grid()  
    plt.show()
