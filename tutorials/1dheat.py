import fenics as fe
import matplotlib.pyplot as plt 

if __name__ == '__main__':
    n_elements = 32
    mesh = fe.UnitIntervalMesh(n_elements)


    #Function_space

    lagrange_first_order = fe.FunctionSpace(
        mesh,
        "Lagrange",
        1
    )

    #Boundary Value
    u_dirichlet = fe.Constant(0.0)

    #A function to return whether we are on the boundary 

    def boundary_boolean_fucntion(x,on_boundary):
        return on_boundary

    #Homogeneous Dirichlet
    boundary_condition = fe.DirichletBC(
        lagrange_first_order,
        u_dirichlet,
        boundary_boolean_fucntion
    )

    #Initial Condition
    initial_condition = fe.Expression(
        "sin(3.141*x[0])",
        degree = 1
    )

    #Discrertize Initial Condition
    u_old = fe.interpolate(
        initial_condition,
        lagrange_first_order
    )
    plt.figure()
    fe.plot(u_old, label = 't=0.0')
    #Time_stepping (IMPLICIT EULER)
    time_step = 0.1
    heat_source = fe.Constant(0.0)

    #FEM Problem
    u_trial = fe.TrialFunction(lagrange_first_order)
    v_test = fe.TestFunction(lagrange_first_order)


    weak_form_residual = (
        u_trial*v_test*fe.dx
        +
        time_step*fe.dot(
            fe.grad(u_trial),
            fe.grad(v_test)
        )*fe.dx
        -
        (
            u_old*v_test*fe.dx
            +
            time_step*heat_source*v_test*fe.dx
        )
    )

    weak_form_lhs = fe.lhs(weak_form_residual)
    weak_form_rhs = fe.rhs(weak_form_residual)

    n_time_Step = 5
    time_current = 0
     
    u_solution = fe.Function(lagrange_first_order)

    for i in range (n_time_Step):
        time_current += time_step

        fe.solve(
            weak_form_lhs == weak_form_rhs,
            u_solution,
            boundary_condition
        ) 

        u_old.assign(u_solution)
        fe.plot(u_solution, label = f't = {time_current:1.1f}')
    
    plt.legend()
    plt.title('HEAT')
    plt.xlabel('x-pos')
    plt.ylabel('Temp')
    plt.grid()  
    plt.show()
