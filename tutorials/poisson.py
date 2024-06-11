import fenics as fe 
import matplotlib.pyplot as plt
N_points = 12
Forcing_mag = 1.0


def main():
    mesh = fe.UnitSquareMesh(N_points,N_points)
    lagrange_polynomial_space = fe.FunctionSpace(mesh,'Lagrange',1)

    def boundary_boolean_fucntion(x,on_boundary):
        return on_boundary

    homo_dirichlet = fe.DirichletBC(lagrange_polynomial_space,
    fe.Constant(0.0),
    boundary_boolean_fucntion
    )

    u_trial = fe.TrialFunction(lagrange_polynomial_space)
    v_test =  fe.TestFunction(lagrange_polynomial_space)

    forcing =  fe.Constant(- Forcing_mag)
    weak_form_lhs = fe.dot(fe.grad(u_trial),fe.grad(v_test))*fe.dx

    weak_form_rhs = forcing*v_test*fe.dx

    #assemble and solve

     
    u_solution = fe.Function(lagrange_polynomial_space)

    fe.solve(
            weak_form_lhs == weak_form_rhs,
            u_solution,
            homo_dirichlet
        ) 
    
    c = fe.plot(u_solution, mode ='color')
    plt.colorbar(c)

    fe.plot(mesh)

    plt.show()


if __name__ == '__main__':
   main()
