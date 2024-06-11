import fenics as fe 
import matplotlib.pyplot as plt

Cantilever_len = 1.0
Cantilever_width = 0.2
N_points_length = 10
N_points_width = 3
Lame_mu = 1.0
Lame_lambda = 1.25
Density = 1.0
Gravity_Acceleration = 0.016

def main(): 
    mesh = fe.BoxMesh(
        fe.Point(0.0,0.0,0.0),
        fe.Point(Cantilever_len,Cantilever_width, Cantilever_width),
        N_points_length,
        N_points_width,
        N_points_width
    )

    lagrange_vector_space = fe.VectorFunctionSpace(
        mesh,
        'Lagrange',
        1
    )
    
    def clamped_boundary(x, on_boundary):
        return on_boundary and x[0] < fe.DOLFIN_EPS 

    dirichlet_clamped_boundary = fe.DirichletBC(
        lagrange_vector_space,
        fe.Constant((0.0,0.0,0.0)),
        clamped_boundary
    )

    #strain and stress

    def epsilon(u):
        engineering_strain = 0.5 * (fe.nabla_grad(u) + fe.nabla_grad(u).T)
        return engineering_strain
    def sigma(u):
        cauchy_stress = (
            Lame_lambda*fe.tr(epsilon(u))*fe.Identity(3)
            +
            2*Lame_mu*epsilon(u)
        )
        return cauchy_stress
    

    u_trial = fe.TrialFunction(lagrange_vector_space)
    v_test =  fe.TestFunction(lagrange_vector_space)
    forcing = fe.Constant((0.0,0.0,- Density*Gravity_Acceleration))
    traction = fe.Constant((0.0,0.0,0.0))

    weak_form_lhs = fe.inner(sigma(u_trial),epsilon(v_test))*fe.dx

    weak_form_rhs =  (fe.dot(forcing,v_test)+ fe.dot(traction,v_test))*fe.dx

    u_solution = fe.Function(lagrange_vector_space)

    fe.solve(
            weak_form_lhs == weak_form_rhs,
            u_solution,
            dirichlet_clamped_boundary
        ) 

    #compute von-mises

    s = sigma(u_solution)- (1/3)*fe.tr(sigma(u_solution))*fe.Identity(3)

    von_mises_stress = fe.sqrt((3/2)*fe.inner(s,s))

    lagrange_scalar_space = fe.FunctionSpace(mesh,'Lagrange',1)

    von_mises_stress = fe.project(von_mises_stress,lagrange_scalar_space)


    #paraview
    u_solution.rename('Displacement Vector','')
    von_mises_stress.rename('Von Mises Stress','')

    beam_deflection_file = fe.XDMFFile('beam_deflection.xdmf')
    beam_deflection_file.parameters["flush_output"]= True
    beam_deflection_file.parameters["functions_share_mesh"]= True
    beam_deflection_file.write(u_solution,0.0)
    beam_deflection_file.write(von_mises_stress,0.0)

if __name__ == '__main__' :
    main()