import matplotlib.pyplot as plt
from dolfin import *
import meshio
mesh_from_file = meshio.read("mesh.msh")


def create_mesh(mesh, cell_type, prune_z=False):
    cells = mesh.get_cells_type(cell_type)
    cell_data = mesh.get_cell_data("gmsh:physical", cell_type)
    points = mesh.points[:, :2] if prune_z else mesh.points
    out_mesh = meshio.Mesh(points=points, cells={cell_type: cells}, cell_data={
                           "name_to_read": [cell_data]})
    return out_mesh


line_mesh = create_mesh(mesh_from_file, "line", prune_z=True)
meshio.write("bound.xdmf", line_mesh)

triangle_mesh = create_mesh(mesh_from_file, "triangle", prune_z=True)
meshio.write("lavender.xdmf", triangle_mesh)


mesh = Mesh()
with XDMFFile("lavender.xdmf") as infile:
    infile.read(mesh)
mvc = MeshValueCollection("size_t", mesh, 2)
with XDMFFile("bound.xdmf") as infile:
    infile.read(mvc, "name_to_read")
sub = cpp.mesh.MeshFunctionSizet(mesh, mvc)

# Strain function


def epsilon(u):
    return 0.5*(grad(u) + grad(u).T)
# return sym(grad(u))

# Stress function


def sigma(u):
    return lmbda_*div(u)*Identity(2) + 2*mu*epsilon(u)

# Define material properties


E = Constant(100000)
nu = Constant(0.3)
mu = E/(2*(1+nu))
lmbda_ = E*nu/((1+nu)*(1-2*nu))
flag_quad = 2
P2 = VectorElement("Lagrange", mesh.ufl_cell(), flag_quad)
TH = P2
W = FunctionSpace(mesh, TH)

# Define traction on the boundary and body forces
T = Constant((0.0, 0.0))
B = Constant((0.0, 0.0))

u = Function(W)
du = TrialFunction(W)
v = TestFunction(W)

# Boundary Conditions


def clamped_boundary(x, on_boundary):
    return on_boundary and x[1] < tol


bc1 = DirichletBC(W, Constant((0, 0)), clamped_boundary)
tol = 1E-5
top = CompiledSubDomain('on_boundary && near(x[1], 30, tol)', tol=1E-5)
applied_disp = 0.5
bc2 = DirichletBC(W.sub(1), Constant((applied_disp/2.0)), top)

bcs = [bc1, bc2]

a = inner(sigma(du), epsilon(v))*dx
l = dot(B, v)*dx + inner(T, v)*ds(1)
A_ass, L_ass = assemble_system(a, l, bcs)
solve(A_ass, u.vector(), L_ass)


plot(u, title="Displacement", mode="displacement")
plt.savefig("u.png")

file_results = XDMFFile("+-elasticity_reet.xdmf")
file_results.write(u, 0.)
file_results.close()