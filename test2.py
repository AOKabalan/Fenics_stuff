from dolfin import *

# Example setup (this part may vary based on your actual use case)
mesh = UnitSquareMesh(8, 8)
V = FunctionSpace(mesh, "Lagrange", 1)
u = Expression("sin(x[0] * x[1])", degree=1)
bc = DirichletBC(V, u, "on_boundary")

# Print some typical attributes of the DirichletBC object
print("DirichletBC object details:")
print(f"  Function space: {bc.function_space()}")
print(f"  Boundary condition: {bc.value()}")
