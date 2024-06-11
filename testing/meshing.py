from dolfin import *

# Read the mesh
mesh = Mesh()
with XDMFFile("rect.xdmf") as infile:
    infile.read(mesh)

# If you have physical groups
mvc = MeshValueCollection("size_t", mesh, 2)  # 2D mesh
with XDMFFile("rect.xdmf") as infile:
    infile.read(mvc, "name_to_read")

# Create the MeshFunction
mf = MeshFunction("size_t", mesh, mvc)

# Now you can use the mesh and mesh functions in your FEniCS simulation