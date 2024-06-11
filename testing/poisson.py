# approach 2
import meshio

import matplotlib.pyplot as plt
import numpy as np
from dolfin import *

def create_mesh(mesh, cell_type, prune_z=False):
    cells = mesh.get_cells_type(cell_type)
    cell_data = mesh.get_cell_data("gmsh:physical", cell_type)
    points = mesh.points[:, :3] if prune_z else mesh.points
    out_mesh = meshio.Mesh(points=points, cells={cell_type: cells}, cell_data={"name_to_read": [cell_data.astype(np.int32)]})
    return out_mesh
# mesh file is helix1.msh


## cell type: tetrahedron and triangle is choosen in Salome, which step before generate .mesh file.
mesh = meshio.read("rect.msh")
meshio.write("helix.xdmf", create_mesh(mesh, "triangle", True))

surface_mesh = create_mesh(mesh, "triangle", prune_z=True)
line_mesh = create_mesh(mesh, "line", prune_z=True)
meshio.write("mesh.xdmf", surface_mesh)
meshio.write("mt.xdmf", line_mesh)

# with XDMFFile("mesh.xdmf", "r") as xdmf:
#     mesh = xdmf.read_mesh(name="Grid")
#     cd = xdmf.read_meshtags(mesh, name="Grid")    # volume mesh
# mesh.topology.create_connectivity(mesh.topology.dim, mesh.topology.dim - 1)
# with XDMFFile( "mt.xdmf", "r") as xdmf:
#     fd = xdmf.read_meshtags(mesh, name="Grid")    # facet element

mesh = Mesh()
mvc = MeshValueCollection("size_t", mesh, mesh.topology().dim())
with XDMFFile("mesh.xdmf") as infile:
    infile.read(mesh)
    infile.read(mvc, "name_to_read")
cd = infile.read_meshtags(mesh, name="Grid") 

print(mvc.values())
mvc = MeshValueCollection("size_t", mesh, mesh.topology().dim()-1)
with XDMFFile("mt.xdmf") as infile:
    infile.read(mvc, "name_to_read")
fd = infile.read_meshtags(mesh, name="Grid") 




# Define boundary and initial conditions

# bottom_cells = ct.find(bottom_narker),  
# bottom_marker is 2 or 3 whaich is defined in gmsh
# these information from tutorials: Defining subdomains for different materials

inlet_cells = fd.find('Inlet')
out_cells = fd.find('Boundary')

inlet_dofs = locate_dofs_topological(V, mesh.topology.dim - 1, inlet_cells)

bc1 = dirichletbc(default_scalar_type(20), inlet_dofs, V)


out_dofs = locate_dofs_topological(V, mesh.topology.dim - 1, out_cells)

bc2 = dirichletbc(default_scalar_type(-20), out_dofs, V)

bc = [bc1, bc2]