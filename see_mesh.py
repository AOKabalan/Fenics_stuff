from fenics import *
import meshio
filename="flowpastcylinder"
msh=meshio.read(filename + ".msh")


def create_mesh(mesh, cell_type, prune_z=False):
     cells=mesh.get_cells_type(cell_type)
     cell_data=mesh.get_cell_data("gmsh:physical", cell_type)
     points = mesh.points[:,:2] if prune_z else mesh.points
     out_mesh=meshio.Mesh(points=points, cells={
     cell_type: cells}, cell_data={"name_to_read": [cell_data]})
    
     return out_mesh

# # tetra_mesh=create_mesh(msh, "tetra")
triangle_mesh=create_mesh(msh,"triangle")
# meshio.write(filename + "_tetra.xdmf", tetra_mesh)
meshio.write(filename + "_triangle.xdmf", triangle_mesh)


mesh = Mesh()
# with XDMFFile(filename + "_tetra.xdmf") as infile:
#     infile.read(mesh)
# mvc_3d = MeshValueCollection("size_t", mesh, 3) 
# with XDMFFile(filename + "_tetra.xdmf") as infile:
#     infile.read(mvc_3d, "name_to_read")
# mf_3d = cpp.mesh.MeshFunctionSizet(mesh, mvc_3d)

mvc_2d = MeshValueCollection("size_t", mesh, 2)
with XDMFFile(filename + "_triangle.xdmf") as infile:
    infile.read(mvc_2d, "name_to_read")
mf_2d = cpp.mesh.MeshFunctionSizet(mesh, mvc_2d)