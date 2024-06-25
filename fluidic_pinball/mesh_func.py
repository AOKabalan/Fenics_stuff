#This file is to change the mesh from gmsh to be compatible with fenics
import meshio
import numpy

def create_mesh(mesh, cell_type, prune_z=False):
    cells = mesh.get_cells_type(cell_type)
    cell_data = mesh.get_cell_data("gmsh:physical", cell_type)
    points = mesh.points[:, :2] if prune_z else mesh.points
    out_mesh = meshio.Mesh(
        points=points, cells={cell_type: cells}, cell_data={"name_to_read": [cell_data]}
    )
    return out_mesh


# def create_mesh(mesh, cell_type, prune_z=False):
#     cells=mesh.get_cells_type(cell_type)
#     cell_data=mesh.get_cell_data("gmsh:physical", cell_type)
#     points = mesh.points[:,:2] if prune_z else mesh.points
#     out_mesh=meshio.Mesh(points=points, cells={
#                          cell_type: cells}, cell_data={"name_to_read": [cell_data]})
    
#     return out_mesh


# msh = meshio.read("rect.msh")
# triangle_mesh = create_mesh(msh, "triangle",True)
# line_mesh = create_mesh(msh, "line",True)
# meshio.write("mesh.xdmf", triangle_mesh)
# meshio.write("mf.xdmf", line_mesh) th
mesh_from_file = meshio.read("pinball.msh")

line_mesh = create_mesh(mesh_from_file, "line", prune_z=True)
meshio.write("mf.xdmf", line_mesh)

triangle_mesh = create_mesh(mesh_from_file, "triangle", prune_z=True)
meshio.write("mesh.xdmf", triangle_mesh)
