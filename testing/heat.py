#approach 1

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
from dolfin import *

# Define the mesh and boundary markers
mesh = Mesh()
mvc = MeshValueCollection("size_t", mesh, mesh.topology().dim())
with XDMFFile("mesh.xdmf") as infile:
    infile.read(mesh)
    infile.read(mvc, "name_to_read")
cf = cpp.mesh.MeshFunctionSizet(mesh, mvc)

mvc = MeshValueCollection("size_t", mesh, mesh.topology().dim()-1)
with XDMFFile("mf.xdmf") as infile:
    infile.read(mvc, "name_to_read")
mf = cpp.mesh.MeshFunctionSizet(mesh, mvc)


cmap = colors.ListedColormap(['blue', 'yellow', 'green', 'red'])
idx = np.unique(mf.array())
print(idx)
labels = dict()

for i, id in enumerate(idx):
    labels[id] = []
    for cell in cells(mesh):
        for facet in facets(cell):
            
            if mf[facet] == id:
                
                for vertex in vertices(facet):
                    labels[id].append([vertex.point().x(), vertex.point().y()])


for i, id in enumerate(idx):
    boundary = np.array(labels[id])
    #
    
    c=id*np.ones_like(boundary[:, 0])
    
    sc = plt.scatter(boundary[:, 0], boundary[:, 1], c=id*np.ones_like(boundary[:, 0]), cmap=cmap, vmin=min(idx)-0.5, vmax=max(idx)+0.5)
plt.colorbar(sc, ticks=[0, 1, 2, 3])


# V = FunctionSpace(mesh, "Lagrange", 1)


# def boundary(x, on_boundary):
#     return on_boundary


# print(mesh.num_edges(), mesh.num_cells(), mesh.num_facets())

# # Define boundary condition

# bc = DirichletBC(FunctionSpace(mesh, "CG", 1), Constant(1), mf, 3)


# c=Expression(("0.0", "0.0"),degree=2)
# uout=DirichletBC(V, c, mf, 1)
# bc=[uout]