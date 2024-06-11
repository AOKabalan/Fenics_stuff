import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from dolfin import *


from dolfin import *
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
ds_custom = Measure("ds", domain=mesh, subdomain_data=mf, subdomain_id=12)



# Extract boundaries
idx = np.unique(mvc.array())

labels = dict()

# Iterate over the labels, cells and facets
for i, id in enumerate(idx):
    labels[id] = []
    for cell in cells(mesh):
        for facet in facets(cell):
            if boundary_mesh_function[facet] == id:
                for vertex in vertices(facet):
                    labels[id].append([vertex.point().x(), vertex.point().y()])

# Plot boundaries
cmap = colors.ListedColormap(['blue', 'yellow', 'green', 'red'])
for i, id in enumerate(idx):
    boundary = np.array(labels[id])
    plt.scatter(boundary[:, 0], boundary[:, 1], c=id*np.ones_like(boundary[:, 0]), cmap=cmap, vmin=min(idx)-0.5, vmax=max(idx)+0.5)

plt.colorbar(plt.cm.ScalarMappable(cmap=cmap), ticks=idx)
plt.show()

