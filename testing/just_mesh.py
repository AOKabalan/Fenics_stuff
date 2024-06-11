from dolfin import *
import matplotlib.pyplot as plt
# Load mesh
mesh = Mesh()
with XDMFFile("mesh.xdmf") as infile:
    infile.read(mesh)

# Plot mesh
plot(mesh, title="Mesh with Boundaries")

# Show plot
plt.show()
