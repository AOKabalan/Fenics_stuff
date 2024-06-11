import meshio

# Read the .msh file
msh = meshio.read("rect.msh")

# Write to XDMF format
meshio.write("rect.xdmf", msh)

