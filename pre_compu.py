import h5py
from dgtd import read_msh, write_mesh, write_field_points

##########################################################################################
##########################################################################################
r0 = 0.00246924
# ###Convert the computation mesh from msh to hdf5
nodes, elements = read_msh(f"computation.msh")
nodes = nodes/r0
write_mesh(nodes, elements, f"computation.hdf5")
