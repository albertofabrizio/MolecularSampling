import numpy as np
import ase.io as aio
import lebedev

#Define Parameters
# Radial distance in Angstrom
rad_min=3
rad_max=10
rad_step=20
# Unit Sphere with Lebedev Grid
nleb=50

# Read dimer
molecule=aio.read("dimer.xyz", format="xyz")

# Compute distance matrix
dist=molecule.get_all_distances()

#Generate Radial Distribution
radial=np.linspace(rad_min,rad_max,rad_step)

# Sample the unit Sphere with Lebedev Positions
grid=np.array(lebedev.Lebedev(nleb))
unit_sphere_xyz=grid[:,:3]

#Generate full grid
pos=[]
for i in range(rad_step):
    pos.append(radial[i]*unit_sphere_xyz)
pos=np.array(pos)

#Original Position
origin=molecule.positions[3,:]

# Modify the positions
for i in range(rad_step):
    for j in range(nleb):
        molecule.positions[3:6,:]+=-origin+pos[i,j,:]
        aio.write(filename="traj.xyz", images=molecule, format='xyz', append=True)
