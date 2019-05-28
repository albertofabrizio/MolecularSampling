import numpy as np
import ase.io as aio

#Define Parameters
# Radial distance in Angstrom
rad_min=3
rad_max=10
rad_step=20

# Read dimer
molecule=aio.read("dimer.xyz", format="xyz")

# Compute distance matrix
dist=molecule.get_all_distances()

#Generate Radial Distribution

radial=np.linspace(rad_min,rad_max,rad_step)

# Modify the distance matrix
for i in range(rad_step):
    molecule.set_distance(0,3,radial[i], fix=0,mask=[0,0,0,1,1,1])
    aio.write(filename="traj.xyz", images=molecule, format='xyz', append=True)
