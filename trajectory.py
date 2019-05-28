import numpy as np
import ase.io as aio
import lebedev
import math

#Define Parameters
# Radial distance in Angstrom
rad_min=3
rad_max=10
rad_step=20
# Unit Sphere with Lebedev Grid
nleb=14
#Euler Rolling
nroll=6

# Read dimer
molecule=aio.read("dimer.xyz", format="xyz")

# Compute distance matrix
dist=molecule.get_all_distances()

#Generate Radial Distribution
radial=np.linspace(rad_min,rad_max,rad_step)

# Sample the unit Sphere with Lebedev Positions
grid=np.array(lebedev.Lebedev(nleb))
unit_sphere_xyz=grid[:,:3]

# Identify symmetry unique points
#C2V
mask = np.ones(len(unit_sphere_xyz), dtype=bool)
mask[1] = False
mask[2] = False
mask[5] = False
mask[7] = False
mask[10] = False
mask[11] = False
mask[12] = False
mask[13] = False
unit_sphere_xyz=unit_sphere_xyz[mask]

#Euler Angles and rotations
alpha=np.arccos(unit_sphere_xyz[:,2])
beta=np.arctan2(unit_sphere_xyz[:,1],unit_sphere_xyz[:,0])
gamma=np.linspace(0, np.pi,nroll)

def eulerAnglesToRotationMatrix(theta):
    R_x = np.array([[1,         0,                  0                   ],
                    [0,         math.cos(theta[0]), -math.sin(theta[0]) ],
                    [0,         math.sin(theta[0]), math.cos(theta[0])  ]
                    ])



    R_y = np.array([[math.cos(theta[1]),    0,      math.sin(theta[1])  ],
                    [0,                     1,      0                   ],
                    [-math.sin(theta[1]),   0,      math.cos(theta[1])  ]
                    ])

    R_z = np.array([[math.cos(theta[2]),    -math.sin(theta[2]),    0],
                    [math.sin(theta[2]),    math.cos(theta[2]),     0],
                    [0,                     0,                      1]
                    ])


    R = np.dot(R_z, np.dot( R_y, R_x ))

    return R

#Generate full grid
pos=[]
for i in range(rad_step):
    pos.append(radial[i]*unit_sphere_xyz)
pos=np.array(pos)

#Original Position
origin=molecule.positions[3,:]


# Modify the positions
for i in range(rad_step):
    for j in range(nleb-8):
        for eul1 in range(len(alpha)):
            for eul2 in range(len(gamma)):
                rotation=eulerAnglesToRotationMatrix(np.array([alpha[eul1],beta[eul1],gamma[eul2]]))
                molecule.positions[3,:]=np.dot(molecule.positions[3,:],rotation)
                molecule.positions[4,:]=np.dot(molecule.positions[4,:],rotation)
                molecule.positions[5,:]=np.dot(molecule.positions[5,:],rotation)
                molecule.positions[3:6,:]+=pos[i,j,:]-origin

                aio.write(filename="traj.xyz", images=molecule, format='xyz', append=True)
