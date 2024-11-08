import numpy as np
import mdtraj as md
from sys import argv

# trajectory file
traj_file = md.formats.HDF5TrajectoryFile(argv[1])
data = traj_file.read()
potE = data.potentialEnergy
time = data.time
nsteps=len(time)
out_file=open('V.dat','w')
for i in range(nsteps):
	out_file.write(str(time[i])+' '+str(potE[i])+'\n')
print('<V> = ',np.mean(potE))
print('<V^2> - <V>^2= ',np.var(potE))
out_file.close()
traj_file.close()