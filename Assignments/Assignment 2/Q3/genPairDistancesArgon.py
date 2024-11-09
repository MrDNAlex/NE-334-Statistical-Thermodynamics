import numpy as np
import mdtraj as md
from sys import argv

# trajectory file
input_data = md.load(argv[1])

##############################
# parameters
nbins = 100
rmin = 0.1
# rmax=2.0594033385430914/2.
rmax = input_data.unitcell_lengths[0, 0] / 2.0
dr = (rmax - rmin) / float(nbins)
# volume=(rmax*2.)**3
volume = (
    input_data.unitcell_lengths[0, 0]
    * input_data.unitcell_lengths[0, 1]
    * input_data.unitcell_lengths[0, 2]
)

#############################

N = int(input_data.n_atoms)
Nsteps = input_data.n_frames
print("There are ", N, " Argon atoms in the trajectory")
pairs = []
for i in range(N):
    for j in range(i + 1, N):
        pairs.append([i, j])
print("There are ", len(pairs), " Ar-Ar pairs")

distances = md.compute_distances(input_data, pairs)

print("There are ", len(distances), " steps in the trajectory")

histo = np.zeros(nbins, float)
# accumulate histograms
n_count = 0
for ArAr in distances:
    for d in ArAr:
        index = int(np.floor((d - rmin) / dr))
        if index < nbins:
            histo[index] += 1.0

# normalize histogram and divide by jacobian
for i in range(nbins):
    r = rmin + i * dr
    histo[i] = histo[i] / (2.0 * np.pi * r * r * dr * N * N / volume) / float(Nsteps)
    
    
Ar_file = open(f"Ar_histo_{N}.txt", "w")
for i in range(nbins):
    r = rmin + i * dr
    Ar_file.write(str(rmin + i * dr) + " " + str(histo[i]) + "\n")


#
# Solution to Neighbors
#

ReinmannSum = 0
for i in range(nbins):
    r = rmin + i * dr
    if r < 0.5:
        ReinmannSum+= histo[i]*r**2*dr
        
Neighbors = (N/volume) * 4 * np.pi * ReinmannSum

print("Number of neighbours = ", Neighbors)

#
#
#
## Solution
#Ar_file = open(f"Ar_histo_{N}.txt", "w")
#NN = 0.0
#for i in range(nbins):
#    r = rmin + i * dr
#    Ar_file.write(str(rmin + i * dr) + " " + str(histo[i]) + "\n")
#    if r < 0.5:
#        NN += histo[i] * r * r
#        
#print("N = ", N)
#print("Volume =", volume)
#print("Delta r = ", r)
#print("Number of neighbours = ", dr * 4.0 * np.pi * (N / volume) * NN)
#Ar_file.close()
