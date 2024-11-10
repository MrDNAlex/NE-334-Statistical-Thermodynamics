import numpy as np
import mdtraj as md
from sys import argv

# trajectory file
input_data = md.load(argv[1])
boxLength = argv[2]

##############################
# parameters
nbins=100
rmin=0.1
rmax=input_data.unitcell_lengths[0,0]/2.
dr=(rmax-rmin)/float(nbins)
volume=input_data.unitcell_lengths[0,0]*input_data.unitcell_lengths[0,1]*input_data.unitcell_lengths[0,2]

#############################

N = int(input_data.n_atoms/3)
Nsteps=input_data.n_frames


print(f"Number of Particles = {N}")
print(boxLength)
OO_pairs = []
OH_pairs = []
for i in range(N):
    for j in range(i+1,N):
        OO_pairs.append([i*3,j*3])

        OH_pairs.append([i*3,j*3+1])
        OH_pairs.append([i*3,j*3+2])
        OH_pairs.append([j*3,i*3+1])
        OH_pairs.append([j*3,i*3+2])

print('There are ', len(OO_pairs), ' O-O pairs')
print('There are ', len(OH_pairs), ' O-H pairs')

OO_distances = md.compute_distances(input_data,OO_pairs)
OH_distances = md.compute_distances(input_data,OH_pairs)

print('There are ', len(OO_distances), ' steps in the trajectory')
print('There are ', len(OH_distances), ' steps in the trajectory')

OO_histo=np.zeros(nbins,float)
#accumulate histograms
for OO in OO_distances:
	for d in OO:
		index_OO=int(np.floor((d-rmin)/dr))
		if index_OO < nbins:
			OO_histo[index_OO]+=1.
OH_histo=np.zeros(nbins,float)
#accumulate histograms
for OH in OH_distances:
	for d in OH:
		index_OH=int(np.floor((d-rmin)/dr))
		if index_OH < nbins:
			OH_histo[index_OH]+=1.

#normalize histogram and divide by jacobian
for i in range(nbins):
        r=rmin+i*dr
        OO_histo[i]=OO_histo[i]/(2.*np.pi*r*r*dr*N*N/volume)/float(Nsteps)
        OH_histo[i]=OH_histo[i]//(2.*np.pi*r*r*dr*N*N*4./volume)/float(Nsteps)

# Save the g(r) function
OO_file=open(f'OO_histo_{Nsteps}_{boxLength}.txt','w')
for i in range(nbins):
	OO_file.write(str(rmin+i*dr)+' '+str(OO_histo[i])+'\n')
OO_file.close()

# Save the g(r) function
OH_file=open(f'OH_histo_{Nsteps}_{boxLength}.txt','w')
for i in range(nbins):
	OH_file.write(str(rmin+i*dr)+' '+str(OH_histo[i])+'\n')
OH_file.close()

#
# Solution to the Neighbors
#

ReinmannSumOO = 0
for i in range(nbins):
    r = rmin + i * dr
    if r < 0.33: # May need to change these
        ReinmannSumOO+= OO_histo[i]*r**2*dr
        
NeighborsOO = (N/volume) * 4 * np.pi * ReinmannSumOO

print("Number of neighbours O-O = ", NeighborsOO)

ReinmannSumOH = 0
for i in range(nbins):
    r = rmin + i * dr
    if r < 0.25: # May need to change these
        ReinmannSumOH+= OH_histo[i]*r**2*dr
        
NeighborsOH = 2*(N/volume) * 4 * np.pi * ReinmannSumOH

print("Number of neighbours O-H = ", NeighborsOH)
