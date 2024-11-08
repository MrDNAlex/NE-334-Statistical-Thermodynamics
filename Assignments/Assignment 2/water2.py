from __future__ import print_function

from openmm import app
import openmm as mm
from openmm import LocalEnergyMinimizer
from openmm import unit
import sys
import mdtraj
import mdtraj.reporters
import numpy as np

#####Parameters
steps = 1000
skipSteps = 1
Temperature=150. # temperature in Kelvin
dt = .1 * unit.femtoseconds 
print(dt)
ensemble='NVE' # NVE or NVT
#ensemble='NVT' # NVE or NVT
#####

#system = mm.System()

pdb = app.PDBFile("water2.pdb")
forcefield = app.ForceField('amber10.xml', 'tip3p.xml')
nonbonded = app.CutoffNonPeriodic
nonbonded = app.NoCutoff

#system = forcefield.createSystem(pdb.topology, nonbondedMethod=nonbonded, nonBondedCutoff=1e3*unit.nanometer, rigidWater=True)

system = forcefield.createSystem(pdb.topology, nonbondedMethod=nonbonded, 
                                     nonbondedCutoff=1e3*unit.nanometer,
                                     constraints=None,rigidWater=False)




if (ensemble == 'NVT'):
	integrator = mm.LangevinIntegrator(Temperature*unit.kelvin, 1.0/unit.picoseconds,dt)
if (ensemble == 'NVE'):
	integrator = mm.VerletIntegrator(dt)

#Use the next line for the Reference platform, slow, easier to read, will only use 1 core
platform = mm.Platform.getPlatformByName('Reference')
#Use the CPU platform, faster, can use multiple cores primarily to do non-bonded interactions (fft's in parallel, etc)
#platform = mm.Platform.getPlatformByName('CPU')
simulation = app.Simulation(pdb.topology, system, integrator, platform)
simulation.context.setPositions(pdb.positions)
simulation.context.computeVirtualSites()
#state = simulation.context.getState(getForces=True, getEnergy=True, getPositions=True)
#potential_energy = state.getPotentialEnergy()
#print(potential_energy)

#minimize the structure
LocalEnergyMinimizer.minimize(simulation.context, 1e-1)
#state = simulation.context.getState(getForces=True, getEnergy=True, getPositions=True)
#potential_energy = state.getPotentialEnergy()
#print(potential_energy)
	
simulation.context.setVelocitiesToTemperature(Temperature*unit.kelvin)

#Outputs progress to command line
simulation.reporters.append(app.StateDataReporter(sys.stdout, skipSteps, step=True, 
    potentialEnergy=True, temperature=True, progress=True, remainingTime=True, 
    speed=True, totalSteps=steps, separator='\t'))

#Saves trajectory to .pdb file that can be opened in VMD
simulation.reporters.append(app.PDBReporter('trajectory.pdb', skipSteps))
#Saves trajectory file to binary format
traj_filename='water2_'+ensemble+'.h5'
simulation.reporters.append(mdtraj.reporters.HDF5Reporter(traj_filename, skipSteps))

#Performs the simulation
simulation.step(steps)

#state = simulation.context.getState(getForces=True, getEnergy=True, getPositions=True)
#potential_energy = state.getPotentialEnergy()
#print(potential_energy)
#Close binary trajectory
simulation.reporters[2].close()
#Read the output file
output_file = mdtraj.formats.HDF5TrajectoryFile(traj_filename)
data = output_file.read()
output_file.close()
potE = data.potentialEnergy
kinE = data.kineticEnergy
positions = data.coordinates
totalE=kinE+potE
time = data.time
nsteps=len(time)
KE_output=open('KE.'+ensemble,'w')
PE_output=open('PE.'+ensemble,'w')
TE_output=open('TE.'+ensemble,'w')
rOO_output=open('rOO.'+ensemble,'w')
for i in range(nsteps):
	KE_output.write(str(time[i])+' '+str(kinE[i])+'\n')
	PE_output.write(str(time[i])+' '+str(potE[i]-potE[0])+'\n')
	TE_output.write(str(time[i])+' '+str(kinE[i]+potE[i]-potE[0])+'\n')
	rOO_output.write(str(time[i])+' '+str(np.linalg.norm(positions[i][0] - positions[i][3]))+'\n')
KE_output.close()
PE_output.close()
TE_output.close()
rOO_output.close()
