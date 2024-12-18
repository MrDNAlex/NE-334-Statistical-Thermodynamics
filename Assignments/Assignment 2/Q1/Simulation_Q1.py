# Imports
from __future__ import print_function
from openmm import app
import openmm as mm
from openmm import LocalEnergyMinimizer
from openmm import unit
import sys
import mdtraj
import mdtraj.reporters
import numpy as np
import pandas as pd

#####Parameters
steps = 1000                        # Number of Simulation Steps/Frames Simulated/Rendered (Unitless) (Integer)
skipSteps = 1                       # The Step Size of each Frame Generated (Unitless) (Integer)
Temperature=225.                    # Temperature at which the simulation is simulated (Kelvin)
dt = 0.01                           # The Change in Time between each Frame (Femtoseconds)
nonbondedCutoff=1e3*unit.nanometer
ensemble='NVE'                      # NVE or NVT Ensemble Type (Unitless) (String)
#####


def Simulation (numSteps, temp, dtVal, ensembleType, nonBondCutoff, title=""):
    #system = mm.System()

    #####Parameters
    steps = numSteps
    skipSteps = 1
    Temperature=temp # temperature in Kelvin
    dt = dtVal * unit.femtoseconds
    ensemble=ensembleType # NVE or NVT
    #####

    pdb = app.PDBFile("water2.pdb")
    forcefield = app.ForceField('amber10.xml', 'tip3p.xml')
    nonbonded = app.CutoffNonPeriodic
    nonbonded = app.NoCutoff

    #system = forcefield.createSystem(pdb.topology, nonbondedMethod=nonbonded, nonBondedCutoff=1e3*unit.nanometer, rigidWater=True)

    system = forcefield.createSystem(pdb.topology, nonbondedMethod=nonbonded, 
                                        nonbondedCutoff=nonBondCutoff,
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
    log_file = open("simulation_log.txt", "w")
    simulation.reporters.append(app.StateDataReporter(log_file, skipSteps, step=True, 
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
    KE_output=open(f'KE{title}.'+ensemble,'w')
    PE_output=open(f'PE{title}.'+ensemble,'w')
    TE_output=open(f'TE{title}.'+ensemble,'w')
    rOO_output=open(f'rOO{title}.'+ensemble,'w')
    for i in range(nsteps):
        KE_output.write(str(time[i])+' '+str(kinE[i])+'\n')
        PE_output.write(str(time[i])+' '+str(potE[i]-potE[0])+'\n')
        TE_output.write(str(time[i])+' '+str(kinE[i]+potE[i]-potE[0])+'\n')
        rOO_output.write(str(time[i])+' '+str(np.linalg.norm(positions[i][0] - positions[i][3]))+'\n')
    KE_output.close()
    PE_output.close()
    TE_output.close()
    rOO_output.close()

# Generate the First 
#Simulation(steps, Temperature, dt, ensemble, nonbondedCutoff)

# Need to Graph the First Settings

#KE = pd.read_csv("KE.NVE", sep=' ', header=None, names=['Time', 'Kinetic Energy (kJ/mol)'])

#print(KE)




