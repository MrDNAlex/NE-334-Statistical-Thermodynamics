from __future__ import print_function

# Imports
import numpy as np
import pandas as pd
from scipy.integrate import solve_bvp, solve_ivp
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

from openmm import app
import openmm as mm
from openmm import LocalEnergyMinimizer
from openmm import unit
import sys
import mdtraj
import mdtraj.reporters
import numpy as np

# Redo Q1 and Q2, try and learn Pandas to do this, increase the Step Size and Dt independent from each other to see their affects
# Honestly a hard reset is probably what needs to be done. Make a single Graph, only for the Real Settings I want to test for have that as it's own system
# Then for all the increase in stuff have it's own section / function to add it to a Data Frame meant for it's own task

# You will need to make sure the pdb file in is the same folder as the file,
# and you also need to create a Avg and Graphs folder for this to work

def Simulation (Numsteps, temp, dtConst, ensembleType):
    #####Parameters - Alexandre Dufresne-Nappert - 20948586
    steps = Numsteps
    skipSteps = 1
    Temperature=temp # temperature in Kelvin
    dt = dtConst * unit.femtoseconds 
    print(dt)
    ensemble=ensembleType


    # NVE or NVT
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


def GenGraphs (index, Numsteps, temp, dtConst, ensembleType, saveLastOnly = False):
    
    # Run the Simulation
    Simulation(Numsteps, temp, dtConst, ensembleType)
    
    # Define File Names and Paths
    fileName = ["Kinetic Energy", "Potential Energy", "Total Energy"]
    files = [f"KE.{ensembleType}", f"PE.{ensembleType}", f"TE.{ensembleType}"]

    Globaltime = 0
    Globalenergies = []
    
    avgLines = []

    # Loop through all File Types
    for i in range(len(fileName)):
        
        # Create List to store Lines
        lines = []
        
        # Extract all lines from file
        with open(files[i], "r") as file:
            lines = file.readlines()
        
        # Create Lists to store values
        time:list[float] = []
        energy:list[float] = []   
        
        time.clear()
        energy.clear()
        
        # Loop through every line and extract the Value for Time and 
        for line in lines:
            # Split the line at the space
            nums = line.split()
            
            # Cast the value and add to the list
            time.append(float(nums[0]))
            energy.append(float(nums[1]))
        
        # Get Average and Variance
        average = np.average(energy)
        variance = np.var(energy)
        
        # Append Lines for File
        avgLines.append(f"{fileName[i]} (Settings : Steps = {Numsteps}, Dt = {dtConst} * femtosec, Ensemble = {ensembleType})")
        avgLines.append("\n")
        avgLines.append(f"Average : {average} (kJ/mole)")
        avgLines.append("\n")
        avgLines.append(f"Variance : {variance} (kJ/mole)")
        avgLines.append("\n")
        avgLines.append("\n")
        
        # Store the Time and add the energy
        Globaltime = time
        Globalenergies.append(energy)

        # Plot individual Graphs
        plt.figure(figsize=(16, 10))
        plt.ylabel(f"{fileName[i]} of Molecules (kJ/mole)")
        plt.xlabel(f"Simulation Time (femtoseconds)")
        plt.plot(time, energy, "--")
        plt.title(f"Measured {fileName[i]} in a Simulation between 2 Water Molecules\n(Settings : Steps = {Numsteps}, Temp = {temp}, Dt = {dtConst} * femtosec, Ensemble = {ensembleType})")
        
        if not saveLastOnly:
            plt.savefig(f"Graphs/{fileName[i]}_Q1_{ensembleType}_{index}.png")
        #plt.show()
    
    # Save Avgs and Vars
    with open(f"Avg/AvgAndVar_{index}_{ensembleType}.txt", "w") as file:
        file.writelines(avgLines)

    # Plot the Comparison of all 3
    plt.figure(figsize=(16, 10))
    plt.ylabel(f"Energy (kJ/mole)")
    plt.xlabel(f"Simulation Time (femtoseconds)")
    for i in range(len(Globalenergies)):
        plt.plot(Globaltime, Globalenergies[i], "--", label=f"{fileName[i]}")

    plt.legend()
    plt.title(f"Comparison of Kinetic Energy, Potential Energy and Total Energy of a System of 2 Water Molecules over the course of a Simulation \n(Settings : Steps = {Numsteps}, Temp = {temp}, Dt = {dtConst} * femtosec, Ensemble = {ensembleType})")
    plt.savefig(f"Graphs/Comparison_Q1_{ensembleType}_{index}.png")
    #plt.show() 
    
    # Create List to store Lines
    lines = []
    
    # Extract all lines from file
    with open(f"rOO.{ensembleType}", "r") as file:
        lines = file.readlines()
    
    "Molecule Distance"
    distance:list[float] = []
    time:list[float] = []
    
    for line in lines:
        # Split the line at the space
        nums = line.split()
        
        # Cast the value and add to the list
        time.append(float(nums[0]))
        distance.append(float(nums[1]))
    
    # Get Average and Variance
    average = np.average(distance)
    variance = np.var(distance)
    
    print(f"\nMolecule Distance (Settings : Steps = {Numsteps}, Dt = {dtConst} * femtosec, Ensemble = {ensembleType})")
    print(f"Average : {average} (Angstrom)")
    print(f"Variance : {variance} (Angstrom)")
    
    
    plt.figure(figsize=(16, 10))
    plt.ylabel(f"Distance Between Molecules (Angstrom)")
    plt.xlabel(f"Simulation Time (femtoseconds)")
    plt.plot(time, distance)
    plt.title(f"Distance between 2 Water Molecules over the course of a Simulation\n(Settings : Steps = {Numsteps}, Temp = {temp}, Dt = {dtConst} * femtosec, Ensemble = {ensembleType})")
    plt.savefig(f"Graphs/MolDistance_{index}_{ensembleType}.png")
    
ensemble = "NVE"    

# Generate Main Settings
GenGraphs(0, 10000, 298, 0.01, ensemble)

# Get the Averages And Variances for new Settings
Steps = 1000
index = 1
for i in range(2):
    DT = 0.001
    for j in range(3):
        GenGraphs(index, Steps, 298, DT, ensemble, True)
        index+=1
        DT *=2
    Steps*=2
    
    