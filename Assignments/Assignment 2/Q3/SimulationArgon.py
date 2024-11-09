from __future__ import print_function
from openmm import app
import openmm as mm
from openmm import unit
from sys import stdout
from openmmtools import testsystems
import mdtraj

# Simulation Parameters
######################
temperature = 300           # Temperature of the System (Kelvin)
steps = 1000                # Number of Steps/Frames to render in the Simulation (Unitless) (Integer)
skipSteps = 10              # Number of Steps to skip initially
equilSteps = 100            # Number of Steps 
N=200                       # Number of Argon-Argon Molecules
######################

def Simulation (temperature, steps, skipSteps, equilSteps, N, displayOutput = False):
    #test_sys = testsystems.LennardJonesFluid(nparticles=200, reduced_density=0.50)
    test_sys = testsystems.LennardJonesFluid(nparticles=N, mass=39.9*unit.dalton,sigma=3.4*unit.angstrom,epsilon=0.238*unit.kilocalories_per_mole, cutoff=None,reduced_density=0.9)
    (system, positions) = test_sys.system, test_sys.positions

    #integrator = mm.VerletIntegrator(.01 * unit.femtoseconds)
    integrator = mm.LangevinIntegrator(temperature*unit.kelvin, 1.0/unit.picoseconds,1.0*unit.femtoseconds)

    platform = mm.Platform.getPlatformByName('Reference')

    simulation = app.Simulation(test_sys.topology, system, integrator, platform)
    simulation.context.setPositions(test_sys.positions)

    if displayOutput:
        print('The size of the periodic box is: ', system.getDefaultPeriodicBoxVectors())
        print('Minimizing...')
        print('Initializing velocities to Boltzmann distribution')
        print('Equilibrating...')
    
    simulation.minimizeEnergy()
    simulation.context.setVelocitiesToTemperature(temperature*unit.kelvin)
    simulation.step(equilSteps*skipSteps)

    if displayOutput:
        simulation.reporters.append(app.StateDataReporter(stdout, skipSteps, step=True,
        potentialEnergy=True, temperature=True, progress=True, remainingTime=True,
        speed=True, totalSteps=steps, separator='\t'))
    else:
        log_file = open(f"simulation_log_argon_{N}.txt", "w")
        simulation.reporters.append(app.StateDataReporter(log_file, skipSteps, step=True,
        potentialEnergy=True, temperature=True, progress=True, remainingTime=True,
        speed=True, totalSteps=steps, separator='\t'))

    simulation.reporters.append(app.PDBReporter('ar_liquid_traj'+ str(N) + '.pdb', skipSteps))
    h5_reporter=mdtraj.reporters.HDF5Reporter('ar_liquid_traj' + str(N) + '.h5', skipSteps)
    simulation.reporters.append(h5_reporter)

    print('Simulation beginning...')
    simulation.step(steps*skipSteps)

    h5_reporter.close()