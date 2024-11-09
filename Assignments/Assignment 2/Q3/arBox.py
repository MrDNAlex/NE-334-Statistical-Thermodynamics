from __future__ import print_function
from openmm import app
import openmm as mm
from openmm import unit
from sys import stdout
from openmmtools import testsystems
import mdtraj

# parameters
######################
temperature = 300.
steps = 1000
skipSteps = 10
equilSteps = 100
N=200
#test_sys = testsystems.LennardJonesFluid(nparticles=200, reduced_density=0.50)
test_sys = testsystems.LennardJonesFluid(nparticles=N, mass=39.9*unit.dalton,sigma=3.4*unit.angstrom,epsilon=0.238*unit.kilocalories_per_mole, cutoff=None,reduced_density=0.9)
(system, positions) = test_sys.system, test_sys.positions

#integrator = mm.VerletIntegrator(.01 * unit.femtoseconds)
integrator = mm.LangevinIntegrator(temperature*unit.kelvin, 1.0/unit.picoseconds,1.0*unit.femtoseconds)

platform = mm.Platform.getPlatformByName('Reference')

simulation = app.Simulation(test_sys.topology, system, integrator, platform)
simulation.context.setPositions(test_sys.positions)

print('The size of the periodic box is: ', system.getDefaultPeriodicBoxVectors())

print('Minimizing...')
simulation.minimizeEnergy()

print('Initializing velocities to Boltzmann distribution')
simulation.context.setVelocitiesToTemperature(temperature*unit.kelvin)

print('Equilibrating...')
simulation.step(equilSteps*skipSteps)

simulation.reporters.append(app.StateDataReporter(stdout, skipSteps, step=True,
    potentialEnergy=True, temperature=True, progress=True, remainingTime=True,
    speed=True, totalSteps=steps, separator='\t'))

simulation.reporters.append(app.PDBReporter('ar_liquid_traj'+ str(N) + '.pdb', skipSteps))
h5_reporter=mdtraj.reporters.HDF5Reporter('ar_liquid_traj' + str(N) + '.h5', skipSteps)
simulation.reporters.append(h5_reporter)

print('Simulation beginning...')
simulation.step(steps*skipSteps)

h5_reporter.close()