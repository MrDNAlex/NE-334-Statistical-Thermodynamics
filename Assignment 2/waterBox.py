import openmm as mm
from openmm import app
from openmm import unit
from openmmtools import testsystems
from sys import stdout
import mdtraj

# parameters
######################
temperature = 300.
steps = 1000
skipSteps = 10
equilSteps = 100
Box_edge=4.*unit.nanometers
#test_sys = testsystems.LennardJonesFluid(nparticles=250, reduced_density=1.)
# test_sys = testsystems.FlexibleWaterBox(box_edge=2.5*unit.nanometers, cutoff=1.*unit.nanometers)
#test_sys = testsystems.WaterBox(box_edge=2.5*unit.nanometers, cutoff=1.0*unit.nanometers, constrained=False)
test_sys = testsystems.WaterBox(box_edge=Box_edge, cutoff=Box_edge/2.)
(system, positions) = test_sys.system, test_sys.positions

print('The size of the periodic box is: ', system.getDefaultPeriodicBoxVectors())

integrator = mm.LangevinIntegrator(temperature*unit.kelvin, 1.0/unit.picoseconds,1.0*unit.femtoseconds)

platform = mm.Platform.getPlatformByName('Reference')
platform = mm.Platform.getPlatformByName('CPU')
platform = mm.Platform.getPlatformByName('OpenCL')
simulation = app.Simulation(test_sys.topology, system, integrator, platform)
simulation.context.setPositions(test_sys.positions)

print('Minimizing...')
simulation.minimizeEnergy()

print('Initializing velocities to Boltzmann distribution')
simulation.context.setVelocitiesToTemperature(temperature*unit.kelvin)

print('Equilibrating...')
simulation.step(equilSteps*skipSteps)

simulation.reporters.append(app.StateDataReporter(stdout, skipSteps, step=True,
    potentialEnergy=True, temperature=True, progress=True, remainingTime=True,
    speed=True, totalSteps=steps, separator='\t'))

simulation.reporters.append(app.PDBReporter('h2o_liquid_traj.pdb', skipSteps))
h5_reporter=mdtraj.reporters.HDF5Reporter('h2o_liquid_traj.h5', skipSteps)
simulation.reporters.append(h5_reporter)

print('Simulation beginning...')
simulation.step(steps*skipSteps)
h5_reporter.close()
