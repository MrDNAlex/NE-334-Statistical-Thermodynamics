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


def Simulation (temperature, steps, skipSteps, equilSteps, Box_edge, displayOutput = False):
    
    temperature = temperature
    steps = steps
    skipSteps = skipSteps
    equilSteps = equilSteps
    Box_edge= Box_edge *unit.nanometers
    
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

    if displayOutput:
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
        log_file = open(f"simulation_log_waterBox_{steps}.txt", "w")
        simulation.reporters.append(app.StateDataReporter(log_file, skipSteps, step=True,
            potentialEnergy=True, temperature=True, progress=True, remainingTime=True,
            speed=True, totalSteps=steps, separator='\t'))

    simulation.reporters.append(app.PDBReporter(f'h2o_liquid_traj_{steps}_{str(Box_edge).strip(" nm")}.pdb', skipSteps))
    h5_reporter=mdtraj.reporters.HDF5Reporter(f'h2o_liquid_traj_{steps}_{str(Box_edge).strip(" nm")}.h5', skipSteps)
    simulation.reporters.append(h5_reporter)

    print('Simulation beginning...')
    simulation.step(steps*skipSteps)
    h5_reporter.close()
