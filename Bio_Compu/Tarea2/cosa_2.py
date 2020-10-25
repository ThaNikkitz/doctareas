from simtk import openmm, unit
from openmmtools import testsystems
import time


L_side = 5.0 #nm
pme_cutoff = 2.5 #nm. <L_side/2 

N_steps = 100

# Try: nonbondedMethod: NoCutoff, Ewald, PME

water_box = testsystems.WaterBox(box_edge=L_side*unit.nanometer, nonbondedMethod = openmm.app.NoCutoff, nonbondedCutoff=pme_cutoff*unit.nanometers)
system = water_box.system

print("Water box with %i atoms"%len(water_box.positions))

integrator = openmm.LangevinIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.002*unit.picoseconds)
simulation = openmm.app.simulation.Simulation(water_box.topology, system, integrator)
simulation.context.setPositions(water_box.positions)

tic = time.time()
simulation.step(N_steps)
tac = time.time()

print ("Time = %1.3fs"%(tac-tic))
