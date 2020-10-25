from simtk import openmm, unit
from openmmtools import testsystems
import numpy
import time

file1 = open('./N-var2.dat', 'w+')

L_side = 6.0 #nm
pme_cutoff = 2.0 #nm, < L_side/2
file1.write('%1.3f \t %1.3f\n\n'%(L_side, pme_cutoff))

for Method in [openmm.app.NoCutoff ,openmm.app.Ewald, openmm.app.PME]:
    
    print(Method)
    file1.write('\t' + str(Method) + '\n\n')

    Ns = numpy.arange(150, 310, 10)
    file1.write('N\t\t time\n')
    for N in Ns:
        print(N)
        N_steps = int(N)

# Try: nonbondedMethod: NoCutoff, Ewald, PME

        water_box = testsystems.WaterBox(box_edge=L_side*unit.nanometer, nonbondedMethod = Method, nonbondedCutoff=pme_cutoff*unit.nanometers)
        system = water_box.system

        integrator = openmm.LangevinIntegrator(300*unit.kelvin, 1/unit.picosecond, 2*unit.femtoseconds)
        simulation = openmm.app.simulation.Simulation(water_box.topology, system, integrator)
        simulation.context.setPositions(water_box.positions)

        if N == 200:
            print("check!")

        tic = time.time()
        simulation.step(N_steps)
        tac = time.time()
        
        file1.write(str(N) + '\t\t' + "time = %1.3fs\n"%(tac-tic))

    file1.write('\n')

file1.write("Water box with %i atoms"%len(water_box.positions))

file1.close()
