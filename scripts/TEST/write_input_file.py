import re
import os
import math

def write_input(Type, Ncycle1, Ncycle2, EpsilonP, Delta, Nsites):

    #RotDif = float(RotFac)*3.0

    f = open('path.inp','w')
    f.write('Simulation\n')
    f.write('sim_type 1\n')
    f.write('start_type %d\n' % Type)
    f.write('graphics 0\n')
    f.write('ncycle1 %d\n' % Ncycle1)
    f.write('ncycle2 %d\n' % Ncycle2)
    f.write('\n')
    f.write('Potential\n')
    f.write('sigma 1.0\n')
    f.write('mobilityT 1.0\n')
    f.write('mobilityR 3.0\n')
    f.write('epsilonP %lf\n' % EpsilonP)
    f.write('delta %lf\n' % Delta)
    f.write('nsites %lf\n' % Nsites)

    f.write('\n')
    f.write('System\n')
    f.write('npart 2\n')
    f.write('beta 1.0\n')
    f.write('boxl 4.0\n')
    f.write('\n')
    f.write('Path\n')
    f.write('timestep 0.0001\n')
    f.write('ninter 1\n')
    f.write('\n')
    f.write('TIS\n')
    f.write('nshoot 10\n')
    f.write('nrepswap 10\n')
    f.write('nstateswap 10\n')
    f.write('nreverse 10\n')
    f.write('stateswapbias 1\n')
    f.write('fixedbias 0\n')
    f.close()
    return













