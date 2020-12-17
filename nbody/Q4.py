# @Author: Killian Sheriff <Killian>
# @Date:   03-Dec-2020
# @Email:  killian.sheriff@gmail.com / killian.sheriff@mail.mcgill.ca
# @Last modified by:   Killian
# @Last modified time: 17-Dec-2020
import numpy as np
import particle as P
import nbody as nb
import pickle
from tools import save_frames, create_gif, plot_energy

npart = 100000
gridsize = 100
size = (gridsize, gridsize, gridsize)
init_v = 0
mass = 40
init_mass = [mass for t in range(npart)]

soft = 10
dt = 400
niter = 1000
freq = 10

s = P.system_init(npart, size, init_mass, specific_r=None, specific_v=init_v, soft=soft, power_law=True)
g = nb.nbody(size, s, dt, soft=soft)

posi = []
energy = []
for i in range(niter):
    p, e = g.evolve()
    posi.append(p)
    energy.append(e)


pickle.dump((posi, energy), open("Q4.p", "wb"))


save_frames("Q4.p", 'Q4', freq, gridsize)
create_gif('Q4', niter, freq, duration=niter / freq)
plot_energy("Q4.p")
