# @Author: Killian Sheriff <Killian>
# @Date:   03-Dec-2020
# @Email:  killian.sheriff@gmail.com / killian.sheriff@mail.mcgill.ca
# @Last modified by:   Killian
# @Last modified time: 17-Dec-2020
import numpy as np
import particle as P
import nbody as nb
import pickle
from tools import save_frames, create_gif, plot_energy, plot_2D_trajectory



npart = 2
gridsize = 100
size = (gridsize, gridsize, gridsize)

mass = 100
init_mas = [mass for t in range(npart)]

specific_r = np.array([[gridsize / 2 - 5, gridsize / 2, gridsize / 2], [gridsize / 2 + 5, gridsize / 2, gridsize / 2]])
specific_v = np.array([[0, 0.5, 0], [0, -0.5, 0]])

r1 = 20
soft = 1
dt = .3
niter = 400
freq = 5

s = P.system_init(npart, size, init_mas, specific_r=specific_r, soft=soft, specific_v=specific_v, boundary='Periodic')
g = nb.nbody(size, s, dt, soft=soft, boundary='Periodic')

posi = []
energy = []
for i in range(niter):
    p, e = g.evolve()
    posi.append(p)
    energy.append(e)


pickle.dump((posi, energy), open("Q2.p", "wb"))


save_frames("Q2.p", 'Q2', freq, gridsize)
create_gif('Q2', niter, freq, duration=niter / freq)
plot_energy("Q2.p")
plot_2D_trajectory("Q2.p", "Q2_Trajectory.png")
