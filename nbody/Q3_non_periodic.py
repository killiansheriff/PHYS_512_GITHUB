# @Author: Killian Sheriff <Killian>
# @Date:   12-Dec-2020
# @Email:  killian.sheriff@gmail.com / killian.sheriff@mail.mcgill.ca
# @Last modified by:   Killian
# @Last modified time: 17-Dec-2020
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

size = (gridsize,gridsize,gridsize)
init_v = 0
mass = 1/npart
init_mas = [mass for t in range(npart)]

soft = 0.8
dt = 5
niter = 1500
freq=10

s = P.system_init(npart,size,init_mas,specific_r=None,specific_v=init_v,boundary='Non-Periodic')
g = nb.nbody(size,s,dt,soft=soft,boundary='Non-Periodic')


posi=[]
energy=[]
for i in range(niter):
    p,e=g.evolve()
    posi.append(p)
    energy.append(e)


pickle.dump( (posi,energy), open( "Q3_non_periodic.p", "wb" ) )



save_frames("Q3_non_periodic.p",'Q3_non_periodic',freq,gridsize)
create_gif('Q3_non_periodic',niter,freq,duration=niter/freq)
plot_energy("Q3_non_periodic.p")
