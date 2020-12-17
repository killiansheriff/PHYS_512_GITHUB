# @Author: Killian Sheriff <Killian>
# @Date:   03-Dec-2020
# @Email:  killian.sheriff@gmail.com / killian.sheriff@mail.mcgill.ca
# @Last modified by:   Killian
# @Last modified time: 17-Dec-2020
import random as rn
import numpy as np


class particle:
    def __init__(self, m, x, y, z, vx=0, vy=0, vz=0):
        '''Intinialize particle with mass, positon and velocity.'''
        self.mass = m
        self.position = (x, y, z)
        self.velocity = (vx, vy, vz)


class system_init:
    def __init__(self, npart, size, init_mass, specific_r=None, specific_v=None, boundary='Periodic', soft=None, power_law=False):
        '''Define system of particles.'''
        self.boundary = boundary
        self.power_law = power_law

        def random_generator(random_func, factor_arr, size, npart, specific_cond, cond=False):
            '''Generate missing intitial conditions from specific_cond by using random_func.'''
            init_cond, xmin, xmax = [], 1, size[0]-1
            if specific_cond is None:
                q = 0
            else:
                q = len(specific_cond)
                for k in range(q):
                    if cond and (specific_cond[k][0].max() > xmax or specific_cond[k][1].max() > xmax or specific_cond[k][2].max() > xmax or specific_cond[k][0].min() < xmin or specific_cond[k][1].min() < xmin or specific_cond[k][2].min() < xmin):
                        raise ValueError(f'Positions must be within 1 and {size[0]-1}')
                    pos = (specific_cond[k][0], specific_cond[k][1], specific_cond[k][2])
                    init_cond.append(pos)
            for p in range(npart - q):
                if not cond:
                    pos = ((random_func()) * factor_arr[0], (random_func()) * factor_arr[1], (random_func()) * factor_arr[2])
                    init_cond.append(pos)
                else:
                    pos = (random_func(1.001, size[0] - 1.001) * factor_arr[0], random_func(1.001, size[1] -
                                                                                            1.001 )* factor_arr[0], random_func(1.001, size[2] - 1.001) * factor_arr[0])
                    init_cond.append(pos)

            return init_cond

        def initial_r(size, npart, specific_cond):
            '''Generate random missing intial positons.'''

            if self.boundary == 'Periodic':
                # random positions
                init_cond = random_generator(rn.random, [size[0] - 1, size[1] - 1, size[2] - 1], size, npart, specific_cond)

            elif self.boundary == 'Non-Periodic':
                # generating between  1 and size[0]-1
                init_cond = random_generator(rn.uniform, [1, 1, 1], size, npart, specific_cond, cond=True)

            return init_cond

        def initial_v(size, npart, max_vel, specific_cond=None):
            '''Generate random missing intial velocities.'''
            init_cond = []
            if np.isscalar(specific_cond) == False:
                # random velocities based on Gaussian with standard deviation of 1
                init_cond = random_generator(np.random.normal, [max_vel, max_vel, max_vel], size, npart, specific_cond)

            else:
                for p in range(npart):
                    # set all velocities to 0
                    vel = (0, 0, 0)
                    init_cond.append(vel)
            return init_cond

        def initial_mass(size, init_mass, r=None, soft=None):
            '''Early Universe mass generations'''
            if self.power_law == False:
                return np.array([init_mass.copy()]).T
            else:
                r = np.rint(r).astype('int') % size[0]

                # find the power spectrum
                k_x = np.real(np.fft.fft(r[:, 0]))
                k_y = np.real(np.fft.fft(r[:, 1]))
                k_z = np.real(np.fft.fft(r[:, 2]))
                k_total = np.sqrt(k_x**2 + k_y**2 + k_z**2)
                k_total[k_total < soft] = soft
                m = init_mass / k_total**3

                return np.array([m.copy()]).T

        self.npart = npart
        init_r = initial_r(size, npart, specific_cond=specific_r)
        init_v = initial_v(size, npart, 1, specific_cond=specific_v)

        self.particles = np.asarray([particle(m, x[0], x[1], x[2], vx=v[0], vy=v[1], vz=v[2]) for m, x, v in zip(init_mass, init_r, init_v)])

        self.velocities = np.asarray([self.particles[i].velocity for i in range(npart)])

        self.pos = np.asarray([self.particles[i].position for i in range(npart)])

        self.masses = initial_mass(size, init_mass, r=self.pos, soft=soft)
