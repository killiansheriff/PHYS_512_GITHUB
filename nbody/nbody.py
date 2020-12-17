# @Author: Killian Sheriff <Killian>
# @Date:   03-Dec-2020
# @Email:  killian.sheriff@gmail.com / killian.sheriff@mail.mcgill.ca
# @Last modified by:   Killian
# @Last modified time: 17-Dec-2020
import tools
import numpy as np


class nbody:
    '''3D Nbody class for simulations.'''

    def __init__(self, size, particles_arr, dt, soft=0.1, G=1, boundary='Periodic'):
        '''Simulation parameters intialization.'''

        self.dt, self.soft, self.G, self.boundary, self.npart = dt, soft, G, boundary, particles_arr.npart
        self.r, self.v, self.mass = particles_arr.pos, particles_arr.velocities, particles_arr.masses
        self.init_size = size
        self.energy = 0

        if self.boundary == 'Non-Periodic':
            self.size = (size[0], size[1], size[2])  # could potentially implement guard cell to avoid wrap around effect in convolution
        else:
            self.size = size

        x, y, z = np.arange(self.size[0], dtype=float), np.arange(self.size[1], dtype=float), np.arange(self.size[2], dtype=float)
        self.mesh = np.array(np.meshgrid(x, y, z))

        self.get_densities()
        self.get_green()

    def get_densities(self):
        '''Assigns grid density according to the Nearest Grid Point (NGP) scheme.'''

        # assigning to nearest grid point
        self.int_r = np.rint(self.r).astype('int')  # np.rint() faster than np.round()
        self.mesh_modified = tuple(self.int_r[:, i] for i in range(3))

        # bin masses to their nearest gridpoint
        E = np.linspace(0, self.size[0] - 1, num=self.size[0] + 1)
        E = np.repeat([E], 3, axis=0)
        hist = np.histogramdd(self.int_r, bins=E, weights=self.mass.flatten())  # could implement a faster version using numba

        self.densities = hist[0]

    def get_green(self):
        'Defines the green function. Will only be relevent for even grid sizes.'
        if self.size[0] % 2 == 0:
            r = np.sum(self.mesh**2, axis=0)
            r[r < self.soft**2] = self.soft**2
            r += self.soft**2
            g = 1 / (4 * np.pi * np.sqrt(r))

            # For 2D, to get periodicity, we flip the corners to get the same behavior around
            #n_x, n_y,n_z = self.size[0], self.size[1],self.size[2]

            # try:
            #    g[n_x:, :n_y] = np.flip(g[:n_x, :n_y], axis=0)
            #    g[:, n_y:] = np.flip(g[:, :n_y], axis=1)
            # except:
            #    g[n_x:, :n_y + 1] = np.flip(g[:n_x + 1, :n_y + 1], axis=0)
            #    g[:, n_y:] = np.flip(g[:, :n_y + 1], axis=1)

            # generalized behvaiour to 3D with np.flip
            g += np.flip(g, 0)
            g += np.flip(g, 1)
            g += np.flip(g, 2)

            self.g = g
        else:
            print('Please use an even grid size!')
            return

    def get_pot(self):
        '''Get the potential by convuluting grid density with grid green function.'''

        ffD = np.fft.rfftn(self.densities)
        ffG = np.fft.rfftn(self.g)
        ffV = ffD * ffG

        V = np.fft.irfftn(ffV)

        #    Unforntunatly with padding for non-periodic BC to run is way longer... so I couldn't run it :(
        #    V=self.pad_conv(self.densities,self.g)
        #    try padding with guard cell
        #    density_padded = np.zeros( (self.size[0]-1,self.size[0]-1,self.size[0]-1))
        #    density_padded[0:self.size[0]-2, self.size[0]-2, 0:self.size[0]-2] = self.densities
        #    fft_greens = np.fft.fftn(self.g)
        #    fft_density = np.fft.fftn(density_padded)
        #    V = np.fft.ifftn(fft_density * fft_greens)
        #    V=self.pad_conv(self.densities,self.g)

        # shift and average the potential to center it back to the particle
        for i in range(3):
            V = 0.5 * (np.roll(V, 1, axis=i) + V)
        # usefull for an eventuel guard cell implementation / padding
        V = V[:self.init_size[0], :self.init_size[1], :self.init_size[2]]

        # set potential to 0 on the boundary of the simualtion box
        if self.boundary == 'Non-Periodic':
            V[0:, 0, 0] = 0
            V[0:, -1, 0] = 0
            V[0:, -1, -1] = 0
            V[0:, 0, -1] = 0
            V[0, 0:, 0] = 0
            V[-1, 0:, 0] = 0
            V[0, 0:, -1] = 0
            V[-1, 0:, -1] = 0
            V[0, 0, 0:] = 0
            V[-1, 0, 0:] = 0
            V[0, -1, 0:] = 0
            V[-1, -1, 0:] = 0
            V[-1, -1, -1] = 0

        return V

    def get_forces_mesh(self):
        '''Using central difference to take the gradiant of the potential and get forces.'''

        V = self.get_pot()
        fgrid = np.zeros([3, self.size[0], self.size[1], self.size[2]])

        for i in range(3):
            fgrid[i] = 0.5 * (np.roll(V, 1, axis=i) - np.roll(V, -1, axis=i))

        fgrid = -fgrid * self.densities * self.G
        return fgrid

    def get_forces_particles(self):
        '''Use the inverse NGP schme to interpolate the force back to the particles.'''
        fxyz = np.moveaxis(self.get_forces_mesh(), 0, -1)
        f_ptcls = fxyz[self.mesh_modified]
        return f_ptcls

    def totalEnergy(self):
        '''Get total energy of the system.'''
        K = np.sum(self.mass * self.v**2)
        P = -0.5 * np.sum(np.sum(self.get_pot()) * self.densities)
        T = K + P
        self.energy = T

    def evolve(self, nsteps=1, file_save=None, file_save_pos=None):
        '''Ecolve the system using a leapfrog solver withg fixed timestep.'''
        F = self.get_forces_particles()
        self.totalEnergy()
        # could implement an higher order scheme with accelerations
        self.r, self.v = tools.evolve(self.r, self.v, self.mass, F, self.dt, self.size[0])

        # boundary conditions check
        if self.boundary == "Periodic":
            self.r = self.r % (self.init_size[0] - 1)
        if self.boundary == "Non-Periodic":
            # delete particle if goes out of the box
            ind_top = np.argwhere((self.r > self.init_size[0] - 1))
            ind_bot = np.argwhere((self.r < 0))
            ind = [i for i in np.append(ind_top, ind_bot)]
            self.v = np.delete(self.v, ind, axis=0)
            self.mass = np.delete(self.mass, ind, axis=0)
            self.r = np.delete(self.r, ind, axis=0)

        self.get_densities()
        self.get_green()
        return self.r, self.energy
