import numpy as np
import math
import random

#--------------------------------------------------------------------
# parameters
#--------------------------------------------------------------------
# constants
L = 1.0              # 1-D conputational domain size
cells = 16           # number of computing cells
particle = 10        # particle number in the box
random = 0           # random initial condition (0:off, 1:on)

# array
if random = 1:       # random initial condition (0:off, 1:on)
  m = np.zeros(particle)          # mass of particles (m)
  r = np.zeros((3, particle))     # coordinates of particles (x, y, z)
  v = np.zeros((3, particle))     # velocity of particles (vx, vy, vz)
  a = np.zeros((3, particle))     # acceleration of particles (ax, ay, az)
else:
  m = np.random.rand(particle)
  r = np.random.rand(3, particle)
  v = np.random.rand(3, particle)
  a = np.random.rand(3, particle)

p = np.zeros((3, particle))       # momentum of particles (px, py, pz)
  
rho = np.zeros((cells, cells, cells))      # empty 3D box of rho (ρ, density)
phi = np.zeros((cells, cells, cells))      # empty 3D box of phi (Φ, potential field)

# derived constants
dx = L/cells         # spatial resolution
for i in range(3):   # p = m*v
  for j in range(particle):
    p[i][j] = m[j]*v[i][j]


# -------------------------------------------------------------------
# define initial condition
# -------------------------------------------------------------------


# -------------------------------------------------------------------
# Deposit particle mass onto grid (NGP, CIC, TSC)
# -------------------------------------------------------------------
def NGP()


def CIC()


def TSC()
    


# -------------------------------------------------------------------
# Poisson solver
# -------------------------------------------------------------------
# ρ to Φ


# FFT


# inter-particle force

# -------------------------------------------------------------------
# Orbit integration (KDK, DKD)
# -------------------------------------------------------------------
def KDK()



def DKD()


# -------------------------------------------------------------------
# Measure the performance scaling
# -------------------------------------------------------------------




# -------------------------------------------------------------------
# Momentum conservation
# -------------------------------------------------------------------
