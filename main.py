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
if random == 1:       # random initial condition (0:off, 1:on)
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
residual = np.zeros((cells, cells, cells))
force = np.zeros((3, cells, cells, cells))

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
  for i in range(particle):
    for j in range(cells):
      for k in range(cells):
        for l in range(cells):
          if abs(r[1][i] - (0.5*dx + j*dx)) < 0.5*dx and abs(r[2][i] - (0.5*dx + k*dx)) < 0.5*dx and abs(r[3][i] - (0.5*dx + l*dx)) < 0.5*dx:
            rho[j][k][l] += m[particle]/dx
            
def CIC()
  for i in range(particle):
    for j in range(cells):
      for k in range(cells):
        for l in range(cells):
          if abs(r[1][i] - (0.5*dx + j*dx)) < dx and abs(r[2][i] - (0.5*dx + k*dx)) < dx and abs(r[3][i] - (0.5*dx + l*dx)) < dx:
            rho[j][k][l] += m[particle] * (1 - abs(r[1][i] - (0.5*dx + j*dx))/dx) * (1 - abs(r[2][i] - (0.5*dx + k*dx))/dx) * (1 - abs(r[3][i] - (0.5*dx + k*dx))/dx) / dx**3

def TSC()
  for i in range(particle):
    for j in range(cells):
      for k in range(cells):
        for l in range(cells):
          if abs(r[1][i] - (0.5*dx + j*dx)) < (1.5*dx) and abs(r[2][i] - (0.5*dx + k*dx)) < (1.5*dx) and abs(r[3][i] - (0.5*dx + l*dx)) < (1.5*dx):
            rho[j][k][l] += m[particle] * ((1 - round(abs(r[1][i] - (0.5*dx + j*dx)) / dx)) * (0.75 - (abs(r[1][i] - (0.5*dx + j*dx))/dx)**2) + round(abs(r[1][i] - (0.5*dx + j*dx)) / dx) * (0.5*(1.5-abs(r[1][i] - (0.5*dx + j*dx))/dx)**2)) \
              * ((1 - round(abs(r[2][i] - (0.5*dx + k*dx)) / dx)) * (0.75 - (abs(r[2][i] - (0.5*dx + k*dx))/dx)**2) + round(abs(r[2][i] - (0.5*dx + k*dx)) / dx) * (0.5*(1.5-abs(r[2][i] - (0.5*dx + k*dx))/dx)**2)) \
              * ((1 - round(abs(r[3][i] - (0.5*dx + l*dx)) / dx)) * (0.75 - (abs(r[3][i] - (0.5*dx + l*dx))/dx)**2) + round(abs(r[3][i] - (0.5*dx + l*dx)) / dx) * (0.5*(1.5-abs(r[3][i] - (0.5*dx + l*dx))/dx)**2)) / dx**3
      
      
# -------------------------------------------------------------------
# Poisson solver
# -------------------------------------------------------------------
# ρ to Φ (parallel later)
relax = 1.6
errorsum = 1

while errorsum > 10**(-12):
  for i in range(cells):
    for j in range(cells):
      for k in range(cells):
        residual[i][j][k] = phi[i-1][j][k] + phi[i+1][j][k] + phi[i][j-1][k] + phi[i][j+1][k] + phi[i][j][k-1] + phi[i][j][k+1] - 6*phi[i][j][k] - rho[i][j][k]*dx*dx
        phi[i][j][k] = phi[i][j][k] + relax * residual[i][j][k] / 6

  errorsum = 0
  for i in range(cells):
    for j in range(cells):
      for k in range(cells):
        errorsum += abs(residual[i][j][k]/phi[i][j][k]) / cells**2
      

# FFT


# inter-particle force
for i in range(cells):
    for j in range(cells):
      for k in range(cells):
        gfield[1][i][j][k] = -(phi[i+1][j][k]-phi[i-1][j][k])*0.5/dx
        gfield[2][i][j][k] = -(phi[i][j+1][k]-phi[i][j-1][k])*0.5/dx
        gfield[3][i][j][k] = -(phi[i][j][k+1]-phi[i][j][k-1])*0.5/dx


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
