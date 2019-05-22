import numpy as np
import math
import random
import time

#--------------------------------------------------------------------
# parameters
#--------------------------------------------------------------------
# constants
L = 1.0              # 1-D conputational domain size
cells = 16           # number of computing cells
particle = 10        # particle number in the box
dt = 1.0e-2          # time interval for data update
end_time = 0.1       # end time

# derived constants
dx = L/cells         # spatial resolution

# -------------------------------------------------------------------
# define initial condition
# -------------------------------------------------------------------
random = 0                                   # random initial condition (0:off, 1:on)
# array
if random == 1:                              # random initial condition (0:off, 1:on)
  m = np.ones(particle)                      # mass of particles (m = 1)
  r = np.zeros((3, particle))                # coordinates of particles (x, y, z = 0)
  v = np.zeros((3, particle))                # velocity of particles (vx, vy, vz = 0)
else:
  m = np.random.rand(particle)               # mass of particles (m = [0,1])
  r = np.random.rand(3, particle)            # coordinates of particles (x, y, z = [0,1])
  v = np.random.normal(size=(3, particle))   # velocity of particles (vx, vy, vz = normal distribution)

a = np.zeros((3, particle))                  # acceleration of particles (ax, ay, az = 0)
p0 = m * v                                   # initial momentum

rho = np.zeros((cells+2, cells+2, cells+2))  # empty 3D box of rho (ρ, density)
phi = np.zeros((cells+2, cells+2, cells+2))  # empty 3D box of phi (Φ, potential field)
residual = np.zeros((cells, cells, cells))
force = np.zeros((3, cells, cells, cells))


tStart = time.time()                         # Start timing
# -------------------------------------------------------------------
# Deposit particle mass onto grid (NGP, CIC, TSC)
# -------------------------------------------------------------------
def NGP()
  for i in range(particle):
    for j in range(cells):
      for k in range(cells):
        for l in range(cells):
          if abs(r[1][i] - (1.5*dx + j*dx)) < 0.5*dx and abs(r[2][i] - (1.5*dx + k*dx)) < 0.5*dx and abs(r[3][i] - (1.5*dx + l*dx)) < 0.5*dx:
            rho[j+1][k+1][l+1] += m[particle]/dx
            
def CIC()
  for i in range(particle):
    for j in range(cells):
      for k in range(cells):
        for l in range(cells):
          if abs(r[1][i] - (1.5*dx + j*dx)) < dx and abs(r[2][i] - (1.5*dx + k*dx)) < dx and abs(r[3][i] - (1.5*dx + l*dx)) < dx:
            rho[j+1][k+1][l+1] += m[particle] * (1 - abs(r[1][i] - (1.5*dx + j*dx))/dx) * (1 - abs(r[2][i] - (1.5*dx + k*dx))/dx) * (1 - abs(r[3][i] - (1.5*dx + k*dx))/dx) / dx**3

def TSC()
  for i in range(particle):
    for j in range(cells):
      for k in range(cells):
        for l in range(cells):
          if abs(r[1][i] - (1.5*dx + j*dx)) < (1.5*dx) and abs(r[2][i] - (1.5*dx + k*dx)) < (1.5*dx) and abs(r[3][i] - (1.5*dx + l*dx)) < (1.5*dx):
            rho[j][k][l] += m[particle] * ((1 - round(abs(r[1][i] - (0.5*dx + j*dx)) / dx)) * (0.75 - (abs(r[1][i] - (0.5*dx + j*dx))/dx)**2) + round(abs(r[1][i] - (0.5*dx + j*dx)) / dx) * (0.5*(1.5-abs(r[1][i] - (1.5*dx + j*dx))/dx)**2)) \
              * ((1 - round(abs(r[2][i] - (1.5*dx + k*dx)) / dx)) * (0.75 - (abs(r[2][i] - (1.5*dx + k*dx))/dx)**2) + round(abs(r[2][i] - (1.5*dx + k*dx)) / dx) * (0.5*(1.5-abs(r[2][i] - (1.5*dx + k*dx))/dx)**2)) \
              * ((1 - round(abs(r[3][i] - (1.5*dx + l*dx)) / dx)) * (0.75 - (abs(r[3][i] - (1.5*dx + l*dx))/dx)**2) + round(abs(r[3][i] - (1.5*dx + l*dx)) / dx) * (0.5*(1.5-abs(r[3][i] - (1.5*dx + l*dx))/dx)**2)) / dx**3
      
      
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
        residual[i][j][k] = phi[i][j+1][k+1] + phi[i+2][j+1][k+1] + phi[i+1][j][k+1] + phi[i+1][j+2][k+1] + phi[i+1][j+1][k] + phi[i+1][j+1][k+2] - 6*phi[i+1][j+1][k+1] - rho[i][j][k]*dx*dx
        phi[i+1][j+1][k+1] = phi[i+1][j+1][k+1] + relax * residual[i][j][k] / 6
  
  errorsum = 0
  for i in range(cells):
    for j in range(cells):
      for k in range(cells):
        errorsum += abs(residual[i][j][k]/phi[i][j][k]) / cells**2
      

# FFT
rhok = np.fft.rfft( rho )
k = 0.5           #!!need to be replaced!!
phiF = np.fft.irfft( -rhok/k**2 )


# inter-particle force
for i in range(cells):
    for j in range(cells):
      for k in range(cells):
        gfield[1][i][j][k] = -(phi[i+2][j+1][k+1]-phi[i][j+1][k+1])*0.5/dx
        gfield[2][i][j][k] = -(phi[i+1][j+2][k+1]-phi[i+1][j][k+1])*0.5/dx
        gfield[3][i][j][k] = -(phi[i+1][j+1][k+2]-phi[i+1][j+1][k])*0.5/dx


# -------------------------------------------------------------------
# Orbit integration (KDK, DKD)
# -------------------------------------------------------------------
def KDK()
  for i in range(particle):
    for j in range(3):
      v[j][i] = v[j][i] + a[j][i]*0.5*dt      # (a) kick: calculate a(t+0.5*dt) and use that to update velocity by 0.5*dt
      r[j][i] = r[j][i] + v[j][i]*dt          # (b) drift: update position by dt
      v[j][i] = v[j][i] + a[j][i]*0.5*dt      # (c) kick: calculate a(t+0.5*dt) and use that to update velocity by 0.5*dt
      
  # update time
  t = t + dt
  if (t >= end_time): 
    break
      
def DKD()
  for i in range(particle):
    for j in range(3):
      r[j][i] = r[j][i] + v[j][i]*0.5*dt  # (a) drift: update position by 0.5*dt
      v[j][i] = v[j][i] + a[j][i]*dt      # (b) kick: calculate a(t+0.5*dt) and use that to update velocity by dt
      r[j][i] = r[j][i] + v[j][i]*0.5*dt  # (c) drift: update position by 0.5*dt
      
  # update time
  t = t + dt
  if (t >= end_time): 
    break
  
tEnd = time.time()                         # End timing
# -------------------------------------------------------------------
# Measure the performance scaling
# -------------------------------------------------------------------
print('program time cost: ' + str(tEnd - tStart) + 's')
print('number of cells/particles: ' + str(cells) + '/' + str(particles))

# -------------------------------------------------------------------
# Momentum conservation
# -------------------------------------------------------------------
pt = m * v
p_diff = 0

for i in range(particle):
  for j in range(3):
    p_diff += pt[j][i] - p0[j][i]
    
print('momentum difference is: ' + str(p_diff))

