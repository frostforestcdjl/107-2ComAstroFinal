import numpy as np
import math
import random
import time
import Poisson_FFT
import matplotlib.animation as animation
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D, get_test_data
import pandas as pd
#--------------------------------------------------------------------
# parameters
#--------------------------------------------------------------------
# constants
L = 1.0              # 1-D conputational domain size
cells = 64           # number of computing cells
N = cells
particle = 10        # particle number in the box
end_time = 0.1       # end time

# derived constants
dx = L/cells         # spatial resolution
dt = dx/(numpy.amax(v))   # time interval for data update (Courant-Friedrichs-Lewy condition)


# -------------------------------------------------------------------
# define initial condition
# -------------------------------------------------------------------
t = 0.0
ran = 0                                      # random initial condition (0:off, 1:on)
# array
if ran == 0:                                 # random initial condition (0:off, 1:on)
  m = np.ones(particle)                      # mass of particles (m = 1)
  r = np.random.rand(3, particle)            # coordinates of particles (x, y, z = 0)
  v = np.zeros((3, particle))                # velocity of particles (vx, vy, vz = 0)
else:
  m = np.random.rand(particle)               # mass of particles (m = [0,1])
  r = np.random.rand(3, particle)            # coordinates of particles (x, y, z = [0,1])
  v = np.random.normal(size=(3, particle))   # velocity of particles (vx, vy, vz = normal distribution)

'''
#Kepler motion
particle = 2
m = np.ones(particle)
r = np.array([[0.8,0.2],[0.5,0.5],[0.5,0.5]])
G = 1.0
M = 2.0
v = np.array([[0.0,0.0],[( G*M/0.3 )**0.5,-( G*M/0.3 )**0.5],[0.0,0.0]])
'''

p0 = m * v                                    # initial momentum
rho = np.zeros((cells+2, cells+2, cells+2))   # empty 3D box of rho (, density)
phi = np.zeros((cells+2, cells+2, cells+2))   # empty 3D box of phi (, potential field)
residual = np.zeros((cells, cells, cells))
force = np.zeros((3, cells, cells, cells))


tStart = time.time()                         # Start timing
# -------------------------------------------------------------------
# Deposit particle mass onto grid (NGP, CIC, TSC)
# -------------------------------------------------------------------
def NGP():
  for i in range(particle):
    for j in range(cells):
      for k in range(cells):
        for l in range(cells):
          if abs(r[0][i] - (1.5*dx + j*dx)) < 0.5*dx and abs(r[1][i] - (1.5*dx + k*dx)) < 0.5*dx and abs(r[2][i] - (1.5*dx + l*dx)) < 0.5*dx:
            rho[j+1][k+1][l+1] += m[i]/dx
            if rho[j+1][k+1][l+1] != 0:
                print j,k,l,rho[j+1][k+1][l+1]

            
def CIC():
  for i in range(particle):
    for j in range(cells):
      for k in range(cells):
        for l in range(cells):
          if abs(r[0][i] - (1.5*dx + j*dx)) < dx and abs(r[1][i] - (1.5*dx + k*dx)) < dx and abs(r[2][i] - (1.5*dx + l*dx)) < dx:
            rho[j+1][k+1][l+1] += m[i] * (1 - abs(r[0][i] - (1.5*dx + j*dx))/dx) * (1 - abs(r[1][i] - (1.5*dx + k*dx))/dx) * (1 - abs(r[2][i] - (1.5*dx + k*dx))/dx) / dx**3
#            if rho[j+1][k+1][l+1] != 0:
#                print j,k,l,rho[j+1][k+1][l+1]

def TSC():
  for i in range(particle):
    for j in range(cells):
      for k in range(cells):
        for l in range(cells):
          if abs(r[0][i] - (1.5*dx + j*dx)) < (1.5*dx) and abs(r[1][i] - (1.5*dx + k*dx)) < (1.5*dx) and abs(r[2][i] - (1.5*dx + l*dx)) < (1.5*dx):
            rho[j][k][l] += m[i] * ((1 - round(abs(r[0][i] - (0.5*dx + j*dx)) / dx)) * (0.75 - (abs(r[0][i] - (0.5*dx + j*dx))/dx)**2) + round(abs(r[0][i] - (0.5*dx + j*dx)) / dx) * (0.5*(1.5-abs(r[0][i] - (1.5*dx + j*dx))/dx)**2)) \
              * ((1 - round(abs(r[1][i] - (1.5*dx + k*dx)) / dx)) * (0.75 - (abs(r[1][i] - (1.5*dx + k*dx))/dx)**2) + round(abs(r[1][i] - (1.5*dx + k*dx)) / dx) * (0.5*(1.5-abs(r[1][i] - (1.5*dx + k*dx))/dx)**2)) \
              * ((1 - round(abs(r[2][i] - (1.5*dx + l*dx)) / dx)) * (0.75 - (abs(r[2][i] - (1.5*dx + l*dx))/dx)**2) + round(abs(r[2][i] - (1.5*dx + l*dx)) / dx) * (0.5*(1.5-abs(r[2][i] - (1.5*dx + l*dx))/dx)**2)) / dx**3
      
NGP()

      
# -------------------------------------------------------------------
# Poisson solver
# -------------------------------------------------------------------
#  (parallel later)
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




# -------------------------------------------------------------------
# Orbit integration (KDK, DKD)
# -------------------------------------------------------------------
def KDK(a):
  global t
  for i in range(particle):
    for j in range(3):
      v[j][i] = v[j][i] + a[j][i]*0.5*dt      # (a) kick: calculate a(t+0.5*dt) and use that to update velocity by 0.5*dt
      r[j][i] = r[j][i] + v[j][i]*dt          # (b) drift: update position by dt
      v[j][i] = v[j][i] + a[j][i]*0.5*dt      # (c) kick: calculate a(t+0.5*dt) and use that to update velocity by 0.5*dt
      v_abs = ( v[0][0]**2 + v[1][0]**2 )**0.5
  E0    = 0.5*v_abs**2 - G*M/0.3
  print 'E0= %s' %E0

      
def DKD(a):
  global t
  for i in range(particle):
    for j in range(3):
      r[j][i] = r[j][i] + v[j][i]*0.5*dt  # (a) drift: update position by 0.5*dt
      v[j][i] = v[j][i] + a[j][i]*dt      # (b) kick: calculate a(t+0.5*dt) and use that to update velocity by dt
      r[j][i] = r[j][i] + v[j][i]*0.5*dt  # (c) drift: update position by 0.5*dt
      

  
tEnd = time.time()                         # End timing

                
# -------------------------------------------------------------------
# Measure the performance scaling
# -------------------------------------------------------------------
print('program time cost: ' + str(tEnd - tStart) + 's')
print('number of cells/particles: ' + str(cells) + '/' + str(particle))

# -------------------------------------------------------------------
# Momentum conservation
# -------------------------------------------------------------------
def momentum_conservation(m,v):
    pt = m * v
    p_diff = 0
    
    for i in range(particle):
        for j in range(3):
            p_diff += pt[j][i] - p0[j][i]
            
    print('momentum difference is: ' + str(p_diff))

# -------------------------------------------------------------------
# Run the simulation
# -------------------------------------------------------------------
def update(m_scheme,v_scheme):
    global rho, t, dt, L, N, dx, particle, cells, r
    rx = np.empty(0)
    ry = np.empty(0)
    rz = np.empty(0)
    vx = np.array([])
    vy = np.array([])
    vz = np.array([])
    while t <= end_time-dt:
        if m_scheme == 'NGP':
            NGP()
            u = Poisson_FFT.FFT_solver(rho,L,N,cells) # potential
            g = Poisson_FFT.gravity(u,L,cells) #gravitational field
            #print g
            a = Poisson_FFT.interpolate_NGP(g,r,dx,particle,cells)
        elif m_scheme == 'CIC':
            CIC()
            u = Poisson_FFT.FFT_solver(rho,L,N,cells) # potential
            g = Poisson_FFT.gravity(u,L,cells) #gravitational field
            a = Poisson_FFT.interpolate_CIC(g,r,dx,particle,cells)
        elif m_scheme == 'TSC':
            TSC()
            u = Poisson_FFT.FFT_solver(rho,L,N,cells) # potential
            g = Poisson_FFT.gravity(u,L,cells) #gravitational field
            a = Poisson_FFT.interpolate_TSC(g,r,dx,particle,cells)
        else:
            print 'Error: invalid scheme name'
            break

        if v_scheme == 'KDK':
            KDK(a)
        elif v_scheme == 'DKD':
            DKD(a)
        else:
            print 'Error: invalid scheme name'
            break
    
        momentum_conservation(m,v)
        print v
        
        t += dt
        print t
        

        
        #save data
        rx = np.concatenate((rx, r[0]),axis=0)
        print r[0]
        ry = np.concatenate((ry, r[1]),axis=0)
        rz = np.concatenate((rz, r[2]),axis=0)
        
    return rx, ry, rz
        

rdata = update('CIC','KDK')
rx = rdata[0]
ry = rdata[1]
rz = rdata[2]

'''
Writer = animation.writers['ffmpeg']
writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)

fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection="3d")
graph = ax.scatter(r[0], r[1], r[2], color='darkblue')
ax.set_xlim3d(0, 1)
ax.set_ylim3d(0, 1)
ax.set_zlim3d(0, 1)

ani = animation.FuncAnimation(fig, update('CIC','KDK'), frames=200, interval=50, blit=False)
#ani.save('PMcode_Ou&Lin.mp4', writer=writer)
plt.show()
'''
r_data = {'rx':rx,
           'ry':ry,
           'rz':rz
        }


df = pd.DataFrame(r_data)
df.to_csv(path_or_buf='test_data.csv',index=True,line_terminator='\n')







# -------------------------------------------------------------------
# Animation and output mp4
# -------------------------------------------------------------------
Writer = animation.writers['ffmpeg']
writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)

def update_r(num):
  global t, r
  KDK()
  graph._offsets3d = (r[0], r[1], r[2])
  return graph,

fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection="3d")
graph = ax.scatter(r[0], r[1], r[2], color='darkblue')
ax.set_xlim3d(-1, 1)
ax.set_ylim3d(-1, 1)
ax.set_zlim3d(-1, 1)

ani = animation.FuncAnimation(fig, update_lines, frames=200, interval=50, blit=False)
ani.save('PMcode_Ou&Lin.mp4', writer=writer)
plt.show()
