import numpy as np
           
def FFT_solver(rho,L,N):
    dk = L/N
    D = np.fft.fftn( rho,s=(N,N,N) )
    kx = 2*np.pi*np.fft.fftfreq( N, dk )
    ky = 2*np.pi*np.fft.fftfreq( N, dk )
    kz = 2*np.pi*np.fft.fftfreq( N, dk )
    DC = D[0][0][0]
    kxv, kyv, kzv = np.meshgrid(kx, ky, kz)
    phik = -D/(kxv**2+kyv**2+kzv**2)
    phik[0][0][0]=DC
   
    phix = np.fft.ifft2(phik)
    u = np.real(phix)
    return u


def gravity(u,L,cells):
    gx = np.zeros((cells, cells, cells))
    gy = np.zeros((cells, cells, cells))
    gz = np.zeros((cells, cells, cells))
    dx = L/cells
    for j in range(cells):
        for k in range(cells):
            for l in range(cells):
                gx[j][k][l] = -(u[j+2][k+1][l+1]-u[j][k+1][l+1])/2/dx
                gy[j][k][l] = -(u[j+1][k+2][l+1]-u[j+1][k][l+1])/2/dx
                gz[j][k][l] = -(u[j+1][k+1][l+2]-u[j+1][k+1][l])/2/dx
    return gx, gy, gz


def interpolate_NGP(g,r,dx,particle,cells):
  a = np.zeros((3, particle)) 
  for i in range(particle):
    for j in range(cells):
      for k in range(cells):
        for l in range(cells):
          if abs(r[0][i] - (1.5*dx + j*dx)) < 0.5*dx and abs(r[1][i] - (1.5*dx + k*dx)) < 0.5*dx and abs(r[2][i] - (1.5*dx + l*dx)) < 0.5*dx:
            a[0][i] += g[0][j][k][l]
            a[1][i] += g[1][j][k][l]
            a[2][i] += g[2][j][k][l]
            
def interpolate_CIC(g,r,dx,particle,cells):
  a = np.zeros((3, particle)) 
  for i in range(particle):
    for j in range(cells):
      for k in range(cells):
        for l in range(cells):
          if abs(r[0][i] - (1.5*dx + j*dx)) < dx and abs(r[1][i] - (1.5*dx + k*dx)) < dx and abs(r[2][i] - (1.5*dx + l*dx)) < dx:
               a[0][i] += g[0][j][k][l]* (1 - abs(r[0][i] - (1.5*dx + j*dx))/dx) * (1 - abs(r[1][i] - (1.5*dx + k*dx))/dx) * (1 - abs(r[2][i] - (1.5*dx + k*dx))/dx) 
               a[1][i] += g[1][j][k][l]* (1 - abs(r[0][i] - (1.5*dx + j*dx))/dx) * (1 - abs(r[1][i] - (1.5*dx + k*dx))/dx) * (1 - abs(r[2][i] - (1.5*dx + k*dx))/dx) 
               a[2][i] += g[2][j][k][l]* (1 - abs(r[0][i] - (1.5*dx + j*dx))/dx) * (1 - abs(r[1][i] - (1.5*dx + k*dx))/dx) * (1 - abs(r[2][i] - (1.5*dx + k*dx))/dx) 

def interpolate_TSC(g,r,dx,particle,cells):
  a = np.zeros((3, particle)) 
  for i in range(particle):
    for j in range(cells):
      for k in range(cells):
        for l in range(cells):
          if abs(r[0][i] - (1.5*dx + j*dx)) < (1.5*dx) and abs(r[1][i] - (1.5*dx + k*dx)) < (1.5*dx) and abs(r[2][i] - (1.5*dx + l*dx)) < (1.5*dx):
            a[0][i] += g[0][j][k][l]* ((1 - round(abs(r[0][i] - (0.5*dx + j*dx)) / dx)) * (0.75 - (abs(r[0][i] - (0.5*dx + j*dx))/dx)**2) + round(abs(r[0][i] - (0.5*dx + j*dx)) / dx) * (0.5*(1.5-abs(r[0][i] - (1.5*dx + j*dx))/dx)**2)) \
              * ((1 - round(abs(r[1][i] - (1.5*dx + k*dx)) / dx)) * (0.75 - (abs(r[1][i] - (1.5*dx + k*dx))/dx)**2) + round(abs(r[1][i] - (1.5*dx + k*dx)) / dx) * (0.5*(1.5-abs(r[1][i] - (1.5*dx + k*dx))/dx)**2)) \
              * ((1 - round(abs(r[2][i] - (1.5*dx + l*dx)) / dx)) * (0.75 - (abs(r[2][i] - (1.5*dx + l*dx))/dx)**2) + round(abs(r[2][i] - (1.5*dx + l*dx)) / dx) * (0.5*(1.5-abs(r[2][i] - (1.5*dx + l*dx))/dx)**2))
            a[1][i] += g[1][j][k][l]* ((1 - round(abs(r[0][i] - (0.5*dx + j*dx)) / dx)) * (0.75 - (abs(r[0][i] - (0.5*dx + j*dx))/dx)**2) + round(abs(r[0][i] - (0.5*dx + j*dx)) / dx) * (0.5*(1.5-abs(r[0][i] - (1.5*dx + j*dx))/dx)**2)) \
              * ((1 - round(abs(r[1][i] - (1.5*dx + k*dx)) / dx)) * (0.75 - (abs(r[1][i] - (1.5*dx + k*dx))/dx)**2) + round(abs(r[1][i] - (1.5*dx + k*dx)) / dx) * (0.5*(1.5-abs(r[1][i] - (1.5*dx + k*dx))/dx)**2)) \
              * ((1 - round(abs(r[2][i] - (1.5*dx + l*dx)) / dx)) * (0.75 - (abs(r[2][i] - (1.5*dx + l*dx))/dx)**2) + round(abs(r[2][i] - (1.5*dx + l*dx)) / dx) * (0.5*(1.5-abs(r[2][i] - (1.5*dx + l*dx))/dx)**2))
            a[2][i] += g[2][j][k][l]* ((1 - round(abs(r[0][i] - (0.5*dx + j*dx)) / dx)) * (0.75 - (abs(r[0][i] - (0.5*dx + j*dx))/dx)**2) + round(abs(r[0][i] - (0.5*dx + j*dx)) / dx) * (0.5*(1.5-abs(r[0][i] - (1.5*dx + j*dx))/dx)**2)) \
              * ((1 - round(abs(r[1][i] - (1.5*dx + k*dx)) / dx)) * (0.75 - (abs(r[1][i] - (1.5*dx + k*dx))/dx)**2) + round(abs(r[1][i] - (1.5*dx + k*dx)) / dx) * (0.5*(1.5-abs(r[1][i] - (1.5*dx + k*dx))/dx)**2)) \
              * ((1 - round(abs(r[2][i] - (1.5*dx + l*dx)) / dx)) * (0.75 - (abs(r[2][i] - (1.5*dx + l*dx))/dx)**2) + round(abs(r[2][i] - (1.5*dx + l*dx)) / dx) * (0.5*(1.5-abs(r[2][i] - (1.5*dx + l*dx))/dx)**2))
                
                
                
                
                
