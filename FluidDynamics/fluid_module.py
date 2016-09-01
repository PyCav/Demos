import numpy as np
from copy import copy

def define_N(Nx,Ny):
    global N_x,N_y
    N_x = Nx
    N_y = Ny

def swap(x,x0):
    x_tmp = copy(x0)
    x0    = copy(x )
    x     = x_tmp
    return x,x0

def initialise():
    v_x  = np.zeros(((N_x+2),(N_y+2)))
    v_y  = np.zeros(((N_x+2),(N_y+2)))
    v_x0 = np.zeros(((N_x+2),(N_y+2)))
    v_y0 = np.zeros(((N_x+2),(N_y+2)))
    
    n    = np.zeros(((N_x+2),(N_y+2)))
    n0   = np.zeros(((N_x+2),(N_y+2)))

    index = int(N_x/2)

    #v_y0[index-1:index+1,index] = 100
    #v_y[index-1:index+1,index] = 100

    return n,n0,v_x,v_y,v_x0,v_y0

def gauss_siedal(x,x0,a):
    if a == 1:
        divisor = 4
    else:
        divisor = (1+4*a)

    for k in range(20):
        x[1:-1,1:-1] = (x0[1:-1,1:-1]+a*(x[:-2,1:-1]+x[2:,1:-1]+x[1:-1,:-2]+x[1:-1,2:]))/divisor
    return x

def forces(n,u0,v0,s):
    index = int(N_x/4)
    n[index-1:index+1,index-1:index+1] += 1
    n[3*index-1:3*index+1,index-1:index+1] += 1
    #index = int(N_x/2)
    v0[1:-1,1:-1] = n[1:-1,1:-1]
    u0[1:-1,1:-1] = s*(np.random.rand(N_x,N_y)-0.5)
    return n,u0,v0

def diffuse(b,n,n0,diff,dt):
    a = dt*diff*N_x*N_y

    n = gauss_siedal(n,n0,a)

    n = set_bnd(b,n)

    return n

def advect(b,d,d0,u,v,dt):
    dt0_x,dt0_y = dt*N_x,dt*N_y

    i,j = np.meshgrid(np.arange(1,N_x+1),np.arange(1,N_y+1), indexing = 'ij')

    x = i-dt0_x*u[1:-1,1:-1]
    y = j-dt0_y*v[1:-1,1:-1]
    
    x[np.where(x < 0.5)] = 0.5
    x[np.where(x > N_x+0.5)] = N_x+0.5

    y[np.where(y < 0.5)] = 0.5
    y[np.where(y > N_y+0.5)] = N_y+0.5
                
    i0,j0 = x.astype('int'),y.astype('int')
    i1,j1 = i0+1,j0+1
            
    s1,t1 = x-i0,y-j0
    s0,t0 = 1-s1,1-t1

    d[1:-1,1:-1] = (s0*(t0*d0[i0,j0]+t1*d0[i0,j1])+s1*(t0*d0[i1,j0]+t1*d0[i1,j1]))

    d = set_bnd(b,d)

    return d

def project(u,v):
    h_x = 1.0/N_x
    h_y = 1.0/N_y
    p   = np.zeros_like(u)
    div = np.zeros_like(u)

    div[1:-1,1:-1] = -0.5*(h_x*(u[2:,1:-1]-u[:-2,1:-1])+h_y*(v[1:-1,2:]-v[1:-1,:-2]))

    set_bnd(0,div)

    p = gauss_siedal(p,div,1)   
             
    set_bnd(0,p)

    u[1:-1,1:-1] -= 0.5*(p[2:,1:-1]-p[:-2,1:-1])/h_x
    v[1:-1,1:-1] -= 0.5*(p[1:-1,2:]-p[1:-1,:-2])/h_y

    set_bnd(1,u)
    set_bnd(2,v)

    return u,v

def dens_step(x,x0,u,v,diff,dt):

    x,x0 = swap(x,x0)
    x = diffuse(0, x, x0, diff, dt )
    
    x,x0 = swap(x,x0)
    x = advect(0, x, x0, u, v, dt )

    return x,x0

def vel_step(u,v,u0,v0,visc,dt):
    u += u0*dt
    v += v0*dt

    u,u0 = swap(u,u0)
    u = diffuse(1, u, u0, visc, dt )

    v,v0 = swap(v,v0)
    v = diffuse(2, v, v0, visc, dt )
    u,v = project(u, v)
                        
    u,u0 = swap(u,u0)
    v,v0 = swap(v,v0)
    u = advect(1, u, u0, u0, v0, dt )
    v = advect(2, v, v0, u0, v0, dt )
                        
    u,v = project(u, v)

    return u,v,u0,v0
    
def set_bnd(b,x):
    x[0,1:-1]    = x[1,1:-1]
    x[-1,1:-1]   = x[-2,1:-1]

    x[1:-1,0]    = x[1:-1,1]
    x[1:-1,-1]   = x[1:-1,-2]

    x[0 ,0 ] = 0.5*(x[1,0]+x[0 ,1])
    x[0 ,-1] = 0.5*(x[1,-1]+x[0 ,-2])
    x[-1,0 ] = 0.5*(x[-2,0]+x[-1,1])
    x[-1,-1] = 0.5*(x[-2,-1]+x[-1,-2])

    return x