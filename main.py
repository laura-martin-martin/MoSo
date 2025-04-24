import numpy as np
import time
from mesh import neighbours, read_mesh
from matrices import matrices
from flux import flux, decomposed_FN_matrix
from initialization import init_sol
from source import source
from output import output_mesh, output_matrix
from RK import RK, rk4a, rk4b
import matplotlib.pyplot as plt
from parameters import num, times, sources, initial_conds
from cfl import dtscale2D

np.set_printoptions(threshold=np.inf, linewidth=np.inf, precision=2, suppress=True)
Numi = num
dim = num.dim
PolyDeg = num.PolyDeg
Nb = num.Nb
N = num.N
FinalTime = times.FinalTime

##Ne = {}                     # Number of elements of highest dimension

# # # Read mesh
#### nodes, elements, Ne = read_msh(f"computation.msh")
nodes, elements, Ne = read_mesh()
# # # Find neighbourgs for comunnications

neighbours(elements)
X,Y, xplot, yplot = output_mesh(nodes, elements, Ne)
nn = 3
XX, YY = np.meshgrid(np.linspace(-1, 1, nn), np.linspace(-1, 1, nn))
(out_index, out_V2D) = output_matrix(nodes, elements, XX, YY)

# Z = (1 - X/2 + X**5 + Y**3) * np.exp(-X**2 - Y**2)
# levels = np.linspace(Z.min(), Z.max(), 7)
u = np.zeros((Ne,N,4))
if (initial_conds.activated=='True'):
    u = init_sol(Ne, X, Y, u)
# storage for low storage RK time stepping
rhsQ = 0*u
resQ = 0*u

print(Ne)
(GF, Bref) = matrices(nodes, elements, Ne)
Fn = decomposed_FN_matrix(nodes, elements, Ne)
dt = dtscale2D(elements,nodes)
if(FinalTime<=0.0):
    if(times.TotalSteps<1):
        print('Error: TotalSteps too small')
        FinalTime = -1
    else:
        FinalTime = dt*times.TotalSteps

plt.figure(figsize=(12, 12))
times = 0
tstep = 0
elapsed_uext = 0
elapsed_rhsQ = 0
while (times<FinalTime):

    print(tstep)
    # check to see if we need to adjust for final time step
    if(times+dt>FinalTime):
        dt = FinalTime-times
    
    for INTRK in range(5):    
        # compute right hand side of compressible Euler equations
        t = time.time()
        uext = flux( elements, u)
        elapsed_uext = elapsed_uext + time.time() - t
        
        t = time.time()
        rhsQ = RK(Bref, GF, Fn, u, uext)
        elapsed_rhsQ = elapsed_rhsQ + time.time() - t

        if (sources.activated=='True'):
            rhsQ = rhsQ + source(Ne, X, Y, times)

        # initiate and increment Runge-Kutta residuals
        resQ = rk4a[INTRK]*resQ + dt*rhsQ
        
        # update fields
        u = u+rk4b[INTRK]*resQ
    
    # Increment time and compute new timestep
    times = times+dt
    tstep = tstep+1
    if(tstep%3==0):
        plt.clf()
        # # plt.subplot(2,2,1)
        print(np.shape(np.einsum('ij,ijk->ik',out_V2D,u[out_index])))
        sol = np.reshape(np.einsum('ij,ijk->ik',out_V2D,u[out_index]),(nn,nn,4),order='C')
        plt.contourf(XX,YY,sol[:,:,0], np.linspace(-2,2,50), cmap='seismic')
        #plt.scatter(xplot,yplot,s=10, cmap='seismic', c=np.reshape(u[:,:,0],(Ne*N,)),vmin=-1, vmax=1)
        # plt.xlim([0, 1])
        # plt.ylim([-0.5, 0.5])
        plt.colorbar()
        # # plt.subplot(2,2,2)
        # # plt.scatter(xplot,yplot,s=10, cmap='seismic', c=np.reshape(u[:,:,1],(Ne*N,)),vmin=-0.01, vmax=0.01)
        # # #plt.xlim([-1, 0])
        # # #plt.ylim([-1.2, 0.2])
        # # plt.colorbar()
        # # plt.subplot(2,2,3)
        # # plt.scatter(xplot,yplot,s=10, cmap='seismic', c=np.reshape(u[:,:,2],(Ne*N,)),vmin=-0.01, vmax=0.01)
        # # #plt.xlim([-1, 0])
        # # #plt.ylim([-1.2, 0.2])
        # # plt.colorbar()
        # # plt.subplot(2,2,4)
        # # plt.scatter(xplot,yplot,s=10, cmap='seismic', c=np.reshape(u[:,:,3],(Ne*N,)),vmin=-1, vmax=1)
        # # #plt.xlim([-1, 0])
        # # #plt.ylim([-1.2, 0.2])
        # # plt.colorbar()
        print('t_uext: ',elapsed_uext/tstep)
        print('t_rhsQ: ',elapsed_rhsQ/tstep)
        plt.pause(1)
    #plt.waitforbuttonpress()

plt.pause(10)

