import math
import numpy as np
from parameters import Baseflow, num, conditions
from mesh import nodes_I2D, nodes_I3D

dim = num.dim
PolyDeg = num.PolyDeg
Nb = num.Nb
N = num.N

sign = lambda x: math.copysign(1, x)

if (dim==2):
    def flux(elements, u):   ## Comunicate the solution at the border of the neighbourging elements 
        (indexB) = nodes_I2D(PolyDeg)
        # if (label=='line2'):
        # elements[label]["neigh_element_type"] = np.chararray((len(tags),1), itemsize=10)
        # elements[label]["neigh_element_type"][:] = 'triangle3'
        # elements[label]["neigh_element"] = np.zeros_like(tags, dtype=int)
        # elements[label]["neigh_side"] = np.zeros_like(tags, dtype=int)
        # else:
        # elements[label]["neigh_element_type"] = np.chararray((len(tags),3), itemsize=10)
        # elements[label]["neigh_element_type"][:] = 'triangle3'
        # elements[label]["neigh_element"] = np.zeros((len(tags),3), dtype=int)
        # elements[label]["neigh_side"] = np.zeros((len(tags),3), dtype=int)

        
        # u: Ne x N x 4
        # uext: Ne x (dim+1) x Nb x 4
        n_1=len(elements['triangle3']['nodes'])
        uext = np.zeros((n_1,dim+1,2,N,4))
        for i in range(n_1):
            for j in range(3):
                # uext[i,j,k,:] = u[i,j,k,:]
                typeelext = elements['triangle3']["neigh_element_type"][i][j].decode("utf-8")
                if (typeelext=='triangle3'):
                    indx=indexB[abs(elements['triangle3']["neigh_side"][i][j])-1]
                    ## uext = [solution in the element i,  solution on the neighbourg element]
                    uext[i,j,0,indexB[j],:] = u[i,indexB[j],:]                                         
                    uext[i,j,1,indexB[j],:] = u[elements['triangle3']["neigh_element"][i][j],indx[::int(sign(elements['triangle3']["neigh_side"][i][j]))],:]

                else:
                    if (conditions.boundary=='NRBC'):             ##non-reflecting boundary condition
                        uext[i,j,0,indexB[j],:] = u[i,indexB[j],:]                                         
                        uext[i,j,1,indexB[j],:] = 0*u[i,indexB[j],:] 
                    elif (conditions.boundary=='wall'):           ## Wall boundary condition
                        uext[i,j,0,indexB[j],0] = u[i,indexB[j],0]  
                        uext[i,j,0,indexB[j],-1] = u[i,indexB[j],-1]  
                        uext[i,j,1,indexB[j],:] = uext[i,j,0,indexB[j],:]
        return (uext)

    def flux_matrices():
        Fx = np.array([[                                   0,                          1,                          0,                                                  0],
                    [               -Baseflow.u0*Baseflow.u0,              2*Baseflow.u0,                          0, Baseflow.c0*Baseflow.c0*Baseflow.rho0/Baseflow.pi0],
                    [               -Baseflow.u0*Baseflow.v0,                Baseflow.v0,                Baseflow.u0,                                                  0],
                    [-Baseflow.pi0*Baseflow.u0/Baseflow.rho0, Baseflow.pi0/Baseflow.rho0,                          0,                                        Baseflow.u0]])
        
        Fy = np.array([[                                   0,                          0,                          1,                                                  0],
                    [               -Baseflow.u0*Baseflow.v0,                Baseflow.v0,                Baseflow.u0,                                                  0],
                    [               -Baseflow.v0*Baseflow.v0,                          0,              2*Baseflow.v0, Baseflow.c0*Baseflow.c0*Baseflow.rho0/Baseflow.pi0],
                    [-Baseflow.pi0*Baseflow.u0/Baseflow.rho0,                          0, Baseflow.pi0/Baseflow.rho0,                                        Baseflow.v0]])
        return (Fx,Fy)
    
    def fluxN_matrix(nx, ny):
        un = Baseflow.u0*nx+Baseflow.v0*ny
        W = [[                          1,      1,     0,                          1],
                    [Baseflow.u0-nx*Baseflow.c0,  un*nx,   -ny, Baseflow.u0+nx*Baseflow.c0],
                    [Baseflow.v0-ny*Baseflow.c0,  un*ny,    nx, Baseflow.v0+ny*Baseflow.c0],
                    [Baseflow.pi0/Baseflow.rho0,      0,     0, Baseflow.pi0/Baseflow.rho0]]
        W = np.array(W, dtype=float)
        if(conditions.flux=="center"):
            D =  np.diag(np.array([un-Baseflow.c0, un, un, un+Baseflow.c0]))
            FnP = np.dot(np.dot(W, D), np.linalg.inv(W))/2  ## so Fn original is the sum of FnP + FnN
            FnN = np.dot(np.dot(W, D), np.linalg.inv(W))/2
        if(conditions.flux=="upwind"):
            eigenvalP = np.array([un-Baseflow.c0, un, un, un+Baseflow.c0])
            eigenvalN = np.array([un-Baseflow.c0, un, un, un+Baseflow.c0])
            eigenvalP[eigenvalP<0]=0
            eigenvalN[eigenvalN>0]=0
            DP =  np.diag(eigenvalP)
            DN =  np.diag(eigenvalN)
            FnP = np.dot(np.dot(W, DP), np.linalg.inv(W)) # W*D*inv(W)
            FnN = np.dot(np.dot(W, DN), np.linalg.inv(W))
        # Fn = np.array([[           0,          nx,          ny,                 0],
        #                [      -Baseflow.u0*un,    un+Baseflow.u0*nx,       Baseflow.u0*ny, Baseflow.rho0*Baseflow.c0*Baseflow.c0*nx/Baseflow.pi0],
        #                [      -Baseflow.v0*un,       Baseflow.v0*nx,    un+Baseflow.v0*ny, Baseflow.rho0*Baseflow.c0*Baseflow.c0*ny/Baseflow.pi0],
        #                [-Baseflow.pi0*un/Baseflow.rho0, Baseflow.pi0*nx/Baseflow.rho0, Baseflow.pi0*ny/Baseflow.rho0,                un]])
        return (FnP, FnN)

    def decomposed_FN_matrix(nodes, elements, Ne):
        Jeb = np.zeros((dim+1))
        normal = np.zeros((dim+1,dim))
        Fn = np.zeros((Ne,dim+1,2,4,4))
        for i in range(Ne):
            ## positions of the elements' vertices
            [v1x, v1y, _] = nodes["coordinates"][int(elements['triangle3']["nodes"][i][0]-1)]
            [v2x, v2y, _] = nodes["coordinates"][int(elements['triangle3']["nodes"][i][1]-1)]
            [v3x, v3y, _] = nodes["coordinates"][int(elements['triangle3']["nodes"][i][2]-1)]
            ## Mtemp = [[dxdr, dxds,
            ##          [dydr, dyds]]
            Mtemp = [[(v2x-v1x)/2, (v3x-v1x)/2],[(v2y-v1y)/2, (v3y-v1y)/2]]
            Je = np.linalg.det(Mtemp)
            normal[0,:] = [-(-(v2y-v1y)), -(v2x-v1x)]/np.linalg.norm([-(v2y-v1y), v2x-v1x])
            normal[1,:] = [-(-(v3y-v2y)), -(v3x-v2x)]/np.linalg.norm([-(v3y-v2y), v3x-v2x])
            normal[2,:] = [-(-(v1y-v3y)), -(v1x-v3x)]/np.linalg.norm([-(v1y-v3y), v1x-v3x])
            Jeb[0] = np.linalg.norm([(v2x-v1x),(v2y-v1y)])/2
            Jeb[1] = np.linalg.norm([(v3x-v2x),(v3y-v2y)])/2
            Jeb[2] = np.linalg.norm([(v1x-v3x),(v1y-v3y)])/2
            for j in range(dim+1):
                [FnP, FnN] = fluxN_matrix(normal[j,0], normal[j,1])
                Fn[i,j,0,:,:] = np.round(Jeb[j]/Je*FnP,10)
                Fn[i,j,1,:,:] = np.round(Jeb[j]/Je*FnN,10)
        return Fn

elif (dim==3):

    def flux(elements, u):   ## Comunicate the solution at the border of the neighbourging elements 
        (indexB) = nodes_I3D(PolyDeg)
        # if (label=='line2'):
        # elements[label]["neigh_element_type"] = np.chararray((len(tags),1), itemsize=10)
        # elements[label]["neigh_element_type"][:] = 'triangle3'
        # elements[label]["neigh_element"] = np.zeros_like(tags, dtype=int)
        # elements[label]["neigh_side"] = np.zeros_like(tags, dtype=int)
        # else:
        # elements[label]["neigh_element_type"] = np.chararray((len(tags),3), itemsize=10)
        # elements[label]["neigh_element_type"][:] = 'triangle3'
        # elements[label]["neigh_element"] = np.zeros((len(tags),3), dtype=int)
        # elements[label]["neigh_side"] = np.zeros((len(tags),3), dtype=int)

        
        # u: Ne x N x 4
        # uext: Ne x (dim+1) x Nb x 4
        n_1=len(elements['tetrahedron4']['nodes'])
        uext = np.zeros((n_1,dim+1,2,N,dim+2))
        for i in range(n_1):
            for j in range(3):
                # uext[i,j,k,:] = u[i,j,k,:]
                typeelext = elements['tetrahedron4']["neigh_element_type"][i][j].decode("utf-8")
                if (typeelext=='tetrahedron4'):
                    indx=indexB[abs(elements['tetrahedron4']["neigh_side"][i][j])-1]
                    ## uext = [solution in the element i,  solution on the neighbourg element]
                    uext[i,j,0,indexB[j],:] = u[i,indexB[j],:]                                         
                    uext[i,j,1,indexB[j],:] = u[elements['tetrahedron4']["neigh_element"][i][j],indx[::int(sign(elements['tetrahedron4']["neigh_side"][i][j]))],:]

                else:
                    if (conditions.boundary=='NRBC'):             ##non-reflecting boundary condition
                        uext[i,j,0,indexB[j],:] = u[i,indexB[j],:]                                         
                        uext[i,j,1,indexB[j],:] = 0*u[i,indexB[j],:] 
                    elif (conditions.boundary=='wall'):           ## Wall boundary condition
                        uext[i,j,0,indexB[j],0] = u[i,indexB[j],0]  
                        uext[i,j,0,indexB[j],-1] = u[i,indexB[j],-1]  
                        uext[i,j,1,indexB[j],:] = uext[i,j,0,indexB[j],:]
        return (uext)

    def flux_matrices():
        Fx = np.array([[                                   0,                          1,                          0,                          0,                                                  0],
                    [               -Baseflow.u0*Baseflow.u0,              2*Baseflow.u0,                          0,                          0, Baseflow.c0*Baseflow.c0*Baseflow.rho0/Baseflow.pi0],
                    [               -Baseflow.u0*Baseflow.v0,                Baseflow.v0,                Baseflow.u0,                          0,                                                  0],
                    [               -Baseflow.u0*Baseflow.w0,                Baseflow.w0,                          0,                Baseflow.u0,                                                  0],
                    [-Baseflow.pi0*Baseflow.u0/Baseflow.rho0, Baseflow.pi0/Baseflow.rho0,                          0,                          0,                                        Baseflow.u0]])
        
        Fy = np.array([[                                   0,                          0,                          1,                          0,                                                  0],
                    [               -Baseflow.v0*Baseflow.u0,                Baseflow.v0,                Baseflow.u0,                          0,                                                  0],
                    [               -Baseflow.v0*Baseflow.v0,                          0,              2*Baseflow.v0,                          0, Baseflow.c0*Baseflow.c0*Baseflow.rho0/Baseflow.pi0],
                    [               -Baseflow.v0*Baseflow.w0,                          0,                Baseflow.w0,                Baseflow.v0,                                                  0],
                    [-Baseflow.pi0*Baseflow.v0/Baseflow.rho0,                          0, Baseflow.pi0/Baseflow.rho0,                          0,                                        Baseflow.v0]])
        
        Fz = np.array([[                                   0,                          0,                          0,                          1,                                                  0],
                    [               -Baseflow.u0*Baseflow.u0,                Baseflow.w0,                          0,                Baseflow.u0,                                                  0],
                    [               -Baseflow.u0*Baseflow.v0,                          0,                Baseflow.w0,                Baseflow.v0,                                                  0],
                    [               -Baseflow.u0*Baseflow.w0,                          0,                          0,              2*Baseflow.w0, Baseflow.c0*Baseflow.c0*Baseflow.rho0/Baseflow.pi0],
                    [-Baseflow.pi0*Baseflow.w0/Baseflow.rho0,                          0,                          0, Baseflow.pi0/Baseflow.rho0,                                        Baseflow.w0]])
        return (Fx,Fy,Fz)    
    
    def fluxN_matrix(nx, ny, nz):
        un = Baseflow.u0*nx+Baseflow.v0*ny+Baseflow.w0*nz
        W = np.array([[               1,      1,     0,     0,                          1],
            [Baseflow.u0-nx*Baseflow.c0,  un*nx,     0,   -nz, Baseflow.u0+nx*Baseflow.c0],
            [Baseflow.v0-ny*Baseflow.c0,  un*ny,   -nz,     0, Baseflow.v0+ny*Baseflow.c0],
            [Baseflow.w0-nz*Baseflow.c0,  un*nz,    ny,    nx, Baseflow.w0+nz*Baseflow.c0],
            [Baseflow.pi0/Baseflow.rho0,      0,     0,     0, Baseflow.pi0/Baseflow.rho0]])
        try:
            np.linalg.inv(W)
        except:
            if(abs(nx)==1):
                W = np.array([[           1,      1,     0,     0,                          1],
                [Baseflow.u0-nx*Baseflow.c0,  un*nx,   -ny,   -nz, Baseflow.u0+nx*Baseflow.c0],
                [Baseflow.v0-ny*Baseflow.c0,  un*ny,    nx,     0, Baseflow.v0+ny*Baseflow.c0],
                [Baseflow.w0-nz*Baseflow.c0,  un*nz,     0,    nx, Baseflow.w0+nz*Baseflow.c0],
                [Baseflow.pi0/Baseflow.rho0,      0,     0,     0, Baseflow.pi0/Baseflow.rho0]])
            elif(abs(ny)==1):
                W = np.array([[           1,      1,     0,     0,                          1],
                [Baseflow.u0-nx*Baseflow.c0,  un*nx,     0,   -ny, Baseflow.u0+nx*Baseflow.c0],
                [Baseflow.v0-ny*Baseflow.c0,  un*ny,   -nz,    nx, Baseflow.v0+ny*Baseflow.c0],
                [Baseflow.w0-nz*Baseflow.c0,  un*nz,    ny,     0, Baseflow.w0+nz*Baseflow.c0],
                [Baseflow.pi0/Baseflow.rho0,      0,     0,     0, Baseflow.pi0/Baseflow.rho0]])
            else:
                print('WARNING: Decomposition of flux matrix not considered')
        if(conditions.flux=="center"):
            D =  np.diag(np.array([un-Baseflow.c0, un, un, un, un+Baseflow.c0]))
            FnP = np.dot(np.dot(W, D), np.linalg.inv(W))/2  ## so Fn original is the sum of FnP + FnN
            FnN = np.dot(np.dot(W, D), np.linalg.inv(W))/2
        if(conditions.flux=="upwind"):
            eigenvalP = np.array([un-Baseflow.c0, un, un, un, un+Baseflow.c0])
            eigenvalN = np.array([un-Baseflow.c0, un, un, un, un+Baseflow.c0])
            eigenvalP[eigenvalP<0]=0
            eigenvalN[eigenvalN>0]=0
            DP =  np.diag(eigenvalP)
            DN =  np.diag(eigenvalN)
            FnP = np.dot(np.dot(W, DP), np.linalg.inv(W)) # W*D*inv(W)
            FnN = np.dot(np.dot(W, DN), np.linalg.inv(W))
            # Fn = np.array([[                             0,                            nx,                              ny,                              nz,                                                       0],
            #                [               -Baseflow.u0*un,             un+Baseflow.u0*nx,                  Baseflow.u0*ny,                  Baseflow.u0*nz,   Baseflow.rho0*Baseflow.c0*Baseflow.c0*nx/Baseflow.pi0],
            #                [               -Baseflow.v0*un,                Baseflow.v0*nx,               un+Baseflow.v0*ny,                  Baseflow.v0*nz,   Baseflow.rho0*Baseflow.c0*Baseflow.c0*ny/Baseflow.pi0],
            #                [               -Baseflow.w0*un,                Baseflow.w0*nx,                  Baseflow.w0*ny,               un+Baseflow.w0*nz,   Baseflow.rho0*Baseflow.c0*Baseflow.c0*nz/Baseflow.pi0],
            #                [-Baseflow.pi0*un/Baseflow.rho0, Baseflow.pi0*nx/Baseflow.rho0,   Baseflow.pi0*ny/Baseflow.rho0,   Baseflow.pi0*nz/Baseflow.rho0,                                                      un]])
        else:
            print('case not available in fluxN_matrix()')
        return (FnP, FnN)

    def decomposed_FN_matrix(nodes, elements, Ne):
        Jeb = np.zeros((dim+1))
        normal = np.zeros((dim+1,dim))
        Fn = np.zeros((Ne,dim+1,2,5,5)) ###(Number of elements, number of faces, Fn positive or Fn negative, (dim Fn) )
        for i in range(Ne):
            ## positions of the elements' vertices
            v1x, v1y, v1z = nodes["coordinates"][int(elements['tetrahedron4']["nodes"][i][0]-1)]
            v2x, v2y, v2z = nodes["coordinates"][int(elements['tetrahedron4']["nodes"][i][1]-1)]
            v3x, v3y, v3z = nodes["coordinates"][int(elements['tetrahedron4']["nodes"][i][2]-1)]
            v4x, v4y, v4z = nodes["coordinates"][int(elements['tetrahedron4']["nodes"][i][3]-1)]
            ## Mtemp = [[dxdr, dxds, dxdt],
            ##          [dydr, dyds, dydt],
            ##          [dzdr, dzds, dzdt]]
            Mtemp = [[(v2x-v1x)/2, (-v3x-v1x)/2, (v4x-v1x)/2],[(v2y-v1y)/2, (-v3y-v1y)/2, (v4y-v1y)/2],[(v2z-v1z)/2, (-v3z-v1z)/2, (v4z-v1z)/2]]
            Je = np.linalg.det(Mtemp)
            normal[0,:] = [(v2y-v1y)*(v4z-v1z)-(v2z-v1z)*(v4y-v1y),(v2z-v1z)*(v4x-v1x)-(v2x-v1x)*(v4z-v1z),(v2x-v1x)*(v4y-v1y)-(v2y-v1y)*(v4x-v1x)]/np.linalg.norm([(v2y-v1y)*(v4z-v1z)-(v2z-v1z)*(v4y-v1y),(v2z-v1z)*(v4x-v1x)-(v2x-v1x)*(v4z-v1z),(v2x-v1x)*(v4y-v1y)-(v2y-v1y)*(v4x-v1x)])
            normal[1,:] = [(v3y-v2y)*(v4z-v2z)-(v3z-v2z)*(v4y-v2y),(v3z-v2z)*(v4x-v2x)-(v3x-v2x)*(v4z-v2z),(v3x-v2x)*(v4y-v2y)-(v3y-v2y)*(v4x-v2x)]/np.linalg.norm([(v3y-v2y)*(v4z-v2z)-(v3z-v2z)*(v4y-v2y),(v3z-v2z)*(v4x-v2x)-(v3x-v2x)*(v4z-v2z),(v3x-v2x)*(v4y-v2y)-(v3y-v2y)*(v4x-v2x)])
            normal[2,:] = [(v1y-v3y)*(v4z-v3z)-(v1z-v3z)*(v4y-v3y),(v1z-v3z)*(v4x-v3x)-(v1x-v3x)*(v4z-v3z),(v1x-v3x)*(v4y-v3y)-(v1y-v3y)*(v4x-v3x)]/np.linalg.norm([(v1y-v3y)*(v4z-v3z)-(v1z-v3z)*(v4y-v3y),(v1z-v3z)*(v4x-v3x)-(v1x-v3x)*(v4z-v3z),(v1x-v3x)*(v4y-v3y)-(v1y-v3y)*(v4x-v3x)])
            normal[3,:] = [(v2y-v1y)*(v3z-v1z)-(v2z-v1z)*(v3y-v1y),(v2z-v1z)*(v3x-v1x)-(v2x-v1x)*(v3z-v1z),(v2x-v1x)*(v3y-v1y)-(v2y-v1y)*(v3x-v1x)]/np.linalg.norm([(v2y-v1y)*(v3z-v1z)-(v2z-v1z)*(v3y-v1y),(v2z-v1z)*(v3x-v1x)-(v2x-v1x)*(v3z-v1z),(v2x-v1x)*(v3y-v1y)-(v2y-v1y)*(v3x-v1x)])
            Jeb[0] = np.linalg.det([[(v2x-v1x)/2, (v4x-v1x)/2],[(v2y-v1y)/2, (v4y-v1y)/2]])
            Jeb[1] = np.linalg.det([[(v3x-v2x)/2, (v4x-v2x)/2],[(v3y-v2y)/2, (v4y-v2y)/2]])
            Jeb[2] = np.linalg.det([[(v1x-v3x)/2, (v4x-v3x)/2],[(v1y-v3y)/2, (v4y-v3y)/2]])
            Jeb[3] = np.linalg.det([[(v2x-v1x)/2, (v3x-v1x)/2],[(v2y-v1y)/2, (v3y-v1y)/2]])
            for j in range(dim+1):
                [FnP, FnN] = fluxN_matrix(normal[j,0], normal[j,1], normal[j,2])
                Fn[i,j,0,:,:] = np.round(Jeb[j]/Je*FnP,10)
                Fn[i,j,1,:,:] = np.round(Jeb[j]/Je*FnN,10)
        return Fn

else:
    print('error: dimension not available')

