from math import sqrt
import numpy as np 
from parameters import CFL, num, Baseflow
from mesh import Nodes2D, Nodes3D

dim = num.dim
if (dim==2):

    def dtscale(elements, nodes):

    # function dtscale = dtscale2D;
    # Purpose : Compute inscribed circle diameter as characteristic
    #           for grid to choose timestep


        #Find vertex nodes
        vx = nodes["coordinates"][(elements['triangle3']["nodes"]-1).astype(int)][:,:,0]
        vy = nodes["coordinates"][(elements['triangle3']["nodes"]-1).astype(int)][:,:,1]

        #Compute semi-perimeter and area
        len1 = np.sqrt((vx[:,0]-vx[:,1])**2+(vy[:,0]-vy[:,1])**2)
        len2 = np.sqrt((vx[:,1]-vx[:,2])**2+(vy[:,1]-vy[:,2])**2)
        len3 = np.sqrt((vx[:,2]-vx[:,0])**2+(vy[:,2]-vy[:,0])**2)
        sper = (len1 + len2 + len3)/2.0
        Area = np.sqrt(sper*(sper-len1)*(sper-len2)*(sper-len3))

        #Compute scale using radius of inscribed circle
        xI,_ = Nodes2D(num.PolyDeg)
        dxi = xI[1]-xI[0]
        dtscale = 2/CFL.cfl*dxi*min(Area/sper)/(4*Baseflow.c0+sqrt(Baseflow.u0**2+Baseflow.v0**2)) #CFL*min(dxi)*min(rD)/v

        return (dtscale)
elif(dim==3):

    def dtscale(elements, nodes):

    # function dtscale = dtscale2D;
    # Purpose : Compute inscribed circle diameter as characteristic
    #           for grid to choose timestep

        #Find vertex nodes
        vx = nodes["coordinates"][(elements['tetrahedron4']["nodes"]-1).astype(int)][:,:,0]
        vy = nodes["coordinates"][(elements['tetrahedron4']["nodes"]-1).astype(int)][:,:,1]
        vz = nodes["coordinates"][(elements['tetrahedron4']["nodes"]-1).astype(int)][:,:,2]

        #Compute semi-perimeter and area
        v01 = [vx[:,0]-vx[:,1],vy[:,0]-vy[:,1],vz[:,0]-vz[:,1]]
        v02 = [vx[:,0]-vx[:,2],vy[:,0]-vy[:,2],vz[:,0]-vz[:,2]]
        v03 = [vx[:,0]-vx[:,3],vy[:,0]-vy[:,3],vz[:,0]-vz[:,3]]
        v12 = [vx[:,1]-vx[:,2],vy[:,1]-vy[:,2],vz[:,1]-vz[:,2]]
        v13 = [vx[:,1]-vx[:,3],vy[:,1]-vy[:,3],vz[:,1]-vz[:,3]]
        Vol = abs(np.vecdot(np.transpose(v01), np.transpose(np.cross(v02, v03, axis=0)))) / 6

        Area = 0.5 *(np.linalg.norm(np.cross(v01, v02, axis=0),axis=0)+np.linalg.norm(np.cross(v01, v03, axis=0),axis=0)+np.linalg.norm(np.cross(v02, v03, axis=0),axis=0)+np.linalg.norm(np.cross(v12, v13, axis=0),axis=0))
        
        #Compute scale using radius of inscribed circle
        xI,_,_ = Nodes3D(num.PolyDeg)
        dxi = xI[1]-xI[0]
        dtscale = float(CFL.cfl*dxi*min(3*Vol/Area)/(Baseflow.c0+sqrt(Baseflow.u0**2+Baseflow.v0**2+Baseflow.w0**2))) #CFL*min(dxi)*min(radius of the inscribed circunference)/v
        return (dtscale)   