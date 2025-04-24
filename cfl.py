from math import sqrt
import numpy as np 
from parameters import CFL, num, Baseflow
from mesh import Nodes2D


def dtscale2D(elements, nodes):

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