import math
import numpy as np
from parameters import Baseflow,sources, num
dim = num.dim

if(dim==2):
    def source(Ne, X, Y, time):
        X0 = sources.positionX(time)
        Y0 = sources.positionY(time)
        s = np.zeros((Ne,num.N,dim+2))
        if(sources.type=='monopole'):
            for i in range(Ne):
                s[i,:,0] = sources.amplitude*np.exp(-((X[i,:]-X0)*(X[i,:]-X0)+(Y[i,:]-Y0)*(Y[i,:]-Y0))/sources.width/sources.width)*math.sin(2*math.pi*sources.freq*time)
                s[i,:,1] = 0
                s[i,:,2] = 0
                s[i,:,3] = Baseflow.pi0*s[i,:,0]/Baseflow.rho0
        return(s)
elif(dim==3):
    def source(Ne, X, Y, Z, time):
        X0 = sources.positionX(time)
        Y0 = sources.positionY(time)
        Z0 = sources.positionZ(time)
        s = np.zeros((Ne,num.N,dim+2))
        if(sources.type=='monopole'):
            for i in range(Ne):
                s[i,:,0] = sources.amplitude*np.exp(-((X[i,:]-X0)*(X[i,:]-X0)+(Y[i,:]-Y0)*(Y[i,:]-Y0)+(Z[i,:]-Z0)*(Z[i,:]-Z0))/sources.width/sources.width)*math.sin(2*math.pi*sources.freq*time)
                s[i,:,1] = 0
                s[i,:,2] = 0
                s[i,:,3] = 0
                s[i,:,4] = Baseflow.pi0*s[i,:,0]/Baseflow.rho0
        return(s)