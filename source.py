import math
import numpy as np
from parameters import Baseflow,sources, num

def source(Ne, X, Y, time):
    X0 = sources.positionX
    Y0 = sources.positionY
    s = np.zeros((Ne,num.N,4))
    if(sources.type=='monopole'):
        for i in range(Ne):
            s[i,:,0] = sources.amplitude*np.exp(-((X[i,:]-X0)*(X[i,:]-X0)+(Y[i,:]-Y0)*(Y[i,:]-Y0))/sources.width/sources.width)*math.sin(2*math.pi*sources.freq*time)
            s[i,:,1] = 0
            s[i,:,2] = 0
            s[i,:,3] = Baseflow.pi0*s[i,:,0]/Baseflow.rho0
    return(s)