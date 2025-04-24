import numpy as np
from parameters import Baseflow, initial_conds

baseflow = Baseflow()
InitCond = initial_conds()

def init_sol(Ne, X, Y, u):
    X0 = InitCond.positionX
    Y0 = InitCond.positionY
    if(InitCond.type=='pulse'):
        for i in range(Ne):
            u[i,:,0] = InitCond.amplitude*np.exp(-((X[i,:]-X0)*(X[i,:]-X0)+(Y[i,:]-Y0)*(Y[i,:]-Y0))/InitCond.width/InitCond.width)
            u[i,:,1] = 0
            u[i,:,2] = 0
            u[i,:,3] = baseflow.pi0*u[i,:,0]/baseflow.rho0
    return(u)