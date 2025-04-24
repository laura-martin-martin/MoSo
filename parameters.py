
class num:
    dim = 2                     # define it also in matrices.py
    PolyDeg = 5                # Polinomial degree per dimension
    Nb = PolyDeg+1              # Number of nodes per face
    N = int(Nb*(Nb+1)/2) 


class Baseflow:
    gamma = 1.4
    u0 = 0. 
    v0 = 0. 
    w0 = 0.
    c0 = 1. #343. 
    rho0 = 1. #1.2 
    p0 = rho0*c0*c0/gamma
    pinf = p0  # otherwise mean value of p0 over the domain
    piinf = 1.
    pi0 = piinf*(p0/pinf)**(1/gamma)

class CFL:
    cfl = 0.5

class times:
    FinalTime = 0.0
    TotalSteps = 801

class initial_conds:
    activated = 'False'
    type = 'pulse'
    positionX = 0.5
    positionY = -0.25
    width = 0.2
    amplitude = 1

class sources:
    activated = 'True'
    type = 'monopole'
    freq = 3
    positionX = 0.5
    positionY = -0.25
    width = 0.15
    amplitude = 5

class conditions:
    boundary = 'NRBC'   # 'wall'
    flux = 'upwind'     # 'center'

class mesh:
    name_file = "computation_fine.msh" # "computation.msh" # 
