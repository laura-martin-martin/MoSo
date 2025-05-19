import numpy as np
import matplotlib.pyplot as plt
import gmsh
import math
import scipy.io
from parameters import num
from scipy.special import gamma
from flux import flux_matrices
from mesh import nodes_I2D, Nodes2D, xytors, rstoab, Nodes3D, xyztorst, rsttoabc, nodes_I3D

dim = num.dim
Nfaces = dim+1
PolyDeg = num.PolyDeg
Nb = num.Nb
N = num.N

################## Vandermonde Matrix #####################
def JacobiP(x,alpha,beta,N):

        # function [P] = JacobiP(x,alpha,beta,N)
        # Purpose: Evaluate Jacobi Polynomial of type (alpha,beta) > -1
        #          (alpha+beta <> -1) at points x for order N and returns P[1:length(xp))]
        # Note   : They are normalized to be orthonormal.

        # Turn points into row if needed.
        #xp = x
        ###dims = size(xp)
        ###if (dims(2)==1) xp = xp'
        #end

        PL = np.zeros((N+1,len(x)))

        # Initial values P_0(x) and P_1(x)
        gamma0 = 2**(alpha+beta+1)/(alpha+beta+1)*gamma(alpha+1)*gamma(beta+1)/gamma(alpha+beta+1)
        PL[0,:] = 1.0/math.sqrt(gamma0)
        if (N==0):
            P=np.transpose(PL)
        else:
            gamma1 = (alpha+1)*(beta+1)/(alpha+beta+3)*gamma0
            PL[1,:] = np.squeeze(((alpha+beta+2)*x/2 + (alpha-beta)/2)/math.sqrt(gamma1))
            if (N==1):
                P=np.transpose(PL[N,:])
            else:
                # Repeat value in recurrence.
                aold = 2/(2+alpha+beta)*math.sqrt((alpha+1)*(beta+1)/(alpha+beta+3))

                # Forward recurrence using the symmetry of the recurrence.
                for i in range(1,N):
                    h1 = 2*i+alpha+beta
                    anew = 2/(h1+2)*math.sqrt( (i+1)*(i+1+alpha+beta)*(i+1+alpha)*(i+1+beta)/(h1+1)/(h1+3))
                    bnew = - (alpha**2-beta**2)/h1/(h1+2)
                    PL[i+1,:] = (1/anew*( -aold*PL[i-1,:] + (np.squeeze(x)-bnew)*PL[i,:]))
                    aold =anew

                P = np.transpose(PL[N,:])

        if (len(np.shape(P))==1):
            P = P[:,np.newaxis]
        return P

################# Differentiation matrices ###################
def GradJacobiP(r, alpha, beta, N):

        # function [dP] = GradJacobiP(r, alpha, beta, N)

        # Purpose: Evaluate the derivative of the Jacobi polynomial of type (alpha,beta)>-1,
        #	       at points r for order N and returns dP[1:length(r))]        

        dP = np.zeros((len(r), 1))
        if(N == 0):
            dP[:,:] = 0.0
        else:
            dP = math.sqrt(N*(N+alpha+beta+1))*JacobiP(r,alpha+1,beta+1, N-1)

        return dP

############# function to generate matrices at the elements borders ##################
  
def Simplex2DP(a,b,i,j):

    # function [P] = Simplex2DP(a,b,i,j)

    # Purpose : Evaluate 2D orthonormal polynomial
    #           on simplex at (a,b) of order (i,j).

    h1 = JacobiP(a,0,0,i)
    h2 = JacobiP(b,2*i+1,0,j)
    P = math.sqrt(2.0)*np.squeeze(h1*h2)*np.pow((1-b),i)
    return P

def Simplex3DP(a,b,c,i,j,k):


    # function [P] = Simplex3DP(a,b,c,i,j,k)

    # Purpose : Evaluate 3D orthonormal polynomial
    #           on simplex at (a,b,c) of order (i,j,k).

    h1 = JacobiP(a,0,0,i)
    h2 = JacobiP(b,2*i+1,0,j)
    h3 = JacobiP(c,2*(i+j)+2,0,k)

    P = 2*math.sqrt(2)*h1*h2*((1-b)**i)*h3*((1-c)**(i+j))
    return (P)

def Vandermonde1D(N,r):

    # function [V1D] = Vandermonde1D(N,r)
    # Purpose : Initialize the 1D Vandermonde Matrix, V_{ij} = phi_j(r_i)


    V1D = np.zeros((len(r),N+1))
    for j in range(N+1):
        V1D[:,j] = np.squeeze(JacobiP(r, 0, 0, j))
    
    return V1D

    ################## Vandermonde Matrix #####################

def Vandermonde2D(PolyDeg, r, s):

    # # function [V2D] = Vandermonde2D(N, r, s)

    # # Purpose : Initialize the 2D Vandermonde Matrix,  V_{ij} = phi_j(r_i, s_i)


    # V2D = np.zeros((len(r),int((N+1)*(N+2)/2)))

    # # Transfer to (a,b) coordinates
    # [a, b] = rstoab(r, s)
    # # build the Vandermonde matrix
    # sk = 0
    # for i in range(N+1):
    #     for j in range(N-i+1):
    #         V2D[:,sk] = np.squeeze(Simplex2DP(a,b,i,j))
    #         sk = sk+1
    if (PolyDeg>0 & PolyDeg<11):
        data = scipy.io.loadmat('Vander2D.mat')
        V = data['V'][0][PolyDeg]

    return V

def Vandermonde3D(PolyDeg, r, s, t):

    # # function [V3D] = Vandermonde3D(N, r, s, t)
    # # Purpose : Initialize the 3D Vandermonde Matrix, V_{ij} = phi_j(r_i, s_i, t_i)

    # V3D = np.zeros((len(r),(N+1)*(N+2)*(N+3)/6))


    # # Transfer to (a,b) coordinates
    # [a, b, c] = rsttoabc(r, s, t)

    # # build the Vandermonde matrix
    # sk = 0

    # for i in range(N+1):       # old ordering
    #     for j in range(N - i+1):
    #         for k in range(N - i - j+1):
    #             V3D[:,sk] = Simplex3DP(a,b,c,i,j,k)
    #             sk = sk+1
    if (PolyDeg>0 & PolyDeg<11):
        data = scipy.io.loadmat('Vander3D.mat')
        V3D = data['V'][0][PolyDeg]

    return (V3D)

################# Differentiation matrices ###################

def GradSimplex2DP(a,b,id,jd):

    # function [dmodedr, dmodeds] = GradSimplex2DP(a,b,id,jd)
    # Purpose: Return the derivatives of the modal basis (id,jd) on the 2D simplex at (a,b).   

    fa = JacobiP(a, 0, 0, id)
    dfa = GradJacobiP(a, 0, 0, id)
    gb = JacobiP(b, 2*id+1,0, jd)
    dgb = GradJacobiP(b, 2*id+1,0, jd)

    # r-derivative
    # d/dr = da/dr d/da + db/dr d/db = (2/(1-s)) d/da = (2/(1-b)) d/da
    dmodedr = dfa*gb
    if(id>0):
        dmodedr = dmodedr*((0.5*(1-b))**(id-1))

    # s-derivative
    # d/ds = ((1+a)/2)/((1-b)/2) d/da + d/db
    dmodeds = dfa*(gb*(0.5*(1+a)))
    if(id>0):
        dmodeds = dmodeds*((0.5*(1-b))**(id-1))

    tmp = dgb*((0.5*(1-b))**id)
    if(id>0):
        tmp = tmp-0.5*id*gb*((0.5*(1-b))**(id-1))
    
    dmodeds = dmodeds+fa*tmp

    # Normalize
    dmodedr = 2**(id+0.5)*dmodedr
    dmodeds = 2**(id+0.5)*dmodeds

    return (dmodedr, dmodeds)

def GradSimplex3DP(a,b,c,id,jd,kd):

    # function [V3Dr, V3Ds, V3Dt] = GradSimplex3DP(a,b,c,id,jd,kd)
    # Purpose: Return the derivatives of the modal basis (id,jd,kd) 
    #          on the 3D simplex at (a,b,c)

    fa = JacobiP(a,0,0,id)
    dfa = GradJacobiP(a,0,0,id)

    gb = JacobiP(b,2*id+1,0,jd)
    dgb = GradJacobiP(b,2*id+1,0,jd)

    hc = JacobiP(c,2*(id+jd)+2,0,kd)
    dhc = GradJacobiP(c,2*(id+jd)+2,0,kd)


    # r-derivative
    V3Dr = dfa*(gb*hc)

    if(id>0):
        V3Dr = V3Dr*((0.5*(1-b))**(id-1))
    if(id+jd>0):
        V3Dr = V3Dr*((0.5*(1-c))**(id+jd-1))

    # s-derivative 
    V3Ds = 0.5*(1+a)*V3Dr

    tmp = dgb*((0.5*(1-b))**id)

    if(id>0):
        tmp = tmp+(-0.5*id)*(gb*(0.5*(1-b))**(id-1))
    if(id+jd>0):
        tmp = tmp*((0.5*(1-c))**(id+jd-1))
    tmp = fa*(tmp*hc)

    V3Ds = V3Ds+tmp


    # t-derivative 
    V3Dt = 0.5*(1+a)*V3Dr+0.5*(1+b)*tmp

    tmp = dhc*((0.5*(1-c))**(id+jd))

    if(id+jd>0):
        tmp = tmp-0.5*(id+jd)*(hc*((0.5*(1-c))**(id+jd-1)))

    tmp = fa*(gb*tmp)
    tmp = tmp*((0.5*(1-b))**id)

    V3Dt = V3Dt+tmp

    # normalize
    V3Dr = V3Dr*(2**(2*id+jd+1.5))

    V3Ds = V3Ds*(2**(2*id+jd+1.5))

    V3Dt = V3Dt*(2**(2*id+jd+1.5))

    return (V3Dr, V3Ds, V3Dt)

def GradVandermonde2D(N,r,s):

    # function [V2Dr,V2Ds] = GradVandermonde2D(N,r,s)
    # Purpose : Initialize the gradient of the modal basis (i,j) at (r,s) at order N	

    V2Dr = np.zeros((len(r),int((N+1)*(N+2)/2)))
    V2Ds = np.zeros((len(r),int((N+1)*(N+2)/2)))

    # find tensor-product coordinates
    [a,b] = rstoab(r,s)

    # Initialize matrices
    sk = 0
    for i in range(N+1):
        for j in range(N-i+1):
            [temp1,temp2] = GradSimplex2DP(a,b,i,j)
            V2Dr[:,sk] = np.squeeze(temp1)
            V2Ds[:,sk] = np.squeeze(temp2) 
            sk = sk+1
        
    return (V2Dr,V2Ds)

def GradVandermonde3D(N,r,s,t):

    # function [V3Dr,V3Ds,V3Dt] = GradVandermonde3D(N,r,s,t)
    # Purpose : Initialize the gradient of the modal basis (i,j,k) at (r,s,t) at order N

    V3Dr = np.zeros((len(r),int((N+1)*(N+2)*(N+3)/6)))
    V3Ds = np.zeros((len(r),int((N+1)*(N+2)*(N+3)/6)))
    V3Dt = np.zeros((len(r),int((N+1)*(N+2)*(N+3)/6)))


    # find tensor-product coordinates
    [a,b,c] = rsttoabc(r,s,t)


    # Initialize matrices

    sk = 0

    for i in range(N+1):
        for j in range(N-i+1):
            for k in range(N-i-j+1):
                [temp1, temp2, temp3] = GradSimplex3DP(a,b,c,i,j,k)
                V3Dr[:,sk] = np.squeeze(temp1)
                V3Ds[:,sk] = np.squeeze(temp2) 
                V3Dt[:,sk] = np.squeeze(temp3)
                sk = sk+1
    return (V3Dr,V3Ds,V3Dt)

def Dmatrices2D(N,r,s,invV):

    # function [Dr,Ds] = Dmatrices2D(N,r,s,V)
    # Purpose : Initialize the (r,s) differentiation matrices
    #	    on the simplex, evaluated at (r,s) at order N

    [Vr, Vs] = GradVandermonde2D(N, r, s)
    Dr = np.dot(Vr,invV)
    Ds = np.dot(Vs,invV)
    return (Dr,Ds)

def Dmatrices3D(N,r,s,t,invV):

    # function [Dr,Ds,Dt] = Dmatrices3D(N,r,s,t,V)
    # Purpose : Initialize the (r,s,t) differentiation matrices
    #	        on the simplex, evaluated at (r,s,t) at order N

    [Vr, Vs, Vt] = GradVandermonde3D(N, r, s, t)

    Dr = np.dot(Vr,invV)
    Ds = np.dot(Vs,invV)
    Dt = np.dot(Vt,invV)

    return (Dr,Ds,Dt)

############# function to generate matrices at the elements borders ##################

def Lift2D(V, r, s):

    # function [LIFT] = Lift2D()
    # Purpose  : Compute surface to volume lift term for DG formulation
    [F1,F2,F3] = nodes_I2D(PolyDeg)

    Emat = np.zeros((N, Nfaces*Nb))
    
    # face 1
    faceR = r[F1]   ## coordinates
    V1D = Vandermonde1D(PolyDeg, faceR)
    massEdge1 = np.linalg.inv(np.dot(V1D,np.transpose(V1D)))
    Emat[F1,0:Nb] = massEdge1

    # face 2
    faceR = r[F2] ## coordinates
    V1D = Vandermonde1D(PolyDeg, faceR)
    massEdge2 = np.linalg.inv(np.dot(V1D,np.transpose(V1D)))
    Emat[F2,Nb:2*Nb] = massEdge2

    # face 3
    faceS = s[F3] ## coordinates
    V1D = Vandermonde1D(PolyDeg, faceS)
    massEdge3 = np.linalg.inv(np.dot(V1D,np.transpose(V1D)))
    Emat[F3,2*Nb:] = massEdge3
    
    # # inv(mass matrix)*\I_n (L_i,L_j)_{edge_n}
    LIFT = np.matmul(V,(np.matmul(np.transpose(V),Emat)))
    return LIFT

def Lift3D(V,R,S,T):

    # function [LIFT] = Lift3D(N, r, s, t)
    # Purpose  : Compute 3D surface to volume lift operator used in DG formulation


    Emat = np.zeros((N, Nfaces*Nb))

    [F1,F2,F3, F4] = nodes_I3D(PolyDeg)

    for face in range(Nfaces):
        # process face
        if(face==0):
            faceR = R[F1]
            faceS = S[F1]
            idr = np.array(F1)

        if(face==1):
            faceR = R[F2]
            faceS = T[F2]
            idr = np.array(F2)

        if(face==2):
            faceR = S[F3]
            faceS = T[F3]
            idr = np.array(F3)

        if(face==3):
            faceR = S[F4]
            faceS = T[F4]
            idr = np.array(F4)

    
        VFace = Vandermonde2D(PolyDeg, faceR, faceS)

        massFace = np.linalg.inv(np.dot(VFace,np.transpose(VFace)))

        idc = np.arange((face)*Nb,(face+1)*Nb)    
        Emat[np.ix_(idr, idc)] = Emat[np.ix_(idr, idc)] + massFace

    # inv(mass matrix)*\I_n (L_i,L_j)_{edge_n}
    LIFT = np.matmul(V,(np.matmul(np.transpose(V),Emat)))

    return (LIFT)

####################################################################################################################################################################################################################################################
################################################################################################################# Dim=2 ############################################################################################################################


    ############## funtion to generate all the mtrices at the same time ##################

def matrices(nodes, elements, Ne):
        
    if (dim==2):
        np.set_printoptions(linewidth=300)
        Gx = np.zeros((Ne,N,N))      #Mref^(-1)*(drdx*Gr+dsdx*Gs+drdy*Gr+dsdy*Gs), Gr=int_e dPHIdr*PHI, Gs=int_e dPHIds*PHI 
        Gy = np.zeros((Ne,N,N))
        GF = np.zeros((Ne,N,N*4, 4))
        drdx = np.zeros((dim,dim))  ######### (2,2) ?
        x,y = Nodes2D(PolyDeg)
        r,s = xytors(x,y)
        # Build reference element matrices 
        V = Vandermonde2D(PolyDeg, r, s)
        # # # invV = np.linalg.inv(V)
        # # # [Dr,Ds] = Dmatrices2D(PolyDeg, r, s, invV)
        # # # Mref = np.dot(np.transpose(invV),invV)
        Bref = Lift2D(V, r, s)  ## dim = (Np, Nfaces, Nb)
        Vr, Vs = GradVandermonde2D(PolyDeg, r, s)
        Drw = np.matmul(np.matmul(V,np.transpose(Vr)),np.linalg.inv(np.matmul(V,np.transpose(V))))
        Dsw = np.matmul(np.matmul(V,np.transpose(Vs)),np.linalg.inv(np.matmul(V,np.transpose(V))))
        for i in range(Ne):
            ## positions of the elements' vertices
            v1x, v1y, _ = nodes["coordinates"][int(elements['triangle3']["nodes"][i][0]-1)]
            v2x, v2y, _ = nodes["coordinates"][int(elements['triangle3']["nodes"][i][1]-1)]
            v3x, v3y, _ = nodes["coordinates"][int(elements['triangle3']["nodes"][i][2]-1)]
            ## Mtemp = [[dxdr, dxds,
            ##          [dydr, dyds]]
            Mtemp = [[(v2x-v1x)/2, (v3x-v1x)/2],[(v2y-v1y)/2, (v3y-v1y)/2]]
            ## inv(Mtemp) = [[drdx, drdy],
            ##               [dsdx, dsdy]]
            drdx[:,:] = np.linalg.inv(Mtemp) 
            # G = Gx*Fx + Gy*Fy = inv(Mref)*[(drdx*Dr+dsdx*Ds)*Fx^t+(drdy*Dr+dsdy*Ds)*Fy^t]
            Gx[i,:,:] = drdx[0,0]*Drw+drdx[1,0]*Dsw #np.dot(np.linalg.inv(Mref),(drdx[0,0]*Drw+drdx[1,0]*Dsw))
            Gy[i,:,:] = drdx[0,1]*Drw+drdx[1,1]*Dsw #np.dot(np.linalg.inv(Mref),(drdx[0,1]*Drw+drdx[1,1]*Dsw))
        
        Fx, Fy = flux_matrices() # here because the baseflow is uniform, otherwise it sould be calculates for each virtual node
        for i in range(4):
            for j in range(4):
                GF[:,:,N*i:N*(i+1),j]=Gx*Fx[j,i]+Gy*Fy[j,i]
    
    elif(dim==3):
    
        np.set_printoptions(linewidth=300)
        Gx = np.zeros((Ne,N,N))      #Mref^(-1)*(drdx*Gr+dsdx*Gs+drdy*Gr+dsdy*Gs), Gr=int_e dPHIdr*PHI, Gs=int_e dPHIds*PHI 
        Gy = np.zeros((Ne,N,N))
        Gz = np.zeros((Ne,N,N))
        GF = np.zeros((Ne,N,N*(dim+2), dim+2))
        drdx = np.zeros((dim,dim))  ######### (2,2) ?
        x,y,z = Nodes3D(PolyDeg)
        r,s,t = xyztorst(x,y,z)
        # Build reference element matrices 
        V = Vandermonde3D(PolyDeg, r, s, t)
        # # # invV = np.linalg.inv(V)
        # # # [Dr,Ds] = Dmatrices2D(PolyDeg, r, s, invV)
        # # # Mref = np.dot(np.transpose(invV),invV)
        Bref = Lift3D(V, r, s, t)  ## dim = (Np, Nfaces, Nb)
        Vr, Vs, Vt = GradVandermonde3D(PolyDeg, r, s, t)
        Drw = np.matmul(np.matmul(V,np.transpose(Vr)),np.linalg.inv(np.matmul(V,np.transpose(V))))
        Dsw = np.matmul(np.matmul(V,np.transpose(Vs)),np.linalg.inv(np.matmul(V,np.transpose(V))))
        Dtw = np.matmul(np.matmul(V,np.transpose(Vt)),np.linalg.inv(np.matmul(V,np.transpose(V))))
        for i in range(Ne):
            ## positions of the elements' vertices
            v1x, v1y, v1z = nodes["coordinates"][int(elements['tetrahedron4']["nodes"][i][0]-1)]
            v2x, v2y, v2z = nodes["coordinates"][int(elements['tetrahedron4']["nodes"][i][1]-1)]
            v3x, v3y, v3z = nodes["coordinates"][int(elements['tetrahedron4']["nodes"][i][2]-1)]
            v4x, v4y, v4z = nodes["coordinates"][int(elements['tetrahedron4']["nodes"][i][3]-1)]
            ## Mtemp = [[dxdr, dxds, dxdt],
            ##          [dydr, dyds, dydt],
            ##          [dzdr, dzds, dzdt]]
            Mtemp = np.array([[(v2x-v1x)/2, (-v3x-v1x)/2, (v4x-v1x)/2],[(v2y-v1y)/2, (-v3y-v1y)/2, (v4y-v1y)/2],[(v2z-v1z)/2, (-v3z-v1z)/2, (v4z-v1z)/2]])
            ## inv(Mtemp) = [[drdx, drdy, drdz],
            ##               [dsdx, dsdy, dsdz],
            ##               [dtdx, dtdy, dtdz]]
            drdx[:,:] = np.linalg.inv(Mtemp) 
            # G = Gx*Fx + Gy*Fy = inv(Mref)*[(drdx*Dr+dsdx*Ds)*Fx^t+(drdy*Dr+dsdy*Ds)*Fy^t]
            Gx[i,:,:] = drdx[0,0]*Drw+drdx[1,0]*Dsw+drdx[2,0]*Dtw 
            Gy[i,:,:] = drdx[0,1]*Drw+drdx[1,1]*Dsw+drdx[2,1]*Dtw 
            Gz[i,:,:] = drdx[0,2]*Drw+drdx[1,2]*Dsw+drdx[2,2]*Dtw 
        
        Fx, Fy, Fz = flux_matrices() # here because the baseflow is uniform, otherwise it sould be calculates for each virtual node

        for i in range(dim+2):
            for j in range(dim+2):
                GF[:,:,N*i:N*(i+1),j]=Gx*Fx[j,i]+Gy*Fy[j,i]+Gz*Fz[j,i]
        # print(Gx[0,0,:])
        # print(Fx[:,4])
        # print(GF[0,0,140:,:])
        # print((Gx*Fx[3,4])[0,0,140:,:])
    else:
        print('ERROR: not correct dimension')
    return (GF, Bref)



####################################################################################################################################################################################################################################################
################################################################################################################# Dim=3 ############################################################################################################################
