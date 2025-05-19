import numpy as np
from mesh import nodes_I2D, nodes_I3D
from parameters import num

dim = num.dim
PolyDeg = num.PolyDeg
Nb = num.Nb
N = num.N

rk4a = [ 0.0, -567301805773.0/1357537059087.0, -2404267990393.0/2016746695238.0, -3550918686646.0/2091501179385.0, -1275806237668.0/842570457699.0]
rk4b = [ 1432997174477.0/9575080441755.0, 5161836677717.0/13612068292357.0, 1720146321549.0/2090206949498.0, 3134564353537.0/4481467310338.0, 2277821191437.0/14882151754819.0]
rk4c = [ 0.0, 1432997174477.0/9575080441755.0, 2526269341429.0/6820363962896.0, 2006345519317.0/3224310063776.0, 2802321613138.0/2924317926251.0]


def r(numero):
    return np.round(numero,10)

if (dim==2):
    def RK(Bref, GF, Fn, u, uext):
        Ne = np.shape(u)[0]
        nablaGF = np.einsum('ijkl,ik->ijl', GF, np.reshape(u,(Ne,N*4),order='F'))
        fn  = np.zeros((Ne, 4, Nb*3)) 
        front = np.zeros((Ne,N,4))
        F1,F2,F3 = nodes_I2D(PolyDeg)
        for i in range(Ne):  ## Ne
            # nabla[i,:,:] = np.dot(Gx[i,:,:],fx[i,:,:])+np.dot(Gy[i,:,:],fy[i,:,:])
            ## fn = FN with positive eigenvalues*solution in element i + FN with negative eigenvalues*solution in neibourgh element
            ## each line calculates the border integral for each of the triangles sides 
            fn[i,:,0:Nb] = np.matmul(Fn[i,0,0,:,:],np.transpose(uext[i,0,0,:,:]))[:,F1] + np.matmul(Fn[i,0,1,:,:],np.transpose(uext[i,0,1,:,:]))[:,F1]
            fn[i,:,Nb:2*Nb] = np.matmul(Fn[i,1,0,:,:],np.transpose(uext[i,1,0,:,:]))[:,F2] + np.matmul(Fn[i,1,1,:,:],np.transpose(uext[i,1,1,:,:]))[:,F2]
            fn[i,:,2*Nb:3*Nb] = np.matmul(Fn[i,2,0,:,:],np.transpose(uext[i,2,0,:,:]))[:,F3] + np.matmul(Fn[i,2,1,:,:],np.transpose(uext[i,2,1,:,:]))[:,F3]
        for j in range(4):
            front[:,:,j] = r(np.matmul(fn[:,j,:],np.transpose(Bref))) ##r(np.matmul(B,np.transpose(fn[:,j,:]))) # (B*fn[j]')' = fn[j]*B'
        ##np.dot(B[l,0,:],fn[i,0,:,:])+np.dot(B[l,1,:],fn[i,1,:,:])+np.dot(B[l,2,:],fn[i,2,:,:]))
        rhsQ = nablaGF-front
        return rhsQ

elif (dim==3):
    def RK(Bref, GF, Fn, u, uext):
        Ne = np.shape(u)[0]
        # print(np.transpose(np.reshape(u,(Ne,N*(dim+2)),order='F')))
        nablaGF = np.einsum('ijkl,ik->ijl', GF, np.reshape(u,(Ne,N*(dim+2)),order='F'))
        # print(nablaGF)
        fn  = np.zeros((Ne, 5, Nb*4)) 
        front = np.zeros((Ne,N,5))
        F1,F2,F3,F4 = nodes_I3D(PolyDeg)
        for i in range(Ne):  ## Ne
            # nabla[i,:,:] = np.dot(Gx[i,:,:],fx[i,:,:])+np.dot(Gy[i,:,:],fy[i,:,:])
            ## fn = FN with positive eigenvalues*solution in element i + FN with negative eigenvalues*solution in neibourgh element
            ## each line calculates the border integral for each of the triangles sides 
            fn[i,:,0:Nb] = np.matmul(Fn[i,0,0,:,:],np.transpose(uext[i,0,0,:,:]))[:,F1] + np.matmul(Fn[i,0,1,:,:],np.transpose(uext[i,0,1,:,:]))[:,F1]
            fn[i,:,Nb:2*Nb] = np.matmul(Fn[i,1,0,:,:],np.transpose(uext[i,1,0,:,:]))[:,F2] + np.matmul(Fn[i,1,1,:,:],np.transpose(uext[i,1,1,:,:]))[:,F2]
            fn[i,:,2*Nb:3*Nb] = np.matmul(Fn[i,2,0,:,:],np.transpose(uext[i,2,0,:,:]))[:,F3] + np.matmul(Fn[i,2,1,:,:],np.transpose(uext[i,2,1,:,:]))[:,F3]
            fn[i,:,3*Nb:4*Nb] = np.matmul(Fn[i,3,0,:,:],np.transpose(uext[i,3,0,:,:]))[:,F4] + np.matmul(Fn[i,3,1,:,:],np.transpose(uext[i,3,1,:,:]))[:,F4]
        for j in range(5):
            front[:,:,j] = r(np.matmul(fn[:,j,:],np.transpose(Bref))) ##r(np.matmul(B,np.transpose(fn[:,j,:]))) # (B*fn[j]')' = fn[j]*B'
        ##np.dot(B[l,0,:],fn[i,0,:,:])+np.dot(B[l,1,:],fn[i,1,:,:])+np.dot(B[l,2,:],fn[i,2,:,:]))
        rhsQ = nablaGF-front
        return rhsQ

else:
    print('error: dimension not available')

