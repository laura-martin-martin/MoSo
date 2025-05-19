import math
import numpy as np
from mesh import rstoe, Nodes2D, xytors, xytoab, Nodes3D, xyztorst, rsttoe, xyztoabc, find_bary
from parameters import num
from matrices import Simplex2DP, Simplex3DP
import matplotlib.tri as mtri

dim = num.dim
PolyDeg = num.PolyDeg
N = num.N

if (dim==2):
    def output_mesh(nodes, elements, Ne):    
        
        [x,y] = Nodes2D(PolyDeg)
        [r,s] = xytors(x,y) 
        X = np.zeros((Ne,N))
        Y = np.zeros((Ne,N))
        xplot = np.zeros((Ne,N))
        yplot = np.zeros((Ne,N))
        for i in range(Ne):
                ## positions of the elements' vertices
                [v1x, v1y, _] = nodes["coordinates"][int(elements['triangle3']["nodes"][i][0]-1)]
                [v2x, v2y, _] = nodes["coordinates"][int(elements['triangle3']["nodes"][i][1]-1)]
                [v3x, v3y, _] = nodes["coordinates"][int(elements['triangle3']["nodes"][i][2]-1)]
                x,y = rstoe(r,s,v1x, v1y, v2x, v2y, v3x, v3y)
                xc=(x[0]+x[PolyDeg]+x[-1])/3
                yc=(y[0]+y[PolyDeg]+y[-1])/3
                X[i,:] = np.squeeze(x)
                Y[i,:] = np.squeeze(y)
                alpha = 0.5
                xplot[i,:] = np.squeeze(x)*(1-alpha)+xc*alpha
                yplot[i,:] = np.squeeze(y)*(1-alpha)+yc*alpha
        return(X,Y, xplot, yplot)

    def output_matrix(nodes, elements, X, Y):
        x = np.array(X).flatten()
        y = np.array(Y).flatten()
        triang = mtri.Triangulation(nodes["coordinates"][:,0], nodes["coordinates"][:,1], elements['triangle3']["nodes"]-1)
        trifinder = triang.get_trifinder()
        tri_index = trifinder(x, y)
        v1x = nodes["coordinates"][elements['triangle3']["nodes"]-1][tri_index,0,0]
        v1y = nodes["coordinates"][elements['triangle3']["nodes"]-1][tri_index,0,1]
        v2x = nodes["coordinates"][elements['triangle3']["nodes"]-1][tri_index,1,0]
        v2y = nodes["coordinates"][elements['triangle3']["nodes"]-1][tri_index,1,1]
        v3x = nodes["coordinates"][elements['triangle3']["nodes"]-1][tri_index,2,0]
        v3y = nodes["coordinates"][elements['triangle3']["nodes"]-1][tri_index,2,1]
        a,b = xytoab(x, y, v1x, v1y, v2x, v2y, v3x, v3y)
        V2D = np.zeros((len(a),N))
        sk = 0
        for i in range(PolyDeg+1):
            for j in range(PolyDeg-i+1):
                V2D[:,sk] = Simplex2DP(a,b,i,j)
                sk = sk +1
        print(V2D)
        return(tri_index,V2D)
    
elif(dim==3):
    def output_mesh(nodes, elements, Ne):    
        
        [x,y,z] = Nodes3D(PolyDeg)
        [r,s,t] = xyztorst(x,y,z) 
        X = np.zeros((Ne,N))
        Y = np.zeros((Ne,N))
        Z = np.zeros((Ne,N))
        xplot = np.zeros((Ne,N))
        yplot = np.zeros((Ne,N))
        zplot = np.zeros((Ne,N))
        for i in range(Ne):
                ## positions of the elements' vertices
                v1x, v1y, v1z = nodes["coordinates"][int(elements['tetrahedron4']["nodes"][i][0]-1)]
                v2x, v2y, v2z = nodes["coordinates"][int(elements['tetrahedron4']["nodes"][i][1]-1)]
                v3x, v3y, v3z = nodes["coordinates"][int(elements['tetrahedron4']["nodes"][i][2]-1)]
                v4x, v4y, v4z = nodes["coordinates"][int(elements['tetrahedron4']["nodes"][i][3]-1)]
                x,y,z = rsttoe(r,s,t,v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z, v4x, v4y, v4z)
                xc=(x[0]+x[PolyDeg]+x[int((PolyDeg+1)*(PolyDeg+2)/2-1)]+x[-1])/4
                yc=(x[0]+y[PolyDeg]+y[int((PolyDeg+1)*(PolyDeg+2)/2-1)]+y[-1])/4
                zc=(z[0]+z[PolyDeg]+z[int((PolyDeg+1)*(PolyDeg+2)/2-1)]+z[-1])/4
                X[i,:] = np.squeeze(x)
                Y[i,:] = np.squeeze(y)
                Z[i,:] = np.squeeze(z)
                alpha = 0.5
                xplot[i,:] = np.squeeze(x)*(1-alpha)+xc*alpha
                yplot[i,:] = np.squeeze(y)*(1-alpha)+yc*alpha
                zplot[i,:] = np.squeeze(z)*(1-alpha)+zc*alpha
        return(X,Y,Z, xplot, yplot, zplot)

    def output_matrix(nodes, elements, X, Y, Z):
        x = np.array(X).flatten()
        y = np.array(Y).flatten()
        z = np.array(Z).flatten()
        v1x = nodes["coordinates"][elements['tetrahedron4']["nodes"]-1][:,0,0]
        v1y = nodes["coordinates"][elements['tetrahedron4']["nodes"]-1][:,0,1]
        v1z = nodes["coordinates"][elements['tetrahedron4']["nodes"]-1][:,0,2]
        v2x = nodes["coordinates"][elements['tetrahedron4']["nodes"]-1][:,1,0]
        v2y = nodes["coordinates"][elements['tetrahedron4']["nodes"]-1][:,1,1]
        v2z = nodes["coordinates"][elements['tetrahedron4']["nodes"]-1][:,1,2]
        v3x = nodes["coordinates"][elements['tetrahedron4']["nodes"]-1][:,2,0]
        v3y = nodes["coordinates"][elements['tetrahedron4']["nodes"]-1][:,2,1]
        v3z = nodes["coordinates"][elements['tetrahedron4']["nodes"]-1][:,2,2]
        v4x = nodes["coordinates"][elements['tetrahedron4']["nodes"]-1][:,3,0]
        v4y = nodes["coordinates"][elements['tetrahedron4']["nodes"]-1][:,3,1]
        v4z = nodes["coordinates"][elements['tetrahedron4']["nodes"]-1][:,3,2]
        l0, l1, l2, l3 = find_bary(x, y, z, v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z, v4x, v4y, v4z)
        V3D = np.zeros((len(X),N))
        for i in range(len(X)):
            tri_index = np.where((0<=l0[:,i]) & (l0[:,i]<=1) & (0<=l1[:,i]) & (l1[:,i]<=1) &(0<=l2[:,i]) & (l2[:,i]<=1) & (0<=l3[:,i]) & (l3[:,i]<=1)) ########### ça marche ????
            v1x = nodes["coordinates"][elements['tetrahedron4']["nodes"]-1][tri_index,0,0]
            v1y = nodes["coordinates"][elements['tetrahedron4']["nodes"]-1][tri_index,0,1]
            v1z = nodes["coordinates"][elements['tetrahedron4']["nodes"]-1][tri_index,0,2]
            v2x = nodes["coordinates"][elements['tetrahedron4']["nodes"]-1][tri_index,1,0]
            v2y = nodes["coordinates"][elements['tetrahedron4']["nodes"]-1][tri_index,1,1]
            v2z = nodes["coordinates"][elements['tetrahedron4']["nodes"]-1][tri_index,1,2]
            v3x = nodes["coordinates"][elements['tetrahedron4']["nodes"]-1][tri_index,2,0]
            v3y = nodes["coordinates"][elements['tetrahedron4']["nodes"]-1][tri_index,2,1]
            v3z = nodes["coordinates"][elements['tetrahedron4']["nodes"]-1][tri_index,2,2]
            v4x = nodes["coordinates"][elements['tetrahedron4']["nodes"]-1][tri_index,3,0]
            v4y = nodes["coordinates"][elements['tetrahedron4']["nodes"]-1][tri_index,3,1]
            v4z = nodes["coordinates"][elements['tetrahedron4']["nodes"]-1][tri_index,3,2]
            a = l0[tri_index,i]*(-1) + l1[tri_index,i]*(1) + l2[tri_index,i]*(0) + l3[tri_index,i]**(0)
            b = l0[tri_index,i]*(-0.5773503) + l1[tri_index,i]*(-0.5773503) + l2[tri_index,i]*(1.154701) + l3[tri_index,i]*(0)
            c = l0[tri_index,i]*(-0.4082483) + l1[tri_index,i]*(-0.4082483) + l2[tri_index,i]*(-0.4082483) + l3[tri_index,i]*(1.224745)
            sk = 0
            for i in range(PolyDeg+1):
                for j in range(PolyDeg-i+1):
                    for k in range(PolyDeg-i-j+1):
                        V3D[:,sk] = Simplex3DP(a,b,c,i,j,k)
                        sk = sk +1
            return(tri_index,V3D)
else:
     print('ERROR: Not correct dimension')     