import math
import numpy as np
from mesh import rstoe, Nodes2D, xytors, xytoab
from parameters import num
from matrices import Simplex2DP
import matplotlib.tri as mtri
import matplotlib.pyplot as plt

dim = num.dim
PolyDeg = num.PolyDeg
N = num.N

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