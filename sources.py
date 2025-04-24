import numpy as np

def monopole(t,s):
    s[:,:,0] = 0
    s[:,:,1] = 0*t
    s[:,:,2] = 0
    s[:,:,3] = 0