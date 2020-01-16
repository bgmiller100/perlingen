# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 10:17:39 2019

@author: Ben
"""
# This function simply builds a basic struct that contains the mesh
# information, to be used by the other functions
import numpy as np

def main(Nx=500,Ny=500,pixel_width=.0075):
    # Create a set of points that fall in centres of pixels
    xv = np.linspace(pixel_width/2, pixel_width*(Nx-1/2), Nx);
    yv = np.linspace(pixel_width/2, pixel_width*(Ny-1/2), Ny);
    Y,X = np.meshgrid(yv,xv);
    points = [np.ravel(X,'F'),np.ravel(Y,'F')]; #2,25x, inverse Matlab
    # Store in mesh struct
    mesh = {'points':points, 'Nx':Nx, 'Ny': Ny};
    return mesh 

