# -*- coding: utf-8 -*-
"""
Created on Sun Nov  3 10:43:12 2019

@author: Ben
"""
# This function takes the provided set of parameter values, and creates the
# requested number of representative patterns.
#
# INPUTS:
#
# params - the parameter values for the generator to use
# density - density of fibrosis in the patterns to be generated
# N_patterns - the number of patterns to generate
# (mesh) - optionally provided mesh to specify size of patterns


#import matplotlib.pyplot.imshow as imshow
import numpy as np
#import matplotlib
#matplotlib.use('TkAgg')
#import matplotlib.pyplot as plt
import random
import lib.lib_create.createFibroPattern as createFibroPattern
from cv2 import imwrite

def main(params, density, N_patterns, mesh, num_imgs):
    
    #Load in the seed data
    seedinfo = np.load("seedinfo.npy",allow_pickle=True);
    Pt = seedinfo.item().get('permute_tables');
    Ot = seedinfo.item().get('offset_tables');
    N_seeds = np.shape(Pt)
    if num_imgs>N_seeds[0]:
        print('Too many images requested.\n')
        print('Max images under current perumtation seed is: %d' % (N_seeds[0]))
        return []
    # Define a 'fibrosis' colormap
    fibroclr = [[0.95, 0.85, 0.55],[0.8, 0.2, 0.2]]
    
    #Create mesh if unimplemented 
    
    #Create the requested number of patterns 
    #to be implemented later
    patterns = []
#getrandom m
    mlist = random.sample(range(N_seeds[0]),num_imgs)
    for ind, m in enumerate(mlist):
        #excludes 'fibre-free' case
        #unrandomized permutations 
        #presence,  O_b, O_d, F = createFibroPattern.main(mesh, density, params, Pt[m,:,:], Ot[m,:,:]);
        presence = createFibroPattern.main(mesh, density, params, Pt[m,:,:], Ot[m,:,:]);
        # storage
        patterns.append(presence)
        outfile = '%s/%s.png' % ('results','out_'+str(ind))
        #plt.imsave(outfile, presence)
        imwrite(outfile,presence)   
    #Initialize a figure to plot some patterns
    #Plot the first 10 or less patterns 
    #####imshow(patterns[0,:,:], cmap=fibroclr)
    last_pattern = presence
    return patterns, last_pattern
