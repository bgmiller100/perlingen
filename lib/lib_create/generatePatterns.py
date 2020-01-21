# -*- coding: utf-8 -*-
"""
@author: Brodie Lawson et al. (see README), Benjamin Miller 

This function takes the provided set of parameter values, and creates the
requested number of representative patterns.

INPUTS:
params - the parameter values for the generator to use
density - density of fibrosis in the patterns to be generated
N_patterns - the number of patterns to generate
mesh - mesh to specify size of patterns

OUTPUTS:
results/input_xxxx.png - image data
"""


import numpy as np
import random
import lib.lib_create.createFibroPattern as createFibroPattern
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

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

    #Generate and store images
    patterns = []
    mlist = random.sample(range(N_seeds[0]),num_imgs)
    for ind, m in enumerate(mlist):
        presence = createFibroPattern.main(mesh, density, params, Pt[m,:,:], Ot[m,:,:]);
        patterns.append(presence)
        outfile = 'results/input_%05d.png' % (ind)
        plt.imsave(outfile, np.squeeze(presence))

    last_pattern = presence
    return patterns, last_pattern
