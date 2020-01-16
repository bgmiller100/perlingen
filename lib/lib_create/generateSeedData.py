# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 10:55:54 2019

@author: Ben
"""
# This function creates seed data and saves it as fibro_seedinfo.mat
# This function only needs to be called when seed data needs to be changed
# or (re)generated. Usage:
#
# generateSeedData( N_seeds, (N_offsets), (seed) );
#
# N_seeds:   number of unique seeds to generate
# N_freqs:   optional argument specifying the maximum number of octave
#            layers that will ever be required (default: 8)
# seed:      an integer that, if provided, will be used for seeding
#            MATLAB's random number generator

import random 
import numpy as np

def main(N_seeds=500, N_freqs=8, seed=None):
    random.seed(seed)
    #Create the requested number of permutation and offset tables
    permute_tables = np.empty([N_seeds,N_freqs,256]);
    offset_tables = np.empty([N_seeds,N_freqs,2]);
    for k in range(N_seeds):
    # Permutation tables for this seed
        for j in range(N_freqs):
            data = np.random.permutation(256);
            permute_tables[k,j,:] = data.astype(np.int32); 
    # Offset table for this seed
        offset_tables[k,:,:] = np.random.rand(N_freqs, 2) - 0.5;
    
    #Store Data
    seedinfo = {'permute_tables': permute_tables, 'offset_tables':offset_tables};
    np.save("seedinfo.npy",seedinfo)

