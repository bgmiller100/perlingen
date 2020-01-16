# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 11:07:02 2019

@author: Ben
"""
# This function takes a list of points and a set of parameters (listed 
# below), and creates a pattern of fibrosis accordingly.
#
# Usage:    presence = CreateFibroPattern(points, direction, density, params, Ps, offsets)
#
# INPUTS:   mesh:      information that defines the mesh over which to create a pattern
#           density:   the density of fibrosis
#           params:    a set of parameters defining the noise pattern (see below)
#           Ps:        an m x n matrix, which is m rows of the numbers 0:n-1 arranged in random order (permutation tables for random assignment of vectors in Perlin noise)
#           offsets:   an m x 2 matrix that specifies grid offsets for each octave in octave noise
#
#           ( m is the maximum number of octaves that will be requested by
#           this function, and n is the number of 
#
# OUTPUTS:  presence:  a presence/absence map of fibrosis
#           (O_b):      the Perlin noise field for fibrosis (optional)
#           (O_d):      the Perlin noise field for density variation (optional)
#           (F):       the sinusoidal field used for fibre-type patterns (optional)
#
# PARAMS:   Parameters are provided as a single vector, for ease of
#           interface with other code.
#           Params vector is specified as:
#           [ fibreness, fibre_separation, patchiness, feature_size, roughness, patch_size, alignment_ratio, direction]

import numpy as np
import itertools
from lib.Octave2D import Octave2D

#WE HAVE INTRODUCED ot[3], mesh[Nz], mesh[points][3], params[direction2]
# so 4 new aspects to consider in this thing

def main(mesh, density, params_old, Pt, Ot):
  if not isinstance(params_old, dict):
      params_old = np.transpose(np.squeeze(np.asarray(params_old)))
      params = {'fibreness':params_old[0], 'fibre_sep':params_old[1], 'patchiness': params_old[2], 'feature_size':params_old[3], 'roughness': params_old[4], 'patch_size':params_old[5], 'fibre_alignment': params_old[6], 'direction':params_old[7], 'n_fibres_similarity':4, 'wiggle_feature_length':4, 'phasefield_strength':5, 'fibre_period':2}
  else:
      params=params_old
  # Create a rotated set of points for the application of anisotropy and
  # creation of fibre-aligned pattern. Stored as two rows for ease of matrix
  # transforms and input into C++ functions  
  R_points = np.matmul(np.array([ [np.cos(params['direction']), np.sin(params['direction'])],[-np.sin(params['direction']), np.cos(params['direction'])]]), (mesh['points']));
  
  # Create new permutation tables from the provided by applying it to itself,
  # and then again
  Pt = np.squeeze(Pt)
  offsets = (np.squeeze(Ot))
  Pt2 = np.zeros(np.shape(Pt))
  Pt3 = np.zeros(np.shape(Pt))
  pshape = Pt.shape
  for k in range(pshape[0]): #this just equates, nothing CHANGED!!!
      inds1 = np.array(Pt[k,:]).astype(np.int)
      Pt2[k,:] = Pt[k,inds1]
      inds2 = np.array(Pt2[k,:]).astype(np.int)
      Pt3[k,:] = Pt[k,inds2]
  Pt = list(itertools.chain.from_iterable(np.transpose(Pt)))
  Pt2 = list(itertools.chain.from_iterable(np.transpose(Pt2)))
  Pt3 = list(itertools.chain.from_iterable(np.transpose(Pt3)))
  offsets = list(itertools.chain.from_iterable((offsets)))
  ## SINUSOIDAL NOISE FOR FIBROUS PATTERNS
  # Define parameters for the phase field (decided via experimentation)
  n_fibres_similarity = params['n_fibres_similarity']
  wiggle_feature_length = params['wiggle_feature_length']
  phasefield_strength = params['phasefield_strength']
      
  # Create an octave noise field with these properties (first transform
  # points, then call Octave2D)
  phasefield_points = [ [R_points[0,:] / wiggle_feature_length] ,[R_points[1,:] / (n_fibres_similarity * params['fibre_sep'])] ]; # Scale dimensions according to anisotropy (multiplication by "D" in paper)
  phasefield_points = np.squeeze(phasefield_points)
  phasefield_points = list(itertools.chain.from_iterable(np.transpose(phasefield_points)))#MODIFIED!!!
  #print(Pt[100]); reasonable
  phasefield = Octave2D.Octave2D(phasefield_points, 4, 0.5, Pt, offsets);    # (4 octaves, low roughness for phasefield)
  #returns a 1D list, convert to array
  ##print("\n\nTHIS IS PRE-ASARRAY, AFTER OCTAVE\n")
  ##print(phasefield[0:20])
  #print("\n\n")
  
  phasefield = np.asarray(phasefield);
  #print(phasefield)#WAY TOO SMALL - -310 power, not 0.4-0.6
  # Create the sinusoidal pattern using cos(RX^T [0; 1] ), then modulate phase using
  # the created phasefield (fibre_period = 2)
	
  F = 0.5 + 0.5 * np.cos(params['fibre_period']*np.pi * (R_points[1,:] / params['fibre_sep'] + phasefield_strength * (phasefield - 0.5) ) );

  # Sharpen this field by powering it up many times (hard-coded for
  # efficiency, but can be easily modified if a different 'sharpening factor'
  # is desired). 
  F = np.multiply(np.multiply(np.multiply(np.multiply(np.multiply(np.multiply(
          np.multiply(np.multiply(np.multiply(np.multiply(np.multiply(np.multiply(
          np.multiply(np.multiply(np.multiply(F,F),F),F),F),F),F),F),F),F),F),F),F),F),F),F)
  
  
  ## CREATE THE MAIN FIBROSIS DEPOSIT EFFECT
  # Transform points according to input parameters, then call Octave2D
  P_f_points = np.array([ [R_points[0,:] / np.sqrt(params['fibre_alignment'])] ,[R_points[1,:] * np.sqrt(params['fibre_alignment'])] ]);
  P_f_points = np.transpose(np.squeeze(P_f_points))
  #print(np.shape(P_f_points))
  P_f_points = list(itertools.chain.from_iterable(P_f_points/params['feature_size']))
  O_b = Octave2D.Octave2D(P_f_points, 4, params['roughness'], Pt2, offsets);
  O_b = np.asarray(O_b);

  ## CREATE A LARGE-SCALE PERLIN NOISE PATTERN FOR DENSITY VARIATION
  # Use Octave2D with scaling of point co-ords to attain desired patch_size
  mesharray = np.transpose(np.squeeze(np.array(mesh['points'])))
  #print(np.shape(mesharray))
  Od_points = list(itertools.chain.from_iterable(mesharray/params['patch_size']))
  O_d = Octave2D.Octave2D(Od_points, 4, 0.5, Pt3, offsets);
  O_d = np.asarray(O_d);

  ## TAKE A COMBINATION OF THESE NOISEFIELDS TO GET THE FINAL PATTERN
  noise = np.multiply(( 1 - params['fibreness'] + params['fibreness'] * F ),O_b) + params['patchiness'] * O_d;

  ## THRESHOLD THIS NOISE TO GET PRESENCE/ABSENCE OF REQUESTED DENSITY
  #presence = thresholdPattern(noise, params['density']);
  presence = noise
  
  ## CONVERT BACK TO MATRICES
  presence = np.reshape(presence,(mesh['Nx'], mesh['Ny']))
  F = np.reshape(F, (mesh['Nx'], mesh['Ny']))
  O_b = np.reshape(O_b, (mesh['Nx'], mesh['Ny']))
  O_d = np.reshape(O_d, (mesh['Nx'], mesh['Ny']))

  ## FLIP TO CONVERT BACK TO IMAGES (starting at top left instead of bottom left)
  presence = np.flipud(presence);
  F = np.flipud(F);
  O_b = np.flipud(O_b);
  O_d = np.flipud(O_d);
  
  #return presence, O_b, O_d, F
  return presence  

if __name__ == '__main__':
    pass
