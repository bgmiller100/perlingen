# -*- coding: utf-8 -*-
"""
@author: Benjamin Miller

Operating interface for launching Perlin noise image generation and Monte-Carlo Bayesian estimation algorithms per [1].  Parses variables, assigns subfunctions, and controls image output. 

IN: Call -h for detailed input options. Includes noise state variables, number of images, number of SMC-ABC particles, tolerable discrepancy, runtime limit, and boolean to run the SMC algorithm after image generation.  

OUT: Input and estimated images as png files.  Text file with corresponding state variables.  Also returns metric sets and discrepancies of each result, not presently incorporated.   
 
[1](Jakes, David & Burrage, Kevin & Drovandi, Christopher & Burrage, Pamela & Bueno-Orovio, Alfonso & Santos, Rodrigo & Rodriguez, Blanca & Lawson, Brodie. (2019). Perlin Noise Generation of Physiologically Realistic Patterns of Fibrosis. 10.1101/668848.)
"""

def check_positive(value):
    ivalue = int(value)
    if ivalue <=0:
    	raise argparse.ArgumentTypeError("%s is an invalid positive int value" % value)
    return ivalue 


def main(num_imgs, feature_size, roughness, fibre_alignment, 
patch_size, fibre_sep, direction,
fibreness, patchiness, density, pixel_width, 
n_fibres_similarity, wiggle_feature_length, phasefield_strength,fibre_period, 
N_parts, desired_D, max_time, run_smcabc):

    """ DICTIONARY (call python3 perlingen.py -h)
    #BASE NOISE (Ob fibrosis field) (impedes lines)
    #   lb = "feature_size" = base obstacle size [.01,2]
    #   gamma = "roughness" = base roughness FIXED[0, 0.99]  
    #   r = "fibre_alignment =  base anisotropy [0.5,50]
    feature_size = .08; #lb
    roughness = 0.1; #gamma
    fibre_alignment = 15; #r
    #DENSITY VARIATION (Od density field)
    #   ld = "patch_size" density feature size [1,8]
    patch_size = 2; #ld
    #FIBRE SELECTION (F sinusoisal field)
    #   L = "fibre_sep" seperation distance [0.3,2]
    #   Phi = "direction"  = fiber orientation [-pi/2, pi/2]
    fibre_sep = 0.6; #L
    direction = 45*np.pi/180; #Phi 
    #MISC
    #   f = "fibreness" = extent of pattern of fibers[0,0.4]
    #   d = "patchiness" magnitude of local density [0, 0.5]
    #   collagen density threshold [0,.99]
    fibreness = .4; #f
    patchiness = .05; #d
    density = 0.2;
	"""
    #IMAGEDATA
    N_patterns = num_imgs;
    Nx = 500;
    Ny = 500;
    #pixel_width = .0075; 
    #SEEDDATA
    N_seeds = 500
    N_freqs = 8; 
    seed = None;
    
    if path.exists('seedinfo.npy')==0:
        generateSeedData.main(N_seeds, N_freqs, seed)
    
    ## RUN   
    #build mesh
    mesh = buildMesh.main(Nx,Ny,pixel_width); #(Nx,Ny,pixel_width);
    params = {'fibreness':fibreness, 'fibre_sep':fibre_sep, 'patchiness':patchiness, 
              'feature_size':feature_size,'roughness':roughness, 'patch_size':patch_size,
              'fibre_alignment':fibre_alignment, 'direction':direction,
	      'n_fibres_similarity':n_fibres_similarity, 'wiggle_feature_length':wiggle_feature_length, 'phasefield_strength':phasefield_strength,'fibre_period':fibre_period}

    #check for results directory
    try:
    	os.makedirs('results')
    except OSError:
    	if not os.path.isdir('results'):
    		raise
    #run generator 
    patterns, img = generatePatterns.main(params, density, N_patterns, mesh, num_imgs);
    target_data = img
    plt.imsave('results/input_target.png',target_data)
    np.save("results/generated_patterns.npy",patterns)
    #run smcabc 
    if run_smcabc==1:
        f_simulate = partial(createFibroPattern.main, mesh, density)
        f_summaries = calculateMetrics.main
        f_discrepancy = ellipseDiscrepancy.main 
        params_mins = np.array([0, 0.3, 0, .01, 0, 1, 0.5, -90]) 
        params_maxs = np.array([0.4, 2, 0.5, 2, .99, 8, 50, 90])
        scale_param = np.array([ 0, 0, 0, 0, 0, 0, 1, 0], dtype=int)
    
        print("running smcabc...\n")
        start_all = time.time()
        part_thetas, part_outputs, part_summaries, part_Ds = performSMCABC.main(N_parts, f_simulate, f_summaries, f_discrepancy, target_data, params_mins, params_maxs, scale_param, desired_D, mesh,max_time)

        print('TOTAL RUNTIME: '+str(time.time()-start_all))
    
        print('printing unique results')
        for i in range(np.shape(part_thetas)[0]):
            outfile = 'results/out_%05d.png' % (i)
            plt.imsave(outfile,np.squeeze(part_outputs[i,:,:]))
        print('...\n')
        print('input: f %.3f,  l %.3f, d %.3f, lb %.3f, g %.3f, ld %.3f, r %.3f, p %.3f'%(params['fibreness'], params['fibre_sep'], params['patchiness'], params['feature_size'], params['roughness'], params['patch_size'], params['fibre_alignment'], params['direction']))

        #save raw data
        np.save("results/part_thetas.npy",part_thetas)
        np.save("results/part_outputs.npy",part_outputs)
        np.save("results/part_summaries",part_summaries)
        np.save("results/part_Ds.npy",part_Ds)

        print('...\n')
        
        for i in range(np.shape(part_thetas)[0]):
            print('out_%05d: f %.3f,  l %.3f, d %.3f, lb %.3f, g %.3f, ld %.3f, r %.3f, p %.3f'%(i, part_thetas[i,0], part_thetas[i,1], part_thetas[i,2], part_thetas[i,3], part_thetas[i,4], part_thetas[i,5], part_thetas[i,6], part_thetas[i,7]))
        print('...\n')
    return 

if __name__=='__main__':
    #from lib.File import Class 
    from os import path
    import numpy as np
    import lib.lib_create.buildMesh as buildMesh
    import lib.lib_create.generateSeedData as generateSeedData
    import lib.lib_create.generatePatterns as generatePatterns
    import lib.lib_create.createFibroPattern as createFibroPattern
    import lib.lib_metrics.calculateMetrics as calculateMetrics
    import lib.lib_metrics.ellipseDiscrepancy as ellipseDiscrepancy
    import lib.lib_smcabc.performSMCABC as performSMCABC
    import os
    import argparse
    from functools import partial 
    import faulthandler
    import time
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    faulthandler.enable()
    print('Launching...')    
    parser = argparse.ArgumentParser(description='Generate Perlin noise')

    parser.add_argument('--lb','--feature_size', type=float,default=.08,help='base obstacle size [.01,2] (.08)',metavar='LB')
    parser.add_argument('--g','--roughness', type=float,default=.1,help='base roughness [0, 0.99] (.1)',metavar='G')
    parser.add_argument('--r','--fibre_alignment', type=float,default=15,help='base anisotropy [.5,50] (15)',metavar='R')

    parser.add_argument('--ld','--patch_size', type=float,default=2,help='density field feature size [1,8] (2)',metavar='LD')

    parser.add_argument('--l','--fibre_sep', type=float,default=.6,help='fibre seperation distance [.3,2] (.6)',metavar='L')
    parser.add_argument('--p','--direction', type=float,default=45,help='fiber orientation [-90,90] (45)',metavar='P')
    parser.add_argument('--fp','--fibre_period', type=float,default=2,help='fibre element widths [x,x] (2)',metavar='FP')

    parser.add_argument('--f','--fibreness', type=float,default=.4,help='mixing extent of fibre pattern [0,0.4] (.4)',metavar='F')
    parser.add_argument('--d','--patchiness', type=float,default=.05,help='mixing local density magnitude [0,.5] (.05)',metavar='D')
    parser.add_argument('--c','--collagen', type=float,default=.2,help='mixing collagen density threshold [0,.99] (.2)',metavar='C')

    parser.add_argument('--num','--number_images',type=check_positive, default=1, help='number of images to generate (1)',metavar='NUM')
    parser.add_argument('--pw','--pixel_width', type=float,default=.0075,help='pixel width, zooming [.001,.05] (.0075)',metavar='PW')
    parser.add_argument('--fs','--n_fibres_similarity', type=float,default=4,help='keep fixed, the similarity off fibres [.5,5] (4)',metavar='FS')
    parser.add_argument('--wl','--wiggle_feature_length', type=float,default=4,help='keep fixed, the length between wave peaks [2,50] (4)',metavar='WL')
    parser.add_argument('--ps','--phasefield_strength', type=float,default=5,help='keep fixed, amplitude of fibre waves [.01,10] (5)',metavar='PS')


    parser.add_argument('--particles','--N_parts', type=check_positive,default=2000,help='number of particles for SMC-ABC [2000 recommended] (50)',metavar='PARTS')
    parser.add_argument('--discrepancy','--desired_D', type=float,default=.5,help='desired discrepancy of error metrics in estimation (.5)',metavar='DISCREP')
    parser.add_argument('--time','--max_time',type=check_positive, default=15*60, help='runtime limit of the smcabc algorithm in seconds, for safety (15*60)',metavar='TIME')
    parser.add_argument('--smcabc','--rumn_smcabc',type=check_positive, default=1, help='0/1 boolean to run the esimation algorithm after target image generation (1)',metavar='SMCABC')

    args = parser.parse_args();
    main(num_imgs = args.num, feature_size = args.lb, roughness = args.g, fibre_alignment = args.r, patch_size = args.ld, fibre_sep = args.l, direction = args.p, fibreness = args.f, patchiness = args.d, density = args.c, pixel_width = args.pw, n_fibres_similarity=args.fs, wiggle_feature_length=args.wl, phasefield_strength=args.ps, fibre_period=args.fp, N_parts = args.particles, desired_D = args.discrepancy, max_time=args.time, run_smcabc=args.smcabc)
    
    
    
