# perlingen
Perlin noise generation and state estimation   
ver 1.2

Benjamin Miller <benjamin.g.miller@utexas.edu>  
The University of Texas at Austin 

01/23/2020

---

## 1. SUMMARY

This algorithm is a Python-based Perlin noise generation and state estimation software developed off of the original MATLAB code by \[1], and is currently implemented on the Stampede2 SKX nodes hosted by the Texas Advanced Computing Center (TACC) at The University of Texas at Austin.  State estimation uses sequential Monte Carlo approximate Bayesian computation following the original software \[2].  As with the original software, processing is currently restricted to 2D datasets, but the error metric calculations have been updated via \[3] for future extension to 3D datasets.  Higher efficacy found when using Effective Sample Size (ESS) adaptive tempering  \[4].    

> \[1] Jakes, David & Burrage, Kevin & Drovandi, Christopher & Burrage, Pamela & Bueno-Orovio, Alfonso & Santos, Rodrigo & Rodriguez, Blanca & Lawson, Brodie. (2019). Perlin Noise Generation of Physiologically Realistic Patterns of Fibrosis. 10.1101/668848. 

> \[2] Drovandi, C. C. and Pettitt, A. N. (2011). Estimation of Parameters for Macroparasite Population Evolution using Approximate Bayesian Computation. Biometrics, 67(1):225-233.

> \[3] Judd, Tom. “Fit Data Points to an Ellipse.” JuddZone, 30 Oct. 2017, www.juddzone.com/ALGORITHMS/least_squares_ellipse.html

> \[4] Gunawan, David & Dang, Khue-Dung & Quiroz, Matias & Kohn, Robert & Rran, Minh-Ngoc. (2018). Subsampling Sequential Monte Carlo for Static Bayesian Models. arXiv: Computation.  

---

## 2. REQUIREMENTS

Python 3.7.0

lhsmdu(0.1)
scipy (1.2.0) /numpy (1.16.1) 
matplotlib (2.2.3)

Tested on Ubuntu 18.04 LTS
Functional on TACC Stampede2 SKX 

---

## 3. VERSION NOTES

v1.0: 01/16/2020: Initial private release for testing and development on Stampede2.

v1.1: 01/23/2020: Updated documentation, argument inputs, and output strategy.  
Confirmed success on Stampede2 using a 2000-particle set on all cores, taking 13:12:40 to reach degeneracy on 96 cores, achieving  .2-.4 discrepancy from an intial 6.7 discrepancy target.

v1.2: 07/07/2020: Adaptive tempering update to performSMCABC.  Initial MPI restructuring in perlingenMPI (under development).  

---

## 4. OPERATION

On personal Unix systems, "python3 perlingen.py > perlingen.out" may be called directly (saving results and process data).  On the Stampede2 system, "sbatch perlingen.job" should be run, and the internal python call may be edited as needed.  

Control of variables for image generation and analysis is handled by optional command line flags, listed by calling "python3 perlingen.py -h" and set in the standard method, i.e. "python3 perlingen.py --particles 1000". 
Please note, on local PCs, the particle should be manually defined small (~50) to avoid memory toubles.  The three most important values to set, both locally and on Stampede2, are "--particles, --discrepancy, --time".   

After running, the ".out" file will contain process information including memory usage, timeing information, convergence rates, etc., as well as state-space information for the initial and unique estimated images.  The unique images will be saved in "/results/", along with several ".npy" data files storing the SMC output. 

---

## 5. LIBRARY TOC
Documentation below taken in part from the original MATLAB code by Brodie Lawson.  Credit provided in code.

### root/

#### perlingen.py
Operating interface for launching Perlin noise image generation and Monte-Carlo Bayesian estimation algorithms per [1].  Parses variables, assigns subfunctions, and controls image output. 

#### perlingen.job
Sbatch file required to run "perlingen.py" on Stampede2.  All input options should be specified in this file.  

#### perlingenMPI.py
Combines perlingen.py and performSMCABC for MPI parallelization, rather than using the multiprocessing module, in order to scale up capability at high particle counts (wide state seach space).  

### lib_create/

#### generateSeedData 
Creates seed data and saves it as seedinfo.npy.  Only needs to be called when seed data needs to be changed
or (re)generated.

#### buildMesh
Builds a basic struct that contains the mesh information, to be used by the other functions.

#### generatePatterns
Takes the provided set of parameter values, and creates the requested number of representative patterns.  Can be run alone without the state estimation algorithm by setting option "--smcabc 0".

#### createFibroPattern
Takes a list of points and a set of parameters and creates a pattern of fibrosis accordingly, using a triple-convolution system of a sinusoidal, density variation, and base noise field.  Field generation is handled in "Octave2D.c".

### lib_metrics/

#### calculateMetrics
Calculates a set of metrics for a provided pattern that quantify the relative positioning of features via the power spectrum. Smoothing is applied to the power spectrum, which in almost all cases results in ellipsoidal shapes when regions containing the top X% of power in the spectrum are taken. The orientation and dimensions of these
ellipses in 10% intervals are the metrics, resulting in a 27D vector.  For N-dimensional analysis in Python, see \[3]. 

#### ellipseDiscrepancy
Calculates a measure of discrepancy between a set of target and estimation pattern metrics. The latter is distinguished because the eccentricity of the ellipses in the target pattern is used to further weight the Mahalnobis distances calculated for angles.

### lib_smcabc/
  
#### performSMCABC
Performs SMC-ABC as laid out by \[2].  Sets of particles satisfying a set of intermediary "discrepancy cutoffs"
are generated by successive resampling and mutation steps in the vein of traditional SMC. The function is designed to be modular, so accepts as input the functions that will be used to simulate the model, and to calculate the discrepancy with the target data.
MCMC step is equipped with 0.8 ESS adaptive tempering to 100% over 20 iterations \[4]. 

### Octave2D

#### Octave2D.c
C extension for Python3. Creates an 'octave noise' pattern (multiple octaves of basic Perlin noise, evaluated at a set of points provided by the user. A grid with unit spacing is used for the Perlin grid - input points should be transformed appropriately to generate features of desired size and  orientation. The user may also specify the number of octaves, and the factor successively applied to each octave. This seeded version of the code also takes as input a permutation table of the numbers 0 to 255 and an m x 2 set of offsets in 2D space to apply to the Perlin grid for each octave. m must therefore be at least as large as N_freq.
 
#### setup
Compiles the C extension into an ".so/.o" file pair.  In current verison, these must be pulled out of subdirectories manually upon compiling.  
	
