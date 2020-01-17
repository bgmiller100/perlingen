# perlingen
Perlin noise generation and state estimation   
ver 1

Benjamin Miller <benjamin.g.miller@utexas.edu>  
The University of Texas at Austin 

01/16/2020

---

## 1. SUMMARY

This algorithm is a Python-based Perlin noise generation and state estimation software developed off of the original MATLAB code by \[1], and is currently implemented on the Stampede2 SKX nodes hosted by the Texas Advanced Computing Center (TACC) at The University of Texas at Austin.  State estimation uses sequential Monte Carlo approximate Bayesian computation following the original software.  As with the original software, processing is currently restricted to 2D datasets.   

> \[1] Jakes, David & Burrage, Kevin & Drovandi, Christopher & Burrage, Pamela & Bueno-Orovio, Alfonso & Santos, Rodrigo & Rodriguez, Blanca & Lawson, Brodie. (2019). Perlin Noise Generation of Physiologically Realistic Patterns of Fibrosis. 10.1101/668848. 

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

v1: Initial private release for testing and development on Stampede2.

---

## 4. OPERATION

On personal Unix systems, "python3 perlingen.py > perlingen.out" may be called directly (saving results and process data).  On the Stampede2 system, "sbatch perlingen.job" should be run, and the internal python call may be edited as needed.  

Control of variables for image generation and analysis is handled by optional command line flags, listed by calling "python3 perlingen.py -h" and set in the standard method, i.e. "python3 perlingen.py --particles 1000". 
Please note, on local PCs, the particle should be manually defined small (~50) to avoid memory toubles. 

After running, the ".out" file will contain process information including memory usage, timeing information, convergence rates, etc., as well as state-space information for the initial and estimated images, which themselves will be stored in the "/results" directory.  

---

## 5. LIBRARY TOC

### root/

#### perlingen.py

#### perlingen.job

### lib_create/

#### generateSeedData 

#### buildMesh

#### generatePatterns

#### createFibroPattern

### lib_metrics/

#### calculateMetrics

#### ellipseDiscrepancy

### lib_smcabc/

#### performSMCABC

### Octave2D

#### Octave2D.c
 
#### setup
	
