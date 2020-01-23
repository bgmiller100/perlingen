# -*- coding: utf-8 -*-
"""
@author: Brodie Lawson et al. (see README), Benjamin Miller 

This function calculates a set of metrics for a provided pattern that
quantify the relative positioning of features via the power spectrum.
Smoothing is applied to the power spectrum, which in almost all cases
results in ellipsoidal shapes when regions containing the top X% of power
in the spectrum are taken. The orientation and dimensions of these
ellipses are the metrics. 

For N-dimensional analysis in Python, see:
Judd, Tom. “Fit Data Points to an Ellipse.” JuddZone, 30 Oct. 2017, www.juddzone.com/ALGORITHMS/least_squares_ellipse.html. 

in: pattern:   the pattern for which metrics are to be calculated 

out: metrics:   a (row) vector of the three metrics relating to the nine ellipses

"""
import numpy as np
from numpy.linalg import eig, inv, pinv, det
from scipy.ndimage.filters import gaussian_filter
from scipy.spatial import ConvexHull
from scipy.ndimage.morphology import binary_dilation
import warnings
#import matplotlib.pyplot as plt


def fitEllipse2(x,y):
    #This function gets fits an ellipse to the hull of a 2D Fourier power spectrum 
    x = x[:,np.newaxis]
    y = y[:,np.newaxis]
    
    #Perform least-squares fitting of data (J) to ellipse polynomial coefficients (Amat)
    J =  np.hstack((x*x, x*y, y*y, x, y)) 
    S = np.dot(J.T,J)
    ABC = np.dot(pinv(S),np.dot(J.T,np.ones_like(x)))
    a = np.append(ABC,-1)
    Amat = np.array([[a[0],a[1]/2.,a[3]/2.], [a[1]/2., a[2], a[4]/2.], [a[3]/2.,a[4]/2.,a[5]]])
    
    #Convert polynomial coefficients to tranformation mapping
    A2 = Amat[0:2,0:2]
    cc = -np.dot(inv(A2),a[3:5]/2.)
    Tofs = np.eye(3)
    Tofs[2,0:2]=cc
    R = np.dot(Tofs,np.dot(Amat,Tofs.T))
    RS = R[0:2,0:2]/(-R[2,2])
    (el,ec) = eig(RS)
    try: 
        recip = 1./np.abs(el)
    except:
        recip = 0.1
        print('recip err')
    #Get axes and orientation data 
    axes = np.sqrt(1./np.abs(el))
    deg = np.degrees(np.arctan2(ec[1,0],ec[0,0]))
    return np.array([deg, axes[0], axes[1]])

def main(pattern):
   
    #list of threshold values for %power in total spectrum
    power_thresholds = [.1,.2,.3,.4,.5,.6,.7,.8,.9]
    ## Calculate the nine FFT ellipse-derived metrics
        
    # Read out number of thresholds supplied
    N_thresholds = len(power_thresholds)
    
    # Calculate the power spectrum of the (mean-subtracted) pattern
    P = np.absolute(np.fft.fftshift(np.fft.fft2(pattern - np.mean(pattern))))
    P = P**2
    
   # Apply smoothing in order to reduce noise
    #S = cv2.GaussianBlur(P,(5,5),4)
    S = gaussian_filter(P, 4)
    # Calculate the total power contained
    P_tot = np.sum(S)
    
    # Using a cumulative sum of power values sorted in descending order,
    # thresholds for values that represent X% of the power can be found
    Ss =np.flip(np.sort(np.ravel(S)))
    Sss = np.cumsum(Ss)
    
    # Loop over thresholds supplied, creating a mask for each, and then using
    # this to calculate ellipse properties that become the metrics
    metrics = []
    for k in range(0,N_thresholds):
        
        # Find the threshold power value that defines cutoff for top X% of
        #power
        allidx = np.squeeze(np.nonzero(Sss>(power_thresholds[k]*P_tot)))
        pow_thresh = Ss[allidx[0]]
        # Create a binary mask denoting which locations do fall within this
        # threshold
        mask = np.zeros(np.shape(S));
        mask[S > pow_thresh] = 1

        # Find properties of the best-fit ellipse for this mask
        pts = np.transpose(np.asarray(np.where(mask==1))) #access pts[0], pts[1]
        hull = ConvexHull(pts)
        #Dilation afety check ensures proper data set size for ellipse interpolation
        if len(pts[hull.vertices,0]) <6: 
            mask = binary_dilation(mask)
            pts = np.transpose(np.asarray(np.where(mask==1))) 
            hull = ConvexHull(pts)
      
        #Test confirms good hull 
        #plt.plot(pts[:,0],pts[:,1],'o')
        #for simplex in hull.simplices:
        #    plt.plot(pts[simplex,0],pts[simplex,1],'k-')
        #plt.show()
        
        #Run ellipse fitting with error handling 
        np.seterr(all='warn')
        with warnings.catch_warnings():
            warnings.filterwarnings('error')
            try: #catch singular matrices 
                metrics.extend(fitEllipse2(pts[hull.vertices,0], pts[hull.vertices,1]))
            except ValueError:
                print('fitting error')
                metrics.extend([0,0,0])
                pass
        
    return metrics


if __name__ == '__main__':
    pass
