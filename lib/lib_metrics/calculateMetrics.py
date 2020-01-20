# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

# This function calculates a set of metrics for a provided pattern that
# quantify the relative positioning of features via the power spectrum.
# Smoothing is applied to the power spectrum, which in almost all cases
# results in ellipsoidal shapes when regions containing the top X% of power
# in the spectrum are taken. The orientation and dimensions of these
# ellipses are the metrics. Usage:
#
# in: pattern:   the pattern for which metrics are to be calculated 
#
# out: metrics:   a (row) vector of the nine metrics relating to the three
#            ellipses

#Credit: http://nicky.vanforeest.com/misc/fitEllipse/fitEllipse.html

import numpy as np
from numpy.linalg import eig, inv, det
#import cv2
from scipy.ndimage.filters import gaussian_filter
from scipy.spatial import ConvexHull
import warnings
#import matplotlib.pyplot as plt

def fitEllipse(x,y):
    x = x[:,np.newaxis]
    y = y[:,np.newaxis]
    D =  np.hstack((x*x, x*y, y*y, x, y, np.ones_like(x)))
    S = np.dot(D.T,D)
    C = np.zeros([6,6])
    C[0,2] = C[2,0] = 2; C[1,1] = -1
    if det(S) == 0:
        raise ValueError
    E, V =  eig(np.dot(inv(S), C))
    n = np.argmax(np.abs(E))
    a = V[:,n]
    return a
def ellipse_angle_of_rotation( a ):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    if b == 0:
        if a > c:
            return 0
        else:
            return np.pi/2
    else:
        if a > c:
            return np.arctan(2*b/(a-c))/2
        else:
            if a-c<10**(-12):
                a=a+10**(-12)
            return np.pi/2 + np.arctan(2*b/(a-c))/2
def ellipse_axis_length( a ):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    up = 2*(a*f*f+c*d*d+g*b*b-2*b*d*f-a*c*g)
    down1=(b*b-a*c)*( (c-a)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    down2=(b*b-a*c)*( (a-c)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    res1=np.sqrt(up/down1)
    res2=np.sqrt(up/down2)
    return np.array([res1, res2])


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
        
        #disp(size(mask))
        #disp(size(S))
        #disp(size(pow_thresh))

        #mask = S > pow_thresh 
        mask[S > pow_thresh] = 1
        # mask = cv2.threshold(np.array(mask*255, dtype=np.uint8), 0,255,cv2.THRESH_BINARY)[1]
        #mask = np.array(mask*255, dtype=np.uint8)
        # Find properties of the best-fit ellipse for this mask
      
        pts = np.transpose(np.asarray(np.where(mask==1))) #access pts[0], pts[1]
        hull = ConvexHull(pts)

        #Test confirms good hull 
        #plt.plot(pts[:,0],pts[:,1],'o')
        #for simplex in hull.simplices:
        #    plt.plot(pts[simplex,0],pts[simplex,1],'k-')
        #plt.show()


        #print('running...')
        np.seterr(all='warn')
        with warnings.catch_warnings():
            warnings.filterwarnings('error')
            try: #catch singular matrices 
                a = fitEllipse(pts[hull.vertices,0], pts[hull.vertices,1])
                angle = ellipse_angle_of_rotation(a)
                if np.isreal(angle): #catch failed angle fitting
                    #print(angle)
                    try: #catch failed axis fitting
                        axes = ellipse_axis_length(a)
                        if all(np.isreal([axes[0],axes[1]])):
                            metrics.extend([np.real(angle)*180/np.pi, np.real(axes[0]), np.real(axes[1])])
                        else:
                            metrics.extend([0,0,0])
                            #print('imaginary axes')
                            #print([angle, axes[0], axes[1]])
                    except Warning:
                        metrics.extend([0,0,0])
                        #print('axis error')
                        pass
                else:
                    #print('imaginary angles')
                    metrics.extend([0,0,0])
            except ValueError:
                #print('fitting error')
                metrics.extend([0,0,0])
                pass
            #print([angle, axes[0], axes[1]])
    return metrics


if __name__ == '__main__':
    pass
