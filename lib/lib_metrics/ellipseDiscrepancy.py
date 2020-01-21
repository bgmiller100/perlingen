# -*- coding: utf-8 -*-
"""
@author: Brodie Lawson et al. (see README), Benjamin Miller 

This function calculates a measure of discrepancy between a set of
metrics associated with a pattern, and the metrics associated with a
target pattern. The latter is distinguished because the eccentricity of
the ellipses in the target pattern is used to further weight the
distances calculated for angles. If the optional third argument is
provided, the distance is the Mahalanobis distance, otherwise it is
Euclidean (except with the period nature of angles taken into account).

INPUTS:
metrics:         a row vector (or matrix) of metrics for a single 
                 (multiple) pattern(s)
target_metrics:  the metrics for the target pattern
(invC):          inverse covariance matrix for Mahalanobis distance

OUTPUT:
D:               discrepancies for each pattern
"""
import numpy as np

def angle_dists(angles1, angle2):
    # This function calculates the distance between two angles, taking into
    # account the periodic nature of alignment angles (i.e. 89 degrees is two
    # degrees away from -89 degrees)
    
    # Store indices of angles where periodicity needs to be taken into account
    p = np.abs( angles1 - angle2 ) > 90
    
    # Also store indices of cases where angle1 is bigger than angle2
    b = angles1 > angle2
    
    # Now calculate distances using these indices
    dists = np.zeros(np.shape(angles1))
    dists[p&b] = angle2 + 180 - angles1[p&b];
    dists[p&(~b)] = angles1[p&(~b)] + 180 - angle2;
    dists[~p] = angle2 - angles1[~p];
    return dists

def main(metrics_old, target_metrics_old, invC):
    # First, calculate the base distances between the metrics
    metrics = np.array(metrics_old)
    target_metrics = np.array(target_metrics_old)
    if len(np.shape(metrics))==1:
        metrics = np.reshape(metrics, (1,target_metrics.size))
        dM = metrics - target_metrics
    else:
        dM = metrics
        for i in range(np.shape(metrics)[0]):
            dM[i,:] = dM[i,:] - target_metrics
   
    # The metrics vectors are of length (3 x <number of ellipses>)
    # Because angle metrics have special handling, a loop is used to modify all
    # of these. Structure is   < orientation, major_axis_length, minor_axis_length >
    for k in range(0,target_metrics.size,3):
        
        # Angles are periodic, so take this into account when calculating the
        # distances between angles
        dM[:,k] = angle_dists( metrics[:,k], target_metrics[k] )
        # Use the eccentricities of ellipses in the target pattern as scaling 
        # factors for the discrepancies in angle metrics
        dM[:,k]  = dM[:,k] * np.absolute( np.log( target_metrics[k+1] / target_metrics[k+2] ) )

    # If an (inverse) covariance matrix has been supplied, calculate the
    # distance using Mahalanobis distance, otherwise use Euclidean
    #if invC != 0:
        # Uses a sum trick so as to only calculate diagonal elements of what would otherwise be two full matrix products
    D = np.sqrt( np.sum(np.multiply(np.matmul(dM , invC) , dM), 1) )
    return D


if __name__ == '__main__':
    pass
