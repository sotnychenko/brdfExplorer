import numpy as np
import scipy as sp
import scipy.io as sio
import os.path as path
import copy
import math
import fnmatch
import os
from fractions import Fraction    
import cmath

from os import listdir
from trainingFunctions import *
from merlFunctions import *
from coordinateFunctions import *
from reconstructionFunctions import *

def rgb2lab(r):
    xyz = np.dot(rgb2xyz,r)
    X = xyz[0,:] / 0.950456
    Y = xyz[1,:]
    Z = xyz[2,:] / 1.088754
    T = 0.008856
     
    idx = [i for i,v in enumerate(X) if v > T]
    fX = (7.787*X + 16.0/116.0)    
    fX[idx] = np.sign(X[idx])*pow(abs(X[idx]),Fraction('1/3'))
    
    Y3 = np.sign(Y)*pow(abs(Y),Fraction('1/3'))    
    idx = [i for i,v in enumerate(Y) if v > T]
      
    
    fY = (7.787*Y + 16.0/116.0)
    L  = (903.3*Y)    
    
    fY[idx] = Y3[idx]
    L[idx]    = (116.0*Y3[idx] - 16.0)
    
    idx = [i for i,v in enumerate(Z) if v > T]
    
    fZ =(7.787*Z + 16.0/116.0)
    fZ[idx] = np.sign(Z[idx])*pow(abs(Z[idx]),Fraction('1/3'))
    a = 500.0*(fX - fY)
    b = 200.0*(fY - fZ)
    # print X
    # print Y
    # print Z
    # print L
    # print a
    # print b
    
    return(L,a,b)

def lab2rgb(L,a,b):    
    T1 = 0.008856
    T2 = 0.206893
    
    fY = pow(((L + 16.0)/116.0),3)
    idx = [i for i,v in enumerate(fY) if v > T1]
    fY = (L/903.3)
    
    fY[idx] = pow(((L[idx] + 16.0)/116.0),3)
    Y = fY

    fY = (7.787*Y + 16.0/116.0)
    fY[idx] = np.sign(Y[idx])*pow(abs(Y[idx]),Fraction('1/3'))
    
    fX = a/500.0 + fY
    idx = [i for i,v in enumerate(fX) if v > T2]
    
    X = (fX - 16.0/116.0)/7.787
    X[idx] = pow(fX[idx],3.0)
    
    fZ = fY - b/200.0;
    
    idx = [i for i,v in enumerate(fZ) if v > T2]
    
    Z = (fZ - 16.0/116.0)/7.787
    Z[idx] = pow(fZ[idx],3.0)
    
    
    X = X * 0.950456
    Z = Z * 1.088754

    rgb=np.dot(xyz2rgb,[X, Y, Z])
    return(rgb)
    
    
dataDir = "data"
##############################################################################################################
#LOAD PRECOMPUTED DATA
##############################################################################################################
maskMap = np.load('%s/MaskMap.npy'%dataDir)      #Indicating valid regions in MERL BRDFs
median = np.load('%s/Median.npy'%dataDir)      #Median, learned from trainingdata
cosMap = np.load('%s/CosineMap.npy'%dataDir)  #Precomputed cosine-term for all BRDF locations (ids)
relativeOffset = np.load('%s/RelativeOffset.npy'%dataDir) #Offset, learned from trainingdata
Q = np.load('%s/Q.npy'%dataDir) #Scaled eigenvectors, learned from trainingdata


rgb2xyz = np.array([[0.412453, 0.357580, 0.180423],
                    [0.212671, 0.715160, 0.072169],
                    [0.019334, 0.119193, 0.950227]])
xyz2rgb = np.array([[3.240479, -1.537150, -0.498535],
                    [-0.969256,     1.875992,    0.041556],
                    [0.055648, -0.204043,  1.057311]])
                    
xyz2lab     = np.array([[ 0.00,    116.00,      0.00,      -16.00],
                      [500.00, -500.00,      0.00,        0.00],
                      [0.00,    200.00,     -200.00,    0.00]])
                      

inputDir = "PCA-projections/"
#Load original projections 
for file in os.listdir(inputDir):
    if fnmatch.fnmatch(file, '*.mat'):
        dictionary = sio.loadmat('%s%s'%(inputDir,file))
        projection = dictionary['proj'] 
        [L,a,b] = rgb2lab(projection.T)
    
#print(projection)
#print(L)
#print(a)
#print(b)

inputDir = "Lab-edits/"

nPCs = 5
for file in os.listdir(inputDir):
    if fnmatch.fnmatch(file, '*.mat'):
        dictionary = sio.loadmat('%s%s'%(inputDir,file),squeeze_me=True)
        L = np.array(dictionary['xnew']) 
        #Use this new edited L together with the original brdf a and b
        rgb = lab2rgb(L,a,b)
        #print(L)
        recon_edit = np.dot(Q[:,0:nPCs],rgb.T)+relativeOffset
        if(cosMap is None):
            unmapped_edit = UnmapBRDF(recon_edit, maskMap, median)
        else:
            unmapped_edit = UnmapBRDF(recon_edit, maskMap, median)
            unmapped_edit[maskMap] /= cosMap
        
        saveMERLBRDF('%s.binary'%(file[0:-4]),unmapped_edit)

        
        
        
