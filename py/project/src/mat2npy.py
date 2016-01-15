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





inputDir = "trained_RBFN/"
outputDir = "trained_RBFN_npz/"
#Load original projections 
for file in os.listdir(inputDir):
    if fnmatch.fnmatch(file, '*.mat'):
        dictionary = sio.loadmat('%s%s'%(inputDir,file))
        Centers=dictionary['Centers']
        Theta=dictionary['Theta']
        betas=dictionary['betas']
        numRBFNeurons=dictionary['numRBFNeurons']
        sigma=dictionary['sigma']
        print(numRBFNeurons.dtype)
        np.savez(outputDir+os.path.splitext(file)[0], Centers=Centers, Theta=Theta,betas=betas,numRBFNeurons=numRBFNeurons,sigma=sigma)
        npzfile = np.load(outputDir+os.path.splitext(file)[0]+'.npz')
        print(npzfile['sigma'])
        
        
        
