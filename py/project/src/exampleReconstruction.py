import numpy as np
from merlFunctions import *
from coordinateFunctions import *
from reconstructionFunctions import *
from os import listdir
#Path of the BRDF you want to project
MERLDir = "test/"

##############################################################################################################
#READ THE BRDF
##############################################################################################################
materials = [brdfFile for i,brdfFile in enumerate(listdir(MERLDir)) \
                if (path.isfile(MERLDir+brdfFile)
                and path.splitext(brdfFile)[1] == ".binary")]          

#Fill observation array
obs = np.zeros((90*90*180, 3*len(materials)),'float32')
#Add each color channel as a single observation
for i in range(0,len(materials)):
    mat = readMERLBRDF("%s/%s"%(MERLDir,materials[i]))
    obs[:,3*i] = np.reshape(mat[:,:,:,0],(-1))
    obs[:,3*i+1] = np.reshape(mat[:,:,:,1],(-1))
    obs[:,3*i+2] = np.reshape(mat[:,:,:,2],(-1))
#Convert to BRDF coordinates

#Convert to IDs (i.e. rows-ids in the PC matrix)
MERLIds = np.array(range(90*90*180))

print type(MERLIds)

#Load precomputed data
dataDir = "data/"
maskMap = np.load('%s/MaskMap.npy'%dataDir)   #Indicating valid regions in MERL BRDFs
median = np.load('%s/Median.npy'%dataDir)     #Median, learned from trainingdata
#cosMap = np.load('%s/CosineMap.npy'%dataDir)  #Precomputed cosine-term for all BRDF locations (ids)
cosMap = None
relativeOffset = np.load('%s/RelativeOffset.npy'%dataDir) #Offset, learned from trainingdata
Q = np.load('%s/ScaledEigenvectors.npy'%dataDir) #Scaled eigenvectors, learned from trainingdata

#Reconstruct BRDF
recon = ReconstructBRDF(obs, MERLIds, maskMap, Q, median, relativeOffset, cosMap, eta=40)
#Save reconstruction as MERL .binary
saveMERLBRDF("reconstruction.binary",recon)