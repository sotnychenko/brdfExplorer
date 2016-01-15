
import copy

from os import listdir
from trainingFunctions import *
from merlFunctions import *
from coordinateFunctions import *
from reconstructionFunctions import *
from numpy import *

OutputDir = "output"
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

##############################################################################################################
#COMPUTE IDs
#In this case we have all the measurements so the IDs vector will span all the possible angles
##############################################################################################################
   

#a=array([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18])
#b = a.reshape(3,2,1,3, order='F')
#print(b.shape)
#print(b)

MERLIds = np.arange(90*90*180)
MERLIds = MERLIds.astype(int)
dataDir = "data/"

##############################################################################################################
#LOAD PRECOMPUTED DATA
##############################################################################################################
maskMap = np.load('%s/MaskMap.npy'%dataDir)	  #Indicating valid regions in MERL BRDFs
median = np.load('%s/Median.npy'%dataDir)	  #Median, learned from trainingdata
cosMap = np.load('%s/CosineMap.npy'%dataDir)  #Precomputed cosine-term for all BRDF locations (ids)
VsU = np.load('%s/VsU.npy'%dataDir)	  #Usv, precomputed pca space for 5 components

relativeOffset = np.load('%s/RelativeOffset.npy'%dataDir) #Offset, learned from trainingdata
Q = np.load('%s/Q.npy'%dataDir) #Scaled eigenvectors, learned from trainingdata

##############################################################################################################
#COMPUTE THE PROJECTION
##############################################################################################################
#Choose which Ids are valid by checking the maskMap

if(MERLIds.dtype != bool):
	sortKeys = np.argsort(MERLIds)	#Sort data so it matches logical indices
	obs = obs[sortKeys]
	MERLIds = MERLIds[sortKeys]
mappedKnownSelector = np.zeros(np.shape(maskMap))	#Create a selector for the mapped version
mappedKnownSelector[MERLIds] = 1
mappedKnownSelector = mappedKnownSelector[maskMap].astype(bool)
obs = obs[maskMap]
MERLIds = MERLIds[maskMap]
print(type(cosMap))
if(cosMap is not None):
	obs *= cosMap[mappedKnownSelector]
#Map the BRDF with their logarithmic mapping
mapped = MapBRDF(obs, maskMap[MERLIds], median[mappedKnownSelector])

#print(mapped.shape)
#Project BRDF to the PCA space
#Note that in the ProjectToPCSpace function they have eta but it has no effect if it's 0 
#In this function they perform SVD again in the basis A (A = PCs[:,0:nPCs]). 
#The purpose here is to get the pseudo-inverse of A
#You were right, if we keep the number of components (nPCs) fixed you can store this value to avoid recomputing it
#proj = VsU.dot(mapped-relativeOffset)  
proj = ProjectToPCSpace(mapped,Q[mappedKnownSelector,:],relativeOffset[mappedKnownSelector], 0) 

#Here I made a copy of the original projections just in case
#new_proj = copy.copy(proj)


##############################################################################################################
#EDITIONS AND COMPUTE THE PROJECTION BACK TO THE ORIGINAL SPACE
##############################################################################################################

#Choose how many components to use
nPCs = 5;

#Perform editions over the projections
proj[0,:] *= 0;

#print(proj[0,:].shape)
#Get the BRDF back from PCA basis to the original space
recon = np.dot(Q[:,0:nPCs],proj)+relativeOffset




#Unmap the bRDF from their logarithmic mapping
if(cosMap is None):
	unm = UnmapBRDF(recon, maskMap, median)
else:
	unm = UnmapBRDF(recon, maskMap, median)
	unm[maskMap] /= cosMap

#Save the BRDF 
saveMERLBRDF('%s/outputName.binary'%OutputDir,unm)
