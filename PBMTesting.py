import numpy as np
from ValidationMetric import ValidationMetric as VM
from scipy.spatial.distance import pdist,squareform

#create simple numpy array for analysis of PBM validation metric
data = np.ones((4,4))
data[0,:]=data[0,:]*3
data[1,:]=data[1,:]*2
data[2,:]=data[2,:]*8
data[3,:]=data[3,:]*5

#find the center of the data "patterns" (i.e., data points)
dPatCenter = np.mean(data,axis=0)

#find the sum of the distances of all points from the pattern center
dataC = np.vstack([dPatCenter,data])


distances = pdist(dataC)
distances = squareform(distances)
centDists = distances[0,:]
Eo = np.sum(centDists)

#dictionary of the clustered data
dists = {0:{0:[0,1],1:[2,3]}}


VM.PBM(dists,data,4,Eo)
