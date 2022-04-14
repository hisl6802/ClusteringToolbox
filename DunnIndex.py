import numpy as np
import time
import statistics as stat
from scipy.spatial.distance import pdist,squareform

def dunnIndex(data,dists,num_groups):
    '''
    Dunn Index (1974?)
    Determine the appropriate number of clusters by threading through the data set one clustering at a time.
    
    Input:
    data - dictionary of the clusters for validation, or raw data for the interpretation of the appropriate number of clusters...
    dists - standardized or submitted non-standardized data. 
    num_groups - number of groups in the data set. 

    Output:

    Array of the number of clusters and the validation index measure. 
 
    '''

    #grab the input dictionary size
    clusterings = len(data)

    startPoint = 0.5*clusterings
    startPoint = int(startPoint)
    numIts = clusterings - startPoint

    numClusters = clusterings-startPoint
    numClusters = int(numClusters)
    val_index = np.zeros((2,clusterings))
    initTime = time.time()

    for i in range(startPoint,clusterings):
        
        #used to analyze the performance of the algorithm
        startTime = time.perf_counter()


        #grab the current set of metabolite clusters
        curClusters = data[i]

        #from the current clusters determine the length in order to determine the next step
        curClustersLength = len(curClusters)

        #sum of intra cluster distances
        sumIntra = 0

        #creating a dictionary of the dispersion of the cluster.
        dispersion = curClusters

        #create a numpy array for cluster centers
        centersNum = np.zeros((curClustersLength,num_groups))

        #create a numpy array for the maximum spread of points in cluster
        maxIntra = np.zeros((curClustersLength,1))

        for j in range(curClustersLength):
            #pull out the current cluster of metabolites
            cluster = curClusters[j]

            #determine whether the current cluster is a list or integer
            if isinstance(cluster,list):
                #check the length of the cluster, get coordinates
                lengthList = len(cluster)
                clustCoordinates = np.zeros((lengthList,num_groups))

                #put the cluster coordinates into the numpy array
                for k in range(lengthList):
                    clustCoordinates[k,:] = dists[cluster[k]]

                #creating a numpy array for the cluster center
                center = np.zeros((1,num_groups))

                for m in range(num_groups):
                    center[0,m] = stat.mean(clustCoordinates[:,m])
                
                #update the numpy array of cluster centers
                centersNum[j,:] = center

                print('staying alive, ah ah ah staying aliiiiiveee')

                #calculate the pairwise distances for the current cluster of interest
                intraDists = pdist(clustCoordinates)

                #input the maximum intra cluster distance for the current cluster
                maxIntra[j] =  max(intraDists)

            elif isinstance(cluster, np.integer) or isinstance(cluster, int):
                #find center and place in dictionary
                center = np.zeros((1,num_groups))
                center[0,:] = dists[cluster]
                centersNum[j,:] = center

                #put the maximum intra distance as 0, only a single point (these will just be place holders)
                maxIntra[j] = 0


        #calculate the max intra cluster spread
        compactness = max(maxIntra)

        #calculate the pairwise distances between cluster centers
        sepDists = pdist(centersNum)

        sepDist = min(sepDists)
        
        if curClustersLength > 1:
            #calculate the validation index
            val_index[0,i] = compactness/sepDist
            val_index[1,i] = clusterings - (i)

        else:
            val_index[0,i] = 0
            val_index[1,i] = clusterings - (i)

        if i == 0:
            # logging.info(': 0% completed')
            print('first done')

        elif i%10==0:
            #calculate the amount of time the algorithm has been running 
            curTime = time.time()
            runTime = (curTime - initTime)
            # runTime = timeConverter(runTime)
            percent = 100*(1 - ((clusterings-i)/numIts))
            message1 = round(percent,2) 
            message1 = str(message1) + '% completed'
            # logging.info(message1)
            message = str(i)+': ' + str(runTime)
            # logging.info(message) 

        endTime = time.perf_counter()
        totalTime = endTime -startTime

    runTime = time.time()-initTime

    val_index = val_index[:,startPoint:clusterings]
    return val_index