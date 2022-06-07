import numpy as np
import time
import statistics as stat
from scipy.spatial.distance import pdist,squareform
import logging

class ValidationMetric:

    def daviesBouldin(data,dists,num_groups):
        '''
        Determine the appropriate number of clusters for the optimal clustering of a input data set. 
        
        Input:
        data - dictionary of the clusters for validation, or raw data for the interpretation of the appropriate number of clusters...
        dists - standardized or submitted non-standardized data. 
        num_groups - number of groups in the data set. 

        Output:

        Array of the number of clusters and the validation index measure. 
    
        '''

        #logging.info(': Starting cluster validation!')
        #grab the input dictionary size
        clusterings = len(data)

        startPoint = 0.5*clusterings
        startPoint = int(startPoint)
        numIts = clusterings - startPoint

        numClusters = clusterings-startPoint
        numClusters = int(numClusters)
        val_index = np.zeros((clusterings,2))
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


            for j in range(curClustersLength):
                #pull out the current cluster of metabolites
                cluster = curClusters[j]

                #current sum of intra cluster distances
                sumIntra = 0

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

                    #an array containing the current intra cluster comparison
                    curDistIntra = np.zeros((2,num_groups))

                    #calculate the intra cluster distance
                    for k in range(lengthList):
                        #calculate the pdist for each metabolite feature against the center
                        curMetab = clustCoordinates[k,:]
                        curDistIntra[0,:] = curMetab
                        curDistIntra[1,:] = center
                        sumIntra += pdist(curDistIntra)


                    #calculating the dispersion of the current cluster of interest
                    dispersion[j] = sumIntra/lengthList

                elif isinstance(cluster, np.integer) or isinstance(cluster, int):
                    #find center and place in dictionary
                    center = np.zeros((1,num_groups))
                    center[0,:] = dists[cluster]
                    centersNum[j,:] = center

                
            ##
            ##------------Calculating the R_i for the current subset of cluster--------------------
            ##
                
            #calculate the euclidean distances between the centers...
            cenDists = pdist(centersNum)
            cenDists = squareform(cenDists)

            #setting the maximum value of the R_i to zero initially, will test each iteration to determine the maximum 
            riMaxes = np.zeros((1,curClustersLength))
            for k in range(curClustersLength):
                curMax = 0
                for j in range(curClustersLength):
                    if k != j:
                        #calculate the R_i and compare to the curMax
                        R_i = (dispersion[k]+dispersion[j])/cenDists[k,j]
                        
                        if R_i > curMax:
                            curMax = R_i

                riMaxes[0,k] = curMax

            #sum up the R_i maxes and divide by the number of clusters in the current partition
            sumRi = np.sum(riMaxes)
            K = curClustersLength

            if K > 1:
                #calculate the validation index
                val_index[i,0] = sumRi/K
                val_index[i,1] = K

            else:
                val_index[i,0] = 10
                val_index[i,1] = K

        return val_index


    def dunnIndex(data,dists,num_groups):
        '''

        Dunn Index (1974)
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
        val_index = np.zeros((clusterings, 2))
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

            if curClustersLength > 1:
                sepDist = min(sepDists)
                #calculate the validation index
                val_index[i,0] = sepDist/compactness
                val_index[i,1] = curClustersLength

            else:
                val_index[i,0] = 0
                val_index[i,1] = curClustersLength

        return val_index


    def PBM(data,dists,num_groups,Eo):
        '''
        Determine the appropriate number of clusters for the optimal clustering of a input data set. 
        
        Input:
        data - dictionary of the clusters for validation.
        dists - standardized or submitted non-standardized data. 
        num_groups - number of groups in the data set. 

        Output:

        Array of the number of clusters and the validation index measure. 

        '''

        logging.info(': Starting cluster validation!')
        #grab the input dictionary size
        clusterings = len(data)

        startPoint = 0.5*clusterings
        startPoint = int(startPoint)

        numIts = clusterings - startPoint


        #number of clusters (K)
        K = len(data[0])
        K = int(K)

        val_index = np.zeros((clusterings,2))
        initTime = time.time()

        for i in range(startPoint,clusterings):
            #**********************************************************************************************************
            #**********************************Threading should occur here*********************************************
            #**********************************************************************************************************
            #**********************************************************************************************************
            
            startTime = time.perf_counter()
            #grab the current set of metabolite clusters
            curClusters = data[i]

            #from the current clusters determine the length in order to determine the next step
            curClustersLength = len(curClusters)

            #sum of intra cluster distances
            sumIntra = 0

            #create a numpy array for that contains the cluster centers for calculation of the inter cluster distance.
            centersNum = np.zeros((curClustersLength,num_groups))
            for j in range(curClustersLength):

                #pull out the current cluster of metabolites
                cluster = curClusters[j]
                
                #check for instances of the clusters imbedded in the dictionary for each clustering outcome
                if isinstance(cluster, list):
                    #check the length of the cluster and then pull in the values for each of the clustered metabolites
                    lengthList = len(cluster)
                    clustCoordinates = np.zeros((lengthList,num_groups))

                    for k in range(lengthList):
                        #grab the cluster coordinates for each metabolite clustered with this particular round of clusters
                        clustCoordinates[k,:] = dists[cluster[k]]

                    #Create a numpy array for the coordinates of the center of the current cluster
                    center = np.zeros((1,num_groups))

                    for m in range(num_groups):
                        #find the mean of each group in the cluster
                        center[0,m] = stat.mean(clustCoordinates[:,m])
                    #update dictionary of the cluster centers for later calculation of inter-cluster distances
                    centersNum[j,:] = center

                    #initialize an array that stores the values that will be sent to the pdist function
                    curDistIntra = np.zeros((2,num_groups))

                    #calculate the intra_cluster distance for each list of metabolites in the 
                    for k in range(lengthList):
                        #grab the first value from the list find its distance and add it to the sum of the Intra cluster distances
                        curMetab = clustCoordinates[k,:]
                        curDistIntra[0,:] = curMetab
                        curDistIntra[1,:] = center
                        sumIntra += pdist(curDistIntra)

                elif isinstance(cluster, np.integer) or isinstance(cluster,int):
                    #find the center and put it into the dictionary
                    center = np.zeros((1,num_groups))
                    center[0,:] = dists[cluster]
                    centersNum[j,:] = center

            # find the distance between the centers
            centerDists = pdist(centersNum)

            
            if curClustersLength > 1:
                D = max(centerDists)

                PBM = ((D/K)*(Eo/sumIntra))**2
            
                val_index[0,0]=PBM
                val_index[0,1]=K
            else:
                val_index[0,0]=0
                val_index[0,1]=K

        return val_index

    def Silhouette(data,dists,num_groups):
        '''
        Determine the appropriate number of clusters for the optimal clustering of a input data set. 
        
        Input:
        data - dictionary of the clusters for validation.
        dists - standardized or submitted non-standardized data. 
        num_groups - number of groups in the data set. 

        Output:

        Array of the number of clusters and the validation index measure. 

        '''

        logging.info(': Starting cluster validation!')
        #grab the input dictionary size
        clusterings = len(data)

        startPoint = 0.5*clusterings
        startPoint = int(startPoint)

        numIts = clusterings - startPoint


        #number of clusters (K)
        K = len(data[0])
        K = int(K)

        val_index = np.zeros((clusterings,2))
        initTime = time.time()

        for i in range(startPoint,clusterings):
            #**********************************************************************************************************
            #**********************************Threading should occur here*********************************************
            #**********************************************************************************************************
            #**********************************************************************************************************
            
            startTime = time.perf_counter()
            #grab the current set of metabolite clusters
            curClusters = data[i]

            #from the current clusters determine the length in order to determine the next step
            curClustersLength = len(curClusters)

            #sum of intra cluster distances
            sumIntra = 0
            
            #create a numpy array for that contains the cluster centers for calculation of the inter cluster distance.
            coordinates = {}
            for j in range(curClustersLength):
                #get the current cluster
                cluster = curClusters[j]

                #find all of the cluster centers and then advance to finding the silhouette Index
                if isinstance(cluster,list):
                    #determine number of metabolites
                    lengthList = len(cluster)

                    #find metabolite coordinates for eventual determination of center
                    clustCoordinates = np.zeros((lengthList,num_groups))
                    for k in range(lengthList):
                        #get coordinates
                        clustCoordinates[k,:] = dists[cluster[k]]


                    #put clust coordinate into the dictionary.
                    coordinates[j] = clustCoordinates

                elif isinstance(cluster,np.integer) or isinstance(cluster,int):
                    clustCoordinates = np.zeros((1,num_groups))
                    clustCoordinates[0,:] = dists[cluster]
                    coordinates[j] = clustCoordinates
            

            #loop thru each data feature within a cluster, and first get average intra-distance.
            SIL = 0
            for j in range(curClustersLength):
                cluster = coordinates[j]
                clusterSIL = 0
                #loop thru the points in cluster k, finding the average intra cluster separation...
                for k in range(cluster.shape[0]):
                    
                    if cluster.shape[0] > 1:
                        #for each value put take numpy array and find the distances
                        curDists = pdist(cluster)
                        #put into squarefomr for ease of use.
                        curDists = squareform(curDists)

                        #sum the row of interest k
                        sumRow = np.sum(curDists[k,:])
                        #compute the average intracluster separation
                        a_x = sumRow/(cluster.shape[0]-1)

                        #number of comparisons is numClustersLength - 1 for calculating b_x
                        curMax = 0
                        #setting b_x prior to its updated assignment to ensure no fail out once all data points are clustered together.
                        b_x = 0
                        for m in range(curClustersLength):
                            if k != m:
                                #get the numpy array of other cluster, and add the data feature of interest from the other cluster to the top of the numpy array (i.e., row 0 for consistency)
                                curArr = coordinates[m]

                                #put cluster coordinate k on top row of curArr
                                curArr = np.vstack([cluster[k,:], curArr])

                                #get the current distances between the data feature of interest and the other cluster data features
                                curInterDists = pdist(curArr)
                                curInterDists = squareform(curInterDists)

                                #get the top row of distances 
                                interDistsAvg = curInterDists[0,:]

                                b_x = np.sum(interDistsAvg)/(curArr.shape[0]-1)

                        #determine the max between a_x and b_x
                        denominator = max([a_x,b_x])

                        #calculate the silhouette for the current data feature.
                        curSIL = (b_x-a_x)/denominator

                        clusterSIL += (curSIL/cluster.shape[0])

                    else:
                        clusterSIL += 0
                
                #add cluster SIL to whole SIL
                SIL += (clusterSIL/K)
                    
            #put into the format for graphing later.
            val_index[0,0]=SIL
            val_index[0,1]=K

        return val_index