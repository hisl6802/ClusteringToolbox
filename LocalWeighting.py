from numpy.lib.arraysetops import isin
from sklearn.impute import SimpleImputer
import GuiBackground as GB
import GUIUtils as GU
import sys, logging
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist,squareform
import math
import numpy as np

class LocallyWeighted:
 
#log that the user called ensemble clustering function
#logging.info(': User called Ensemble Clustering function.')

#optimum number of clusters from validation index.
# sys.setrecursionlimit(10**8)

#Make sure file can be read in. 
# metab_data = GB.fileCheck()
# if metab_data is None:
#     logging.error(': File to not meet input requirements.')

#List for the use in creating and plotting the clustering results
# linkageList = ['single','complete','average']
# distList = ['euclidean','sqeuclidean','chebyshev','seuclidean']  #,'cosine']

#calculate the number of clusterings based upon the size of the lists and an additional term for the ward-euclidean run. 
# numClusterings = (len(linkageList)*len(distList))+1

#read in the data
# data = GB.readInColumns(metab_data)

#Standardize the data before clustering the results
#logging.info(': Standardizing data.')
# for i in range(data.shape[0]):
#     data[i,:] = GB.standardize(data[i,:])

#creates empty dictionary for clusterings
# clusters = {}

#performs first 12 base clusterings and populates clusters dictionary
# for i in range(len(linkageList)):
#      for j in range(len(distList)):
#          linkCur = linkage(data,linkageList[i],distList[j])
#          valid = GB.clustConnectLink(linkCur)
#          index = str(linkageList[i] + '_' + distList[j])
#          clusters.update({index:valid})
#          logging.info(str(linkageList[i])+'-'+str(distList[j]) +' done!')

#performs 13th base clustering and populates clusters dictionary
# linkCur = linkage(data, 'ward', 'euclidean')
# valid = GB.clustConnectLink(linkCur)
# logging.info(str('ward-euclidean done!'))
# clusters.update({'ward_euclidean':valid})


    def clustFinder(data=0,optNum=2,clusters = {}):
        '''
        Used to locate the clusters of interest.

        Input: Optimum number of clusters as determined from MST. Same value as used in ensembleClustering function.

        Output: Object containing the clusters of interest.
        '''
        print(data)
        numMetabs = data.shape[0]
        dictLoc = numMetabs - optNum - 1
        
        refClust = {}

        for key in clusters:
            refClust.update({key:clusters[key][dictLoc]})


        return refClust#, optNum 


    # refClust = clustFinder()


    def probClust(keys, key, localClust, refClust):
        '''
        
        '''
        #Brady 08/10/21
        keys.remove(key)

        #checks if the local cluster is a list or not. 
        if isinstance(localClust, list):
            #local cluster is a list so we need to loop over the list.
            #loop over the length of the localClust list
            probNums = {}
            for i in keys:
                curCompClusts = refClust[i]
                nonZeroComps = []
                for key in curCompClusts:
                    curClust = curCompClusts[key]
                    tally = 0
                    if isinstance(curClust, list):
                        for j in range(len(localClust)):
                            if localClust[j] in curClust:
                                tally += 1
                    else:
                        if curClust in localClust:
                            tally += 1

                    if tally > 0:
                        nonZeroComps.append(tally)
                
                #put the list into the appropriate 
                probNums[i] = nonZeroComps

        else:
            probNums = {}
            for i in keys:
                curCompClusts = refClust[i]
                nonZeroComps = []
                for key in curCompClusts:
                    curClust = curCompClusts[key]
                    tally = 0 
                    if isinstance(curClust,list):
                        if localClust in curClust:
                            tally += 1
                    else:
                        if localClust == curClust:
                            tally += 1

                    if tally > 0:
                        nonZeroComps.append(tally)
                probNums[i] = nonZeroComps
        
        return probNums

            


    def clustCompare(refClust):
        '''
        Used to compute the proportion of elements in local cluster i shared with elements of local cluster j within each base clustering. 

        Input: Object containing the clusters of interest as determined via clustFinder().

        Output: p(C_i, (C_j)^m)
        '''

        refClust = clustFinder()
        ECI = {}

        #populating the ECI dictionary with the keys from refClust dictionary, which point to a dictionary with the optimal number of keys (i.e., number of ECI for each) 
        for key in refClust:
            ######
            ###### Need to update to accomadate different number of optimal number of clusters
            ######
            ECI.update({key:{0:0,1:1}})


        for key in refClust:
            baseClust = refClust[key]

            for i in baseClust:
                keys = list(refClust.keys())
                localClust = baseClust[i]
                if isinstance(localClust, list):
                    denom = len(localClust)
                else: 
                    denom = 1

                #print the denominator prior to sending it to the function.
                #print("Probability being calculated...")
                p = probClust(keys, key, localClust, refClust) 
                entropyKeys = list(p.keys())
                entropy = 0
                for k in entropyKeys:
                    #grab the list from the dictionary to calculate the entropy
                    curEntropy = p[k]
                    for j in range(len(curEntropy)):
                        #calculate the entropy given the list
                        entropy += -(curEntropy[j]/denom)*math.log(curEntropy[j]/denom)
                        

                ECI[key][i] = math.exp(-entropy/(0.5*13))
        
        return ECI

    # ECI = clustCompare(refClust) 

    def consensus(ECI, refClust):
        '''
        Consensus matrix calculation

        Inputs:
        ECI: clustering index (weight of summation)
        refClust: dictionary containing the clusters of each base clustering which each metabolite/entity belongs

        Output:
        Consensus matrix
        '''

        ####
        #### Currently setting the M to 13 based upon 13 base clusterings being performed.
        ####
        M = 13

        #create consensus matrix of zeros to start
        numMetabs = data.shape[0]
        consensusMat = np.zeros((numMetabs,numMetabs))


        keys = list(refClust.keys())
        #popultate the diagonals with ones, since the evidence should always suggest that 
        for i in range(numMetabs):
            consensusMat[i,i] = 1

        #loop over the number of metabolites to determine the appropriate value for each matrix input
        for i in range(numMetabs):
            for j in range(i+1,numMetabs):
                curSum = 0
                for key in refClust:
                    #locate the cluster which contains the ith metabolite, in the current base clustering
                    for k in range(len(refClust[key])):
                        #search each list. 
                        curClust = refClust[key][k]
                        if isinstance(curClust, list):
                            #check for entity in list
                            if i in curClust and j in curClust:
                                curSum += ECI[key][k]
                                continue

                #populate the i,j position and the j,i position
                consensusMat[i,j] = curSum/M
                consensusMat[j,i] = curSum/M

        
        return consensusMat



    #determine the consensus
    # consensusMat = consensus(ECI,refClust)


    def regions(similarity):
        '''
        Calculating the regions/clusters

        Input:
        Consensus matrix - similarity matrix

        Output final matrix of regions. 
        '''

        #set up regions dictionary to contain all regions/clusters which are in the current clustering
        regions = {}


        #trick computer to think that the same elements have no similarity, by making the diagonal elements zeros. 
        for i in range(similarity.shape[0]):
            similarity[i,i] = 0

        #create matrix matching similarity that is current similarity matrix. 
        curSim = similarity

        #loop over N-1 times to determine regions each time
        for i in range(similarity.shape[0]-1):
            #fill dictionary with N-1 keys, with the key matching teh values at first, before being updated each iteration
            regions[i] = {}

            #fill first entry with a dictionary containing integers in each location. 
            if i == 0:
                for j in range(similarity.shape[0]):
                    #fill with initial region with dictionary of metabolites
                    regions[i][j] = j

            else:
                #search the upper matrix for the max value (i.e., highest similarity)
                curMax = np.amax(curSim)
                curMaxLoc = np.where(curSim == curMax)

                #for the indicies combine the values from the locations in the dictionary into the first entry to the next dictionary
                #grab teh first key value/object, and the second key value/object
                first = regions[i-1][curMaxLoc[0][0]]
                second = regions[i-1][curMaxLoc[1][0]]

                #check first, and second for type of variable, list or integer
                if isinstance(first,list):
                    if isinstance(second,list):
                        #combine lists
                        regions[i][0] = first + second

                    else:
                        #append second(which is an integer) to list
                        regions[i][0] = first + [second]

                elif isinstance(second,list):
                    #append first (which is an integer) to list
                    regions[i][0] = second + [first]
                
                else:
                    #combine integers into a list
                    regions[i][0] = [first, second]
                
                #updating to contain full list for when the first entry doesn't match
                for j in range(0,len(regions[i-1])):
                    #input all previous dictionary entries unless those indicies match those seen in 
                    if j != curMaxLoc[0][0]:
                        if j != curMaxLoc[1][0]:
                            regions[i][len(regions[i])] = regions[i-1][j]
                    elif j != curMaxLoc[1][0]:
                        if j != curMaxLoc[0][0]:
                            regions[i][len(regions[i])] = regions[i-1][j]
                

                #grab shape of current similarity matrix
                numRows = curSim.shape[0]
                #overwrite curSim
                curSimNew = np.zeros((numRows-1,numRows-1))

                #calculate the new similarity matrix between the new region and the existing regions. 
                for j in range(len(regions[i])):
                    for k in range(j+1,len(regions[i])):
                        #get the current regions of interest
                        region1 = regions[i][j]
                        region2 = regions[i][k]

                        if isinstance(region1,list):
                            if isinstance(region2,list):
                                #find the length of each list
                                len1 = len(region1)
                                len2 = len(region2)
                                curSum = 0

                                for l in range(len(region1)):
                                    for m in range(len(region2)):
                                        curSum += similarity[region1[l],region2[m]]
                                
                                #input the new similarity
                                curSimNew[j,k] = curSum/(len1*len2)
                                curSimNew[k,j] = curSum/(len1*len2)
                            else:
                                #find the length of the first region, second region has length 1
                                len1 = len(region1)
                                len2 = 1
                                curSum = 0

                                for l in range(len(region1)):
                                    curSum += similarity[region1[l],region2]

                                #input the new similarity
                                curSimNew[j,k] = curSum/(len1*len2)
                                curSimNew[k,j] = curSum/(len1*len2)

                        elif isinstance(region2,list):
                            #find the length of region2, region1 has length 1
                            len1 = 1
                            len2 = len(region2)
                            curSum = 0
                            
                            for l in range(len(region2)):
                                curSum += similarity[region1,region2[l]]

                            curSimNew[j,k] = curSum/(len1*len2)
                            curSimNew[k,j] = curSum/(len1*len2)

                        else:
                            curSimNew[j,k] = similarity[region1,region2]
                            curSimNew[k,j] = similarity[region2,region2]
                curSim = curSimNew

        return regions


    # regionsOut = regions(consensusMat)

