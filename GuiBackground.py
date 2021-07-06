import numpy as np
import statistics as stat
from scipy.cluster.hierarchy import dendrogram
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist,squareform
from scipy.sparse import csr_matrix
from matplotlib import pyplot as plt
import pandas as pd
import seaborn as sns
from tkinter import filedialog
import logging, time, glob,sys,os,ast
from PIL import Image

def fileCheck(file=''):
    '''
    Check that the selected file is of the appropriate file extension and is able to be read in. 

    Input:

    Does not require an input, but does accept a full file path which can be used to check for an excel sheet with 'Medians' as 
    the sheet name. 

    Output:
    
    Returns the data from the excel file. 
    '''
    #log that the user called the Create Clustergram function
    logging.warning(': Checking file for the appropriate input.')
    #ask the user to select the file that they would like to create a clustergram for.
    if file =='':
        file = filedialog.askopenfilename()
    
    if file == '': 
        logging.error(': Failed to select a file')
        return

    #Open the excel file that the user to looking to use for their clustering analysis
    try:
        data = pd.read_excel(file)
        #data = pd.read_excel(file, sheet_name='Medians')
        del(file)
    except:
        logging.error(': Failed to read in the excel sheet check the sheet name is Medians')
        return
    logging.info(': User file opened and submitted to the function.')

    return data
        
#Standardizing the data that is input to python.
def standardize(data):
    '''
    Standardize the input data for best clustering results.
    
    Input:

    Single row of data.

    Output:

    Standardized single row of data. 

    '''
    #I would like to add the ability to check for the number of rows or columns that are given in a numpy array.
    #find the mean and standard deviation of the given row(eventually will need to move this to allow for the entire table to input.)
    mean_data = stat.mean(data)
    std_data = stat.stdev(data)
    dataCur = np.zeros(data.shape[0])

    if std_data == 0:
        for i in range(data.shape[0]):
            dataCur[i] = 0
    else:
        for i in range(data.shape[0]):
            dataCur[i] = (data[i]-mean_data)/std_data

    return dataCur

#standardizing data after it has been normalized to a specific biological profile
def normStandardize(data,leaves):
    '''
    Normalize input data that has been standardized for normalized clustergrams.

    Input:

    data - all data. 
    leaves - leaves of the groups you would like to normalize. 

    ***Currently the only possible normalization is the first column of data. Will be updated soon. 

    Output:
    
    Normlized standardardized data. 

    '''
    logging.info(': Normalizing strandardized data.')
    #The mean for a normalized data set should always be 1 for all of the metabolites.
    mean = data[leaves[0]]

    #initialize the needed arrays
    dataCur = np.zeros(data.shape[0])

    #Calculate the residuals of the data
    residuals = 0
    for j in leaves:
        residuals += (data[j]-mean)**2
    #calculate the standard deviation
    std_data = (residuals/(data.shape[0]-1))**0.5
    #standardize each row of data
    for i in leaves:
        #Input the data into an array of standardized values
        dataCur[i] = (data[i]-mean)/std_data

    return dataCur

#ensemble dendrogram function
def createEnsemDendrogram(data,metab_data,norm=1,link='ward',dist='euclidean',func="ensemble"):
    '''
    Create ensemble dendrogram
    
    Input:

    Required:

    data - standardized data
    metab_data - original raw data

    Optional:

    norm - 0 or 1 for the normalization of ensemble (should always be zero)
    link - input the linkage function you would like the program to use for the ensemble clustering (or final clustering)
    dist - distance measure you would like to use in the final clustering of the data. 
    func - should always be ensemble. 


    Output:
    
    Outputs the ensemble clustergram and sends the output to the function which generates the output files. 

    '''
    #generate linkage function
    linkageMetabOut = linkage(data,link,dist)
    #create the dendrogram
    metaboliteDend = dendrogram(linkageMetabOut, orientation='left',no_plot=True)
    metaboliteDendLeaves = metaboliteDend['leaves']

    #tranpose the data and then run an analysis on the groups.
    groupCluster = np.transpose(data)
    logging.info(': Ensemble clustergram being generated.')

    #Create a linkage matrix for the data
    linkageGroupOut = linkage(groupCluster,link,dist)
    #create the dendrogram of the output
    groupDend = dendrogram(linkageGroupOut,orientation='top',no_plot=True)
    #get the leaves of the dendrogram to allow for the appropriate regrouping of the data
    groupDendLeaves = groupDend['leaves']

    #Rework the data to create the clustergram
    dataFinal = np.zeros((data.shape[0],data.shape[1]))
    for i in range(data.shape[1]):
        #rearranging the data for heatmap
        for j in range(data.shape[0]):
            #going through the metabolites
            dataFinal[j,i] = data[metaboliteDendLeaves[j],groupDendLeaves[i]]
    g = sns.clustermap(data, figsize=(7, 5), yticklabels=False, xticklabels=False, row_linkage=linkageMetabOut, col_linkage=linkageGroupOut, cmap="viridis", cbar_pos=(0.01, 0.8, 0.025, 0.175))
    ax = g.ax_heatmap

    
    recClusters(dataFinal,ax,groupDendLeaves,metab_data)
    plt.savefig('Clustergram01.png',dpi=600)
    plt.show()

#dendrogram function
def create_dendrogram(data, norm=1,link='ward',dist='euclidean'):
    '''
    Create dendrogram for either the ensemble or clustergram functions

    Input:

    Required:

    data - standardized data recommended. 

    
    Optional:
    norm - 0 or 1 (0 is for standard clustergram, 1 for normalized clustergram)
    link - input string of wanted linkage function. 
    dist - input string of wanted distance measure.
    func - 'clustergram' -> do not change the argument. 

    Output:

    This function outputs a clustergram using the seaborn clustermap function 
    '''
    #increase recursion limit
    sys.setrecursionlimit(10**8)
    
    #Create the linkage matrix
    linkageMetabOut = linkage(data,link,dist)

    #tranpose the data and then run an analysis on the groups.
    try:
        groupCluster = np.transpose(data)
        logging.info(': Clustergram being generated.')
    except:
        logging.error(': Inappropriate entry for the create_dendrogram function. Check docs.')

    #Create a linkage matrix for the data
    linkageGroupOut = linkage(groupCluster,link,dist)

    if norm == 0:
        g = sns.clustermap(data, figsize=(7, 5), yticklabels=False, row_linkage=linkageMetabOut, col_linkage=linkageGroupOut, cmap="viridis", cbar_pos=(0.01, 0.8, 0.025, 0.175))
        plt.show()

    elif norm == 1:
        for i in range(dataFinal.shape[0]):
            data[i,:] = normStandardize(data[i,:],groupDendLeaves)
        #create array that will store final version of reorganized data.
        dataFinalNorm = np.zeros((data.shape[0],data.shape[1]))

        for i in range(data.shape[1]):
            #rearranging the data for heatmap
            for j in range(data.shape[0]):
                #going through the metabolites
                dataFinalNorm[j,i] = data[metaboliteDendLeaves[j],groupDendLeaves[i]]
         
        #create the axes in which the heatmap will be mapped upon
        heatmapAxes = [0.15, 0.05, 0.8, 0.8]
        heatmapAxes = fig.add_axes(heatmapAxes)  
        #output the normalized heatmap 
        heatmapAxes.matshow(dataFinalNorm,aspect='auto',origin='upper',cmap="hot")
        plotting(link=link,dist=dist)

def cooccurrence(data):
    '''
    Creation of the cooccurrence matrix for the determination of the number times each set of metabolites is clustered together
    in a set of N clusterings. 
    
    Input:
    data - standardized data or data that in general you would like to have an ensemble clustering performed on. 

    Output:

    NxN numpy array of zeros for the initial cooccurrence matrix. 

    '''
    logging.info(': Creating the co-occurrence matrix.')

    #find the number of metabolites
    numMetabolites = data.shape[0]

    #create co-occurrence matrix of zeros
    coOcc = np.zeros((numMetabolites,numMetabolites))

    for i in range(numMetabolites):
        #populate the diagonal matrix with ones
        coOcc[i,i] = 1

    return coOcc

def popCooccurrence(clusters,coOcc,numClusterings):
    '''
    Populate the cooccurrence matrix with the connections for each matrix based upon the number of clusterings that occur. 
    
    Input:

    clusters - dictionary containing the clusters from clustering (link-dist)
    coOcc - pass the current version of the cooccurence matrix. 
    numClusterings - pass integer of the optimal number clusters. 

    Output:
    
    Updated NxN cooccurence matrix. 
    '''
    logging.info(': Populating the co-occurrence matrix.')
    dictKeys = list(clusters.keys())
    dictKeys = len(dictKeys)

    #weight to be added to the matrix
    weight = 1/numClusterings

    #populate the coOccurence matrix with the appropriate values based upon the keys using the binary operation +=
    for i in range(dictKeys):
        #get the first entry from the dictionary 
        curCluster = clusters[i]

        #check for a list of metabolites clustered together.
        if isinstance(curCluster, list):
            #determine the length of the list
            lenList = len(curCluster)
            
            for j in range(lenList):
                #populate the Cooccurrence matrix
                remList = lenList-j
                if remList > 1:
                    for k in range(j+1,lenList):
                        #add the weight to the coOcc matrix. 
                        coOcc[curCluster[j],curCluster[k]] += weight
                        coOcc[curCluster[k],curCluster[j]] += weight
    return coOcc

def clustConnectLink(linkageCheck):
    '''
    Determine the connections from the scipy linkage function output.

    Input:

    linkageCheck - a scipy linkage output 

    Output:

    validationClusters, the metabolites (other identity) output into a dictionary with M keys (N(metabolites)-1) 

    '''
    logging.info(': Starting to determine metabolite clusters.')
    #determine the clusters using a new alogorithm based on linkage method.
    #create a dictonary to store the linkage functions. 
    metabolites = np.zeros((linkageCheck.shape[0]+1,2))
    #fill the array with the metabolite identifiers
    for i in range(linkageCheck.shape[0]+1):
        metabolites[i,0] = i
        metabolites[i,1] = i

    metabolites = metabolites.astype(int)
    linkageCheck = linkageCheck.astype(int)
    #set the limit for the metabolites such that if the metabolites match up with the 
    metabLimit = linkageCheck.shape[0]
    limit = linkageCheck.shape[0]

    #create dictionary to store the clusters as they are created.
    clusters = {}

    #create an empty dictionary to store all of the cluster configurations.
    validationClusters = {}
    for i in range(linkageCheck.shape[0]):
        #check the linkage connections to determine the appropriate 
        curCon1 = linkageCheck[i,0]
        curCon2 = linkageCheck[i,1]

        curCon1 = int(curCon1)
        curCon2 = int(curCon2)

        if i == 0:
            #place first connection into the list 
            firstConnect = [curCon1, curCon2]
            metabolites[curCon1,1] = limit + 1
            metabolites[curCon2,1] = limit + 1
            limit += 2

            for j in range(linkageCheck.shape[0]-(i+1)):
                clusters.update({j:j})

            #populate the first part of the dictionary with the first list containing the clustered metabolites
            clusters[0] = firstConnect 

            #search the array for values of the array which are equal to the one plus the number of metabolites studied
            unClustered = np.where(metabolites[:,1] != (metabLimit+1))
            for j in range(1,linkageCheck.shape[0]):
                clusters[j] = unClustered[0][j-1]
  
        else:
            #save the previous clusters dictionary prior to deleting it. 
            clusterPrevious = clusters
            del(clusters)
            clusters = {}

            for j in range(linkageCheck.shape[0]-(i+1)):
                clusters.update({j:j})

            curCon1 = linkageCheck[i,0]
            curCon2 = linkageCheck[i,1]

            curCon1 = int(curCon1)
            curCon2 = int(curCon2)

            if curCon1 > metabLimit and curCon2 <= metabLimit:
                #search the second column of the list for the appropriate value.
                oneCon = np.where(metabolites[:,1] == curCon1)
                connect1 = oneCon[0][0]
                curCon1 = int(connect1)
                curCon2 = int(curCon2)
            elif curCon1 > metabLimit and curCon2 > metabLimit:
                #search the second column of the reference list for the appropriate value.
                oneCon = np.where(metabolites[:,1] == curCon1)
                connect1 = oneCon[0][0]
                curCon1 = int(connect1)
                twoCon = np.where(metabolites[:,1] == curCon2)
                connect2 = twoCon[0][0]
                curCon2 = int(connect2)
            elif curCon1 <= metabLimit and curCon2 > metabLimit:
                #search the second column of the list for the appropriate value. 
                twoCon = np.where(metabolites[:,1] == curCon2)
                connect2 = twoCon[0][0]
                curCon2 = int(connect2)
                curCon1 = int(curCon1)

            curCon1Connect = 0;
            curCon2Connect = 0;
            unchanged = []

            previousKey = list(clusterPrevious.keys())
            previousKey = len(previousKey)

            for k in range(previousKey):
                #Determine matching metabolites in the dictionary.
                curCheck = clusterPrevious[k]

                if isinstance(curCheck, list):
                    #check list for the first connection
                    curCon1Check = curCon1 in curCheck

                    #check list for the second connection
                    curCon2Check = curCon2 in curCheck

                    if curCon1Check == True and curCon2Check == False:
                        curCon1Connect = 1
                        curCon1Location = k
                    elif curCon1Check == False and curCon2Check == True:
                        curCon2Connect = 1
                        curCon2Location = k
                    elif curCon1Check == False and curCon2Check == False:
                        unchanged.append(k)
                    elif curCon1Check == True and curCon2Check == True:
                        logging.warning(': Issue clustering the data a duplication has been discovered.')


                elif isinstance(curCheck, np.integer) or isinstance(curCheck, int):
                    #check list for the first connection
                    curCon1Check = curCon1 == curCheck

                    #check list for the second connection
                    curCon2Check = curCon2 == curCheck

                    if curCon1Check == True and curCon2Check == False:
                        curCon1Connect = 1
                        curCon1Location = k
                    elif curCon1Check == False and curCon2Check == True:
                        curCon2Connect = 1
                        curCon2Location = k
                    elif curCon1Check == False and curCon2Check == False:
                        unchanged.append(k)
                    elif curCon1Check == True and curCon2Check == True:
                        loggging.warning(': Issue clustering the data a duplication has been discovered.')

                else:
                    logCheck = type(curCheck)
                    logging.error(logCheck)
                    logging.error(': Inappropriate data type for the for clustering ID algorithm.')
                    return

                if curCon1Connect == 1 and curCon2Connect == 1 and k == previousKey-1:
                    for m in range(1,len(unchanged)+1):
                        #cluster the appropriate remaining clusters together
                        clusters[m] = clusterPrevious[unchanged[m-1]]

                    newConnect1 = clusterPrevious[curCon1Location]
                    newConnect2 = clusterPrevious[curCon2Location]
                    
                    if isinstance(newConnect1,list) and isinstance(newConnect2, list):
                        newConnect = newConnect1 + newConnect2
                        clusters[0] = newConnect
                        for m in newConnect:
                            metabolites[m,1] = limit
                        limit += 1
                    elif isinstance(newConnect1,list) and isinstance(newConnect2, np.integer):
                        newConnect = newConnect1
                        intList = newConnect[:] + [newConnect2]
                        clusters[0] = intList
                        for m in intList:
                            metabolites[m,1] = limit
                        limit += 1
                    elif isinstance(newConnect1,list) and isinstance(newConnect2, int):
                        newConnect = newConnect1
                        intList = newConnect[:] + [newConnect2]
                        clusters[0] = intList
                        for m in intList:
                            metabolites[m,1] = limit
                        limit += 1
                    elif isinstance(newConnect1,np.integer) and isinstance(newConnect2, list):
                        newConnect = newConnect2
                        intList = newConnect[:] + [newConnect1]
                        clusters[0] = intList
                        for m in intList:
                            metabolites[m,1] = limit
                        limit += 1
                    elif isinstance(newConnect1,int) and isinstance(newConnect2, list):
                        newConnect = newConnect2
                        intList = newConnect[:] + [newConnect1]
                        clusters[0] = intList
                        for m in intList:
                            metabolites[m,1] = limit
                        limit += 1
                    elif isinstance(newConnect1,np.integer) and isinstance(newConnect2,np.integer):
                        newConnect = [newConnect1, newConnect2]
                        clusters[0] = newConnect
                        for m in newConnect:
                            metabolites[m,1] = limit
                        limit += 1
                    elif isinstance(newConnect1,int) and isinstance(newConnect2,np.integer):
                        newConnect = [newConnect1, newConnect2]
                        clusters[0] = newConnect
                        for m in newConnect:
                            metabolites[m,1] = limit
                        limit += 1
                    elif isinstance(newConnect1,np.integer) and isinstance(newConnect2,np.integer):
                        newConnect = [newConnect1, newConnect2]
                        clusters[0] = newConnect
                        for m in newConnect:
                            metabolites[m,1] = limit
                        limit += 1
                    elif isinstance(newConnect1,int) and isinstance(newConnect2,int):
                        newConnect = [newConnect1, newConnect2]
                        clusters[0] = newConnect
                        for m in newConnect:
                            metabolites[m,1] = limit
                        limit += 1

        validationClusters.update({i:clusters})
    logging.info(': Success! Metabolite clusters determined.')
    
    return validationClusters

def clustConnect(dataMST,mstOutNp):
    '''
    Clustering connections determination from the minimum spanning tree output. 

    Input:
    dataMST - direct output from the MST
    mstOutNp - numpy array of the output. 

    Output:
    dictionary of all possible clusterings combinations for the data. 

    '''
    logging.info(': Determining the metabolite cluster for MST.')
    #Create an empty dictionary that will contain the clusters
    clusters = {}

    #create an empty dictionary that allows us to save the clusters for validation.
    validationClusters = {}

    # create initial list of metabolites that serves as the initial list of metabolites that will be clustered.
    metabolites = np.ones((dataMST.shape[0]+1,1))

    for i in range(dataMST.shape[0]):
        #pull out the connections for the current iteration
        curCon1 = mstOutNp[i,0]
        curCon2 = mstOutNp[i,1]

        curCon1 = int(curCon1)
        curCon2 = int(curCon2)

        if i == 0:
            #convert the metabolites to string for easier comparison
            firstConnect = [curCon1, curCon2]

            #set the metabolites equal to zero in the initial metabolite list
            metabolites[curCon1,0] = 0 
            metabolites[curCon2,0] = 0

            #create dictionary of clusters 
            for j in range(dataMST.shape[0]-(i+1)):
                clusters.update({j:j})
            
            #find the metabolites that are ones and were not clustered
            unClustered = np.where(metabolites == 1)

            #input the connection
            clusters[0] = firstConnect
            
            for j in range(1,dataMST.shape[0]):
                #input the cluster values into the dictionary
                clusters[j] = unClustered[0][j-1]
        else:
            #save the previous dictionary 
            clusterPrevious = clusters
            del(clusters)
            clusters = {}

            #create a new dictionary 
            for j in range(dataMST.shape[0]-(i+1)):
                clusters.update({j:j})

            #grab the latest connections
            curCon1 = mstOutNp[i,0]
            curCon2 = mstOutNp[i,1]

            curCon1 = int(curCon1)
            curCon2 = int(curCon2)

            curCon1Connect = 0;
            curCon2Connect = 0;
            unchanged = []

            previousKey = list(clusterPrevious.keys())
            previousKey = len(previousKey)

            for k in range(previousKey):
                #Determine if any of the new clustered meatbolites 
                curCheck = clusterPrevious[k]

                if isinstance(curCheck, list):
                    #check list for the first connection
                    curCon1Check = curCon1 in curCheck

                    #check list for the second connection
                    curCon2Check = curCon2 in curCheck

                    if curCon1Check == True and curCon2Check == False:
                        curCon1Connect = 1
                        curCon1Location = k
                    elif curCon1Check == False and curCon2Check == True:
                        curCon2Connect = 1
                        curCon2Location = k
                    elif curCon1Check == False and curCon2Check == False:
                        unchanged.append(k)
                    elif curCon1Check == True and curCon2Check == True:
                        logging.warning(': Issue clustering the data a duplication has been discovered.')


                elif isinstance(curCheck, np.integer) or isinstance(curCheck, int):
                    #check list for the first connection
                    curCon1Check = curCon1 == curCheck

                    #check list for the second connection
                    curCon2Check = curCon2 == curCheck

                    if curCon1Check == True and curCon2Check == False:
                        curCon1Connect = 1
                        curCon1Location = k
                    elif curCon1Check == False and curCon2Check == True:
                        curCon2Connect = 1
                        curCon2Location = k
                    elif curCon1Check == False and curCon2Check == False:
                        unchanged.append(k)
                    elif curCon1Check == True and curCon2Check == True:
                        loggging.warning(': Issue clustering the data a duplication has been discovered.')


                else:
                    logCheck = type(curCheck)
                    logging.error(logCheck)
                    logging.error(': Inappropriate data type for the for clustering ID algorithm.')
                    return

                if curCon1Connect == 1 and curCon2Connect == 1 and k == previousKey-1:
                    for m in range(1,len(unchanged)+1):
                        #cluster the appropriate remaining clusters together
                        clusters[m] = clusterPrevious[unchanged[m-1]]

                    newConnect1 = clusterPrevious[curCon1Location]
                    newConnect2 = clusterPrevious[curCon2Location]
                    
                    if isinstance(newConnect1,list) and isinstance(newConnect2, list):
                        newConnect = newConnect1 + newConnect2
                        clusters[0] = newConnect
                    elif isinstance(newConnect1,list) and isinstance(newConnect2, np.integer):
                        newConnect = newConnect1
                        intList = newConnect[:] + [newConnect2]
                        clusters[0] = intList
                    elif isinstance(newConnect1,np.integer) and isinstance(newConnect2, list):
                        newConnect = newConnect2
                        intList = newConnect[:] + [newConnect1]
                        clusters[0] = intList
                    elif isinstance(newConnect1,np.integer) and isinstance(newConnect2,np.integer):
                        newConnect = [newConnect1, newConnect2]
                        clusters[0] =newConnect

        validationClusters.update({i:clusters})
    logging.info(': Success! MST clusters determined.')
    return validationClusters
    
#initialize the plot
def plotting(link='',dist=''):
    '''
    Plotting the clustergram from the above create_dendrogram function. 

    Input:

    link - list of the linkage functions used to cluster your data. ex. createClustergram calls plotting after performing a ward-euclidean clustering link accepts the ward argument in a string. 
    dist - list of the distance measures used to cluster your data. See link for example. 

    Output:
    
    plots the clustegram. 

    '''
    logging.info(': Plotting Clustergram!')
    plt.xlabel('Clustered Metabolites')
    logging.info(': Saving...')

    sep = ''
    #check length of the link and dist inputs
    if len(link)>0 and len(dist)>0:
        #make the separator _
        sep += '_'

    clustPre = 'Clustergram'
    clustSuf = '.png'
    firstCheck = clustPre+sep+link+sep+dist+sep+'01'+clustSuf

    chkBuffer = glob.glob("*.png")
    count = 1
    if firstCheck in chkBuffer:
        checkVal = False
        while checkVal == False:
            count += 1
            #search the "buffer" for ensemble cluster
            if count < 10:
                #determine if the file has already been made
                curFileCheck = clustPre + sep + link + sep + dist + sep + '0' + str(count) + clustSuf
                if curFileCheck not in chkBuffer:
                    checkVal = True
                    clustFile = curFileCheck

            else:
                curFileCheck = clustPre + sep + link + sep + dist + sep + str(count) + clustSuf
                if curFileCheck not in chkBuffer:
                    checkVal = True
                    clustFile = curFileCheck

        plt.savefig(clustFile,dpi=600)
    else:
        clustFile = firstCheck  
        plt.savefig(clustFile,dpi=600)

    logging.info(': Success!')
    plt.show()

#convert seconds to HH:MM:SS
def timeConverter(runTime):
    '''
    Convert run time in seconds to a measure of Hours:Minutes:Seconds.

    Input:
    runTime - seconds of run time. 

    Output:

    runTime a string in HH:MM:SS format. 

    '''
    #immediately determine the number of full hours that were consumed
    hrs = runTime/3600
    hrs = int(hrs)

    #subtract the number of seconds that hrs took up and then determine the number of full minutes left in the
    #remaining seconds
    runTime -= (hrs*3600)
    mins = runTime/60
    mins = int(mins)
    runTime -= (mins*60)
    secs = int(runTime)
    runTime -= secs

    runTime = str(runTime)
    runTime = runTime.strip("0")
    runTime = runTime[0:3]

    if hrs < 10:
        hrsStr = '0' + str(hrs)
    else:
        hrsStr = str(hrs)

    if mins < 10:
        minsStr = "0" + str(mins)
    else:
        minsStr = str(mins)

    if secs < 10:
        secsStr = "0" + str(secs)
    else:
        secsStr = str(secs)

    runTime = hrsStr +':'+ minsStr +':'+ secsStr

    return runTime

def Validate(data,dists,num_groups):
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

    startPoint = 0.6*clusterings
    startPoint = int(startPoint)
    numIts = clusterings - startPoint

    numClusters = clusterings-startPoint
    numClusters = int(numClusters)
    val_index = np.zeros((2,clusterings))
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

        #calculate the average compactness of the clusters
        intraDist = sumIntra/(clusterings+1)
        
        # find the distance between the centers
        centerDists = pdist(centersNum)
        lenCenters = len(centerDists)
        centerDists = squareform(centerDists)
        centerDists = csr_matrix(centerDists)
        centerDistsInd = centerDists.nonzero()

        dataMST = np.zeros([lenCenters,1])

        for k in range(lenCenters):
            #input the values of each matrix element to the numpy array for saving to csv file.
            curVal0 = centerDistsInd[0][k]
            curVal1 = centerDistsInd[1][k]
            dataMST[k,0] = centerDists[curVal0,curVal1]
        

        #calculate the inter-cluster distance for the current cluster set-up
        if len(dataMST[:,0]) > 0:
            interDist = np.min(dataMST[:,0])
            #calculate the validation index
            val_index[0,i] = intraDist/interDist
            val_index[1,i] = clusterings - (i)

        elif len(dataMST[:,0]) == 0:
            interDist = 0
            #set the validation index to a large number since denominator would be zero in current config.
            val_index[0,i] = 1
            val_index[1,i] = clusterings - (i)

        if i == 0:
            logging.info(': 0% completed')

        elif i%10==0:
            #calculate the amount of time the algorithm has been running 
            curTime = time.time()
            runTime = (curTime - initTime)
            runTime = timeConverter(runTime)
            percent = 100*(1 - ((clusterings-i)/numIts))
            message1 = round(percent,2) 
            message1 = str(message1) + '% completed'
            logging.info(message1)
            message = str(i)+': ' + str(runTime)
            logging.info(message) 


        endTime = time.perf_counter()
        totalTime = endTime -startTime

        print(totalTime, curClustersLength) 

    runTime = time.time()-initTime
    runTime = timeConverter(runTime)
    logging.info(runTime)

    val_index = val_index[:,startPoint:clusterings]
    logging.info(': Cluster Validation complete!')
    return val_index

def who(curUser):
    '''
    Determines who is currently using the GUI for ease of saving of log files, and outputs. 

    Input:
    curUser - the current MSU NetID for current user. 

    Output:

    Name of current user. 
    '''
    #define dictionary of the NetID's for identification of the user for the naming of PDF
    file = open('C:/Users/Public/Documents/June_Lab_ClusteringGUI-master/Users.txt','r')
    contents = file.read()
    users = ast.literal_eval(contents)
    
    try:
        curUser = users[curUser]
    except:
        logging.error(': Unknown User')
        curUser = 'Anonyomous'
    return curUser 

def files(directory):
    #******************* Do we need this or should we update this? 
    '''
    Checks input directory for pngs.
    '''
    #checks the given path for files ending in .png
    directory = '*.png'
    chk_buffer = glob.glob(directory)

    return chk_buffer

def imageSize(file):
    '''
    Function checks for image size ratio for appropriate sizing in the pdf generator. 

    Input:
    file - either .png/.jpg/.jpeg (make sure to include full path)

    Output:
    ratio of the image size. 

    '''
    #determines the size of the image being imported to a pdf document
    im = Image.open(file)

    size = im.size
    ratio = size[0]/size[1]

    return ratio

def pdfHeader(file):
    '''
    Determine the appropriate header based upon the name of the file. 

    Input: 
    file - input full file path to the .png/.jpg/.jpeg image

    Output:

    string containing the appropriate header for the image. 
    '''
    #define a list with the shortened file names
    identifiers = ['pca_pair','pca_scree','pca_loading','pca_biplot','pls_pair','pls_score2d','pls_score3d','pls_loading','pls_cv','pls_imp','tree','fc','tt','volcano','PCA','pca_score2d','Peak_Intensity']

    #define a list with names for shortened file names
    names = ['PCA pairs','PCA Scree plots','PCA Loading','PCA biplot','PLS-DA Pairs','PLS-DA 2D score','PLS-DA 3D plot','PLS-DA Loading','PLS-DA CV','PLS-DA Imp','Dendrogram','Fold Change','T-test','Volcano Plot','PCA','PCA 2D score','Peak Intensities']

    for i in range(len(identifiers)):
        #find the length of the current evaluation list
        curList = identifiers[i]
        curLen = len(curList)

        #grab the first curLen characters from the file string
        curFile = file[0:curLen]

        #see if first part matches second part
        check = curList == curFile

        if check == True:
            header = names[i]
            return header
        if i == len(identifiers)-1 and check == False:
            header = 'Unknown test, have Brady add the new test to the text file.'
            return header

def recClusters(dataFinal,heatmapAxes,groupDendLeaves,metab_data):
    '''
    Determines the areas containing 100% clusters metabolites for the 13 separate clusterings performed. 

    Input:
    dataFinal - final formatted data (reorganized into the ensemble heatmap)
    heatmapAxes - input matplotlib axes for plotting of dashed lines.
    groupDendLeaves - leaves of current dendrogram being studied.
    metab_data - raw data.

    Output:
    sends data to the ensembleClustersOut function to export ensemble clusters. 
    '''

    #determine the appropriate number ensemble clusters and there location
    ensemMetabs = dataFinal.shape[0]

    #set the current row to check in the ensembles output
    j = 0
    while j < ensemMetabs:
        #look for the number of ones indicating the total number of 
        #metabolites in the current cluster.
        #
        #
        #
        #*******************************************************************
        # Need to make the number of clusterings available for dynamic analysis here. 
        #
        #
        #
        found = np.where(dataFinal[j,:]>12/13)

        #determine the length of the array found
        if len(found[0])==1:
            ensembleClustersOut(found[0],groupDendLeaves,metab_data)
            #then simply take the current j value and give it limits of j-0.5 to j+0.5, in the x and y directions.
            arrays = {0:np.linspace(j-0.5, j+0.5, num=5),1:np.linspace(j-0.5, j-0.5, num=5),2:np.linspace(j+0.5, j+0.5, num=5)}
            if len(found[0]) > 40:
                heatmapAxes.plot(arrays[0], arrays[1], 'r--', linewidth=3, markersize=3)
                heatmapAxes.plot(arrays[0], arrays[2], 'r--', linewidth=3, markersize=3)
                heatmapAxes.plot(arrays[1], arrays[0], 'r--', linewidth=3, markersize=3)
                heatmapAxes.plot(arrays[2], arrays[0], 'r--', linewidth=3, markersize=3)
                heatmapAxes.text(j,j,"1",color='red',fontsize=7)
            j += 1
            del(arrays)
        else:
            #use j and the maximum value from the connection search to locate create the outline
            maxMetab = max(found[0])
            ensembleClustersOut(found[0],groupDendLeaves,metab_data)
            arrays = {0:np.linspace(j-0.5, maxMetab+0.5, num=5),1:np.linspace(j-0.5, j-0.5, num=5),2:np.linspace(maxMetab+0.5, maxMetab+0.5, num=5)}
            if len(found[0]) > 40:

                heatmapAxes.plot(arrays[0], arrays[1], 'r--', linewidth=3, markersize=3)
                heatmapAxes.plot(arrays[0], arrays[2], 'r--', linewidth=3, markersize=3)
                heatmapAxes.plot(arrays[1], arrays[0], 'r--', linewidth=3, markersize=3)
                heatmapAxes.plot(arrays[2], arrays[0], 'r--', linewidth=3, markersize=3)
                numMetabs = str(maxMetab - j + 1)
                midPoint = (maxMetab + j)/2
                heatmapAxes.text(midPoint,midPoint,numMetabs,color='red',fontsize=7)
            j = maxMetab + 1
            del(arrays)
    return

def ensembleClustersOut(found,groupDendLeaves,metab_data):
    '''
    Intake the connected clusters of all ones and output a file with the appropriate title. This function should not be called on
    it's own but should be called through the recClusters function (which automatically recommends clusters for the user).


    Input:
    found - metabolites which have clustred together all 13 times. 
    groupDendLeaves - leaves of the current clustering of interest. 
    metab_data - raw data

    Output:
    csv files of the found metabolite clusters. 
    '''
    logging.info(': Creating ensemble clusters output files.')
    #take the found metabolites/biomarkers/etc. and grab the indicies which are matching in the Leaves of the dendrogram. 
    lenFound = len(found)
    foundMetabs = np.zeros((lenFound,metab_data.shape[1]-1))
    columnHeaders = foundMetabs.shape[1]

    #check the data type for the first column of the DataFrame
    columnsData = list(metab_data.columns)

    #create a dictionary of values for the indexing of the output of the ensemble fits.
    textList = metab_data[columnsData[0]]
    textList = textList.to_dict()

    #Get the raw data columns from the DataFrame
    columnsData = columnsData[1:len(columnsData)]#-1]
    rawData = metab_data[columnsData]
    rawData = rawData.to_numpy()

    #set prefix that will be used for all files
    ensemPre = 'EnsembleCluster'
    ensemSuf = '.csv'

    try:
        del(idents)
    except:
        logging.warning(': Deleting variable prior to its creation is not advised!!')

    if lenFound > 40:
        idents = []
        for i in range(lenFound):
            #find where in the row of metabData I need to extract from by determining which metabolite was clustered where from the groupDendLeaves
            curMetab = np.where(groupDendLeaves==found[i])
            
            foundMetabs[i,:] = rawData[curMetab[0][0],:]

            #save a list that contains identities for the study
            idents.append(textList[curMetab[0][0]])

        #create column headers for the data frame
        columns = []
        for i in range(columnHeaders-1):
            columns.append("M"+str(i+1))
        columns.append("rt_med")

        foundMetabs = pd.DataFrame(foundMetabs,columns=columns)

        #add identities to the first column of the data that will be output
        foundMetabs.insert(0, "Identities", idents, True)

        chkBuffer = glob.glob("*.csv")
        count = 1
        if 'EnsembleCluster01.csv' in chkBuffer:
            checkVal = False
            while checkVal == False:
                count += 1
                #search the "buffer" for ensemble cluster
                if count < 10:
                    #determine if the file has already been made
                    curFileCheck = ensemPre + '0' + str(count) + ensemSuf
                    if curFileCheck not in chkBuffer:
                        checkVal = True
                        ensemFile = curFileCheck
                else:
                    curFileCheck = ensemPre + str(count) + ensemSuf
                    if curFileCheck not in chkBuffer:
                        checkVal = True
                        ensemFile = curFileCheck
            foundMetabs.to_csv(ensemFile, index=False)
        else:
            ensemFile = ensemPre + '0'+ str(count) + ensemSuf 
            foundMetabs.to_csv(ensemFile, index=False)
        logging.info(':Success!')

def readInColumns(metab_data):
    '''
    Robust excel reading in tool. 

    Input:
    metab_data -raw excel file

    Output:
    Columns of metabolites, this goes with the convention that the submitted files contain the metabolite intensities in the 2->(N-1) columns.
    '''
    #creating a numpy array that is the size of the data that is being read in.
    data = np.zeros((metab_data.shape[0],metab_data.shape[1]-2))

    columnsData = list(metab_data.columns)

    for i in range(len(columnsData)):
        if i > 0 and i < len(columnsData)-1:
            #try to add the values to the current medians values
            try:
                medianCur = metab_data[columnsData[i]]
            except:
                logging.error(': Unable to read in column headers. ')

            #add the medians data to the array to be clustered
            data[:,i-1] = medianCur

    return data