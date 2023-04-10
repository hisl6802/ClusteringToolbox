import numpy as np
import itertools
import statistics as stat
from scipy.cluster.hierarchy import dendrogram
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist,squareform
from scipy.sparse import csr_matrix
from matplotlib import pyplot as plt
from matplotlib import colors as mcolors
import pandas as pd
import seaborn as sns
import tkinter as tk
from tkinter import HORIZONTAL, filedialog, messagebox
from tkinter.ttk import Button
from tkinter.ttk import Progressbar
import logging, time, glob,sys,os,ast
from PIL import Image
from seaborn.matrix import clustermap
import math
import config

from fpdf import FPDF

from Bio.KEGG import REST
from Bio.KEGG import Compound
from Bio.KEGG import Enzyme
from tqdm import tqdm


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
        messagebox.showerror(title ='Error' ,message='Please select file to continue!')
        return

    #Open the excel file that the user to looking to use for their clustering analysis
    try:
        data = pd.read_excel(file)
        del(file)
    except:
        logging.error(': Failed to read in the excel sheet, it is recommended that an excel workbook with a single sheet is uploaded!')
        messagebox.showerror(title='Error',message='Failed to read in the excel sheet, it is recommended that an excel workbook with a single sheet is uploaded!')
        return
    logging.info(': User file opened and submitted to the function.')


    #send the data to the data checker

    # data = dataCheck(data)
    return data



###------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
###-------------------------------------------------------------- DATA SCALING---------------------------------------------------------------------------------------------------------
###------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


#Standardizing the data that is input to python. (Auto-scaling in Metaboanalyst)
def standardize(data):
    '''
    Standardize the input data for best clustering results. This should be the same as Auto-scaling
    
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
def normStandardize(data,first):
    '''
    Normalize input data that has been standardized for normalized clustergrams.

    Input:

    data - all data. 
    leaves - leaves of the groups you would like to normalize. 

    ***Currently the only possible normalization is the first column of data. Will be updated soon. 

    Output:
    
    Normlized standardardized data. 

    '''

    first = int(first)-1
    first = int(first)

    logging.info(': Normalizing strandardized data.')
    #The mean for a normalized data set should always be 1 for all of the metabolites.
    mean = data[first]

    #initialize the needed arrays
    dataCur = np.zeros(data.shape[0])
    #Calculate the residuals of the data
    residuals = 0
    for j in range(data.shape[0]):
        residuals += (data[j]-mean)**2

    #calculate the standard deviation
    std_data = (residuals/(data.shape[0]-1))**0.5 
    #add jitter to the standard deviation of the data to ensure it does not get inf in result
    std_data += 0.001
    #standardize each row of data
    for i in range(data.shape[0]):
        #Input the data into an array of standardized (auto-scaling in metaboanalyst) values
        dataCur[i] = (data[i]-mean)/std_data

    return dataCur


def meanCentering(data):
    '''
    Mean centering the data, subtracting the mean from all data values

    Input:

    data - raw data

    Output:

    mean centered data
    '''

    logging.info(': Mean-centering the data.')
    mean_data = stat.mean(data)

    dataCur = np.zeros(data.shape[0])

    for i in range(data.shape[0]):
        dataCur[i] = (data[i]-mean_data)
    
    return dataCur

def paretoScaling(data):
    '''
    Pareto scaling is similar to Auto-scaling but with data divided by the square root of the standard deviation.

    Input:

    data - raw data

    Output:

    Pareto scaled data
    '''

    logging.info(': Pareto scaling the data.')
    mean_data = stat.mean(data)
    std_data = stat.stdev(data)

    dataCur = np.zeros(data.shape[0])

    if std_data == 0:
        for i in range(data.shape[0]):
            dataCur[i] = data[i]-mean_data
    else:
        for i in range(data.shape[0]):
            dataCur[i] = (data[i]-mean_data)/(std_data**.5)

    return dataCur

def rangeScaling(data):
    '''
    Range scaling is similar to Auto-scaling but with the data divided by the range of the data.

    Input: 

    data - raw data

    Output:

    Range scaled data
    '''
    logging.info(': Range scaling the data.')
    mean_data = stat.mean(data)
    max = np.max(data)
    min = np.min(data)
    rangeData = max-min

    dataCur = np.zeros(data.shape[0])
    if rangeData == 0:
        for i in range(data.shape[0]):
            dataCur[i] = data[i]-mean_data
    else:
        for i in range(data.shape[0]):
            dataCur[i] = (data[i]-mean_data)/rangeData

    return dataCur



###------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
###------------------------------------------------------------------------- DATA TRANSFORMATIONS -------------------------------------------------------------------------------------
###------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

### Log Transformation
def logTrans(data):
    '''
    Log transform

    Input:

    Required 
    data - raw data

    Output:

    Outputs the input data log transformed
    '''
    #given a row of data from the data take a log transform 
    #add small amount amount to get past zeros
    data +=0.1
    dataUpdate = np.log10(data) 

    
    return dataUpdate

### Square Root Transformation
def sqrtTrans(data):
    '''
    Square root transformation

    Input:

    Require
    data - raw data

    Output:

    Outputs the input data square root transformed
    '''

    dataUpdate = np.sqrt(data)

    return dataUpdate

### Cube Root Transformation
def cubeRtTrans(data):
    '''
    Cubic root transformation

    Input:

    Required
    data - raw data

    Output:

    Outputs the input data cube root transformed
    '''

    dataUpdate = np.cbrt(data)

    return dataUpdate


#ensemble dendrogram function
def createEnsemDendrogram(data,metab_data,norm=1,minMetabs=0, numClusts = 13, link='ward',dist='euclidean',func="ensemble",colMap ='viridis'):
    '''
    Create ensemble dendrogram
    
    Input:

    Required:

    data - co-occurence matrix
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



    g = sns.clustermap(data, figsize=(7, 5), yticklabels=False, xticklabels=False, row_linkage=linkageMetabOut, col_linkage=linkageGroupOut, cmap=colMap, cbar_pos=(0.01, 0.8, 0.025, 0.175))
    ax = g.ax_heatmap
    axD = g.ax_row_dendrogram

    
    recClusters(dataFinal,ax,groupDendLeaves,metab_data,minMetabs,numClusts)
    plt.savefig('EnsembleClustergram01.png',dpi=600,transparent=True)
    plt.show()

#dendrogram function
def create_dendrogram(data, col_groups, norm=1,colOrder =[], link='ward',dist='euclidean', color='viridis'):
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

    #set out color options and map to the groups. 
    colorOpts = ('b','y','m','r','k','#929292')
    
    #find the unique groups
    col_groupsUni = col_groups.unique()

    #create a dictionary for mapping color options
    colRefDict = {}
    for i in range(len(col_groupsUni)):
        colRefDict[col_groupsUni[i]] = colorOpts[i]

    colSeries = col_groups.map(colRefDict)
    col_groups = colSeries.to_list()
    

    #Create the linkage matrix
    linkageMetabOut = linkage(data,link,dist)

    #tranpose the data and then run an analysis on the groups.
    try:
        groupCluster = np.transpose(data)
        logging.info(': Clustergram being generated.')
    except:
        logging.error(': Inappropriate entry for the create_dendrogram function. Check docs.')
        messagebox.showerror(title ='Error' ,message='Inappropriate entry for the create_dendrogram function. Check docs.')

    #Create a linkage matrix for the data
    linkageGroupOut = linkage(groupCluster,link,dist)

    
    if norm == 0:
        #g = sns.clustermap(data, figsize=(7, 5), row_linkage=linkageMetabOut, col_linkage=linkageGroupOut, cmap=color, cbar_pos=(0.01, 0.8, 0.025, 0.175))
        g = sns.clustermap(data, method=link,metric=dist, figsize=(7, 5), col_cluster=True,col_colors=col_groups,cmap=color,yticklabels=False,xticklabels=True)
        plt.savefig('Clustergram.png',dpi=600,transparent=True)
        plt.show()

    elif norm == 1:
        #reordering data for users based upon preference
        colOrder = [int(i) for i in colOrder]
        colOrder = [i-1 for i in colOrder]
        col_groups = [col_groups[i] for i in colOrder]

        data[:,:] = data[:,colOrder]
        g = sns.clustermap(data,method=link,metric=dist,figsize=(7,5), col_cluster=False,cmap=color,col_colors=col_groups,yticklabels=False,xticklabels=False)
        plt.savefig('Clustergram.png',dpi=600,transparent=True)
        plt.show()

    elif norm == 2:
        #reordering the data for user based upon preference
        colOrder = [int(i) for i in colOrder]
        colOrder = [i-1 for i in colOrder]
        col_groups = [col_groups[i] for i in colOrder]

        data[:,:] = data[:,colOrder]
        g = sns.clustermap(data,method=link,metric=dist,figsize=(7,5), col_cluster=False,cmap=color,col_colors=col_groups,yticklabels=False,xticklabels=False)
        plt.savefig('Clustergram.png',dpi=600,transparent=True)
        plt.show()

        

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

            curCon1Connect = 0
            curCon2Connect = 0
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
                        logging.warning(': Issue clustering the data a duplication has been discovered.')

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
                        logging.warning(': Issue clustering the data a duplication has been discovered.')


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
    
    plots the clustergram. 

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

        plt.savefig(clustFile,dpi=600,transparent=True)
    else:
        clustFile = firstCheck  
        plt.savefig(clustFile,dpi=600,transparent=True)

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
        intraDist = sumIntra/(dists.shape[0])
        
        # find the distance between the centers
        interDist = 1000
        centerDists = pdist(centersNum)
        for j in range(len(centerDists)):
            #determine the minimum non-zero value in the distance array
            #index does not matter for this algorithm. 
            if centerDists[j] < interDist and centerDists[j] != 0:
                interDist = centerDists[j]

        #calculate the inter-cluster distance for the current cluster set-up
        if len(centerDists) > 0:
            #calculate the validation index
            val_index[0,i] = intraDist/interDist
            val_index[1,i] = clusterings - (i)
            val_index[1,i] = len(data[0])

        # elif len(dataMST[:,0]) == 0:
        elif len(centerDists) == 0:
            interDist = 0
            #set the validation index to a large number since denominator would be zero in current config.
            val_index[0,i] = 1
            val_index[1,i] = clusterings - (i)
            val_index[1,i] = len(data[0])

    val_index = val_index[:,startPoint:clusterings]

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
        messagebox.showinfo(title='New User?', message = 'Are you a new user on the server? Contact Brady to have your name added!')
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

def recClusters(dataFinal,heatmapAxes,groupDendLeaves,metab_data,minMetabs,numClusts):
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
        found = np.where(dataFinal[j,:]>(numClusts-1)/numClusts)

        #determine the length of the array found
        if len(found[0])==1:
            maxMetab = found[0]
            ensembleClustersOut(found[0],groupDendLeaves,metab_data,minMetabs)
            #then simply take the current j value and give it limits of j-0.5 to j+0.5, in the x and y directions.
            arrays = {0:np.linspace(j, maxMetab+1, num=5),1:np.linspace(j, j, num=5),2:np.linspace(maxMetab+1, maxMetab+1, num=5)}

            if len(found[0]) >= minMetabs:
                heatmapAxes.plot(arrays[0], arrays[1], 'r--', linewidth=3, markersize=3)
                heatmapAxes.plot(arrays[0], arrays[2], 'r--', linewidth=3, markersize=3)
                heatmapAxes.plot(arrays[1], arrays[0], 'r--', linewidth=3, markersize=3)
                heatmapAxes.plot(arrays[2], arrays[0], 'r--', linewidth=3, markersize=3)
                heatmapAxes.text(j+0.5,j+0.5,"1",color='red',fontsize=7)
            j += 1
            del(arrays)

        else:
            #use j and the maximum value from the connection search to locate create the outline
            maxMetab = max(found[0])
            ensembleClustersOut(found[0],groupDendLeaves,metab_data,minMetabs)
            arrays = {0:np.linspace(j, maxMetab+1, num=5),1:np.linspace(j, j, num=5),2:np.linspace(maxMetab+1, maxMetab+1, num=5)}

            if len(found[0]) >= minMetabs:

                heatmapAxes.plot(arrays[0], arrays[1], 'r--', linewidth=3, markersize=3)
                heatmapAxes.plot(arrays[0], arrays[2], 'r--', linewidth=3, markersize=3)
                heatmapAxes.plot(arrays[1], arrays[0], 'r--', linewidth=3, markersize=3)
                heatmapAxes.plot(arrays[2], arrays[0], 'r--', linewidth=3, markersize=3)
                numMetabs = str(maxMetab - j + 1)
                midPoint = ((maxMetab+1)+j)/2
                heatmapAxes.text(midPoint,midPoint,numMetabs,color='red',fontsize=7)
            j = maxMetab + 1
            del(arrays)
    messagebox.showinfo(title="Success",message="Successfully created ensemble clustergram and cluster files!!")
    return

def ensembleClustersOut(found,groupDendLeaves,metab_data,minMetabs):
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
    ensemSuf = '.xlsx'

    try:
        del(idents)
    except:
        logging.warning(': Deleting variable prior to its creation is not advised!!')


    #check for an EnsembleOutputFiles directory in the current directory
    if os.path.isdir('EnsembleOutputFiles'):
        os.chdir('EnsembleOutputFiles')
    else:
        os.mkdir('EnsembleOutputFiles')
        os.chdir('EnsembleOutputFiles')

    if lenFound >= minMetabs:
        idents = []
        for i in range(lenFound):
            #match the 
            foundMetabs[i,:] = rawData[groupDendLeaves[found[i]],:]

            #save a list that contains identities for the study
            idents.append(textList[groupDendLeaves[found[i]]])

        #create column headers for the data frame
        columns = []
        for i in range(columnHeaders-1):
            columns.append("M"+str(i+1))
        columns.append("rt_med")

        foundMetabs = pd.DataFrame(foundMetabs,columns=columns)

        #add identities to the first column of the data that will be output
        foundMetabs.insert(0, "Identities", idents, True)

        chkBuffer = glob.glob("*.xlsx")
        count = 1
        if 'EnsembleCluster01.xlsx' in chkBuffer:
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
            foundMetabs.to_excel(ensemFile, index=False)
        else:
            ensemFile = ensemPre + '0'+ str(count) + ensemSuf 
            foundMetabs.to_excel(ensemFile, index=False)
        

    #change directory back to the original directory
    os.chdir('..')
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
                messagebox.showerror(title='Error',message='Program unable to read column headers, check submitted file!')

            #add the medians data to the array to be clustered
            data[:,i-1] = medianCur
    
    return data
    
def select(index,dend,link,linkDir,linkClust,data_orig):
    '''
    Function responsible for the coloring clustergram based upon users selection. This function is also responsible for saving the selected cluster
    

    Input:
    Index of dendrogram
    dendrogram
    linkage
    linkage directory 
    linkage clustering
    original data

    Output:
    excel workbook with selected values

    '''
    
    #getting the current color number
    colSel = config.colorNum

    if colSel == 0:
        colSel = 1
        config.colorNum = 1

    elif colSel == 1:
        colSel = 2
        config.colorNum = 2

    else:
        colSel = 0
        config.colorNum = 0


    #grab the first value of the list. 
    dCord = -index[0]
    iCord = -index[1]
    #value which needs to be subtracted from each of the indicies
    sub = linkDir[0][0]
    try:
        countLink = 0
        curSelection = -1
        #find the list index which contains this value found.
        for i in range(len(link)):
            if abs(link[i][2] - dCord) < 0.0000001:
                curSelection = i

                countLink += 1

        clustMetabs = linkClust[curSelection][0]

        colors = ['steelblue','darkkhaki','darkorchid']
        lenColors = len(colors)
        if curSelection != -1:
            #if current selection is valid grab the curSelection list from the linkage directory
            curList = linkDir[curSelection]

            #if it's length is greater than 1 loop over the selections
            if len(curList) > 1:
                count = []
                for i in range(len(curList)):
                    curCord = curList[i]-sub
                    curDist = link[curCord][2]
                    for j in range(len(dend['dcoord'])):

                        if abs(dend['dcoord'][j][1] - curDist) < 0.0000001:
                            curCord = j
                            count.append(j)

                    #determine the number of counts in the current list
                    if len(count) > 0:
                        for j in count:
                            # if abs(iCord - dend['icoord'][j][1]) < 15.001 and abs(dend['icoord'][j][2]- iCord) < 15.001:
                            curCord = j
                    #plot the appropriate linkages        
                    x = np.array(dend['icoord'][curCord])
                    y = np.array(dend['dcoord'][curCord])
                    plt.plot(-y,-x,colors[colSel])
                    plt.draw()
            else:
                count = []
                curCord = curList[0]-sub
                curDist = link[curCord][2]
                
                for j in range(len(dend['dcoord'])):
                        if abs(dend['dcoord'][j][1] - curDist) < 0.0000001:
                            curCord = j
                            count.append(j)
                            
                #determine the number of counts in the current list
                if len(count) > 0:
                    for j in count:
                        # if abs(iCord - dend['icoord'][j][1]) < 15.001 and abs(dend['icoord'][j][2]- iCord) < 15.001:
                        curCord = j
                
                #plot the appropriate linkages
                x = np.array(dend['icoord'][curCord])
                y = np.array(dend['dcoord'][curCord])
                #plt.plot(-y,-x,colors[bSel])
                plt.plot(-y,-x, colors[colSel])
                plt.draw()
    except:
        logging.error('Cannot select the side of a linkage!')
        messagebox.showerror(title="Error",message='Make sure to select vertical lines only, selecting horizontal lines will continue to result in this error! YOU MAY NEED TO RESTART THE CLUSTER SELECTION!')
    
    
    with open('ClusterReference.txt') as f:
        lines = f.read()



    
    #make the current selection into a .csv file to be submitted to the peaks to Pathways function.
    selectedMetabs = np.zeros([len(clustMetabs),data_orig.shape[1]])
    p2pMetabs = np.ones([data_orig.shape[0],3])

    p2pMetabs[:,0] = data_orig[:,0]
    p2pMetabs[:,2] = data_orig[:,data_orig.shape[1]-1]

    for j in clustMetabs:
        p2pMetabs[j,1] = 0.04
    
    #look for matches in the list
    dendrogramLeaveLoc = []
    dendLeaves = dend['leaves']
    for i in range(len(clustMetabs)):
        #for each clustMetabs in the list put the value into a numpy array
        selectedMetabs[i,:] = data_orig[clustMetabs[i],:]
        dendrogramLeaveLoc.append(dendLeaves.index(clustMetabs[i]))


    lines = lines.split("\n")
    lenTxt = len(lines)
    newCluster = "\n" 
    open('ClusterReference.txt','a').write(newCluster)
    open('ClusterReference.txt','a').write(str(len(dendrogramLeaveLoc)))

    clustPre = 'Cluster'
    clustSuf = '.xlsx'
    #create column headers for the data frame
    columnHeaders = selectedMetabs.shape[1]
    columns = []
    for i in range(columnHeaders-1):
        if i == 0:
            columns.append("Identities")
        else:
            columns.append("M"+str(i))
    columns.append("rt_med")

    foundMetabs = pd.DataFrame(selectedMetabs,columns=columns)
    p2pFile = pd.DataFrame(p2pMetabs,columns=["m.z","p.value","r.t"])
    p2pFile = p2pFile.sort_values(by=['p.value'], ascending=True)

    chkBuffer = glob.glob("*.xlsx")
    count = 1
    if 'Cluster01.xlsx' in chkBuffer:
        checkVal = False
        while checkVal == False:
            count += 1
            #search the "buffer" for ensemble cluster
            if count < 10:
                #determine if the file has already been made
                curFileCheck = clustPre + '0' + str(count) + clustSuf
                if curFileCheck not in chkBuffer:
                    checkVal = True
                    clustFile = curFileCheck
            else:
                curFileCheck = clustPre + str(count) + clustSuf
                if curFileCheck not in chkBuffer:
                    checkVal = True
                    clustFile = curFileCheck
        foundMetabs.to_excel(clustFile, index=False)
        p2pF = "P2P_"+ clustFile.rstrip('.xlsx') + '.csv'
        p2pFile.to_csv(p2pF,index=False)
    else:
        clustFile = clustPre + '0'+ str(count) + clustSuf 
        foundMetabs.to_excel(clustFile, index=False)
        p2pF = "P2P_"+ clustFile.rstrip('.xlsx') +'.csv'
        p2pFile.to_csv(p2pF,index=False)
    logging.info(':Success!')


def linkDir(linkageOne,maxIndex):
    '''
    Creates dictionary containing all of the indicies and corresponding parameter which goes with the index.

    Input:

    linkage function output
    maxIndex??

    Output:

    Dictionary containing all of the iterations of the creation of the clustering solution

    '''
    #initializing dictionary for storage of the linkage names.
    linkageDir = {}
    #create list with the linkage names embedded. 
    linkNums = []

    #fill the array with the linkage identifiers
    for i in range(linkageOne.shape[0]):
        curList = []
        curList.append(int(linkageOne[i][0]))
        curList.append(int(linkageOne[i][1]))
        linkNums.append(curList)

    for i in range(len(linkageOne)):
        if i == 0:
            #with the initial linkage just put the first index in.
            linkageDir[i] = [maxIndex+(i+1)]
        else:
            curList = [maxIndex+(i+1)]
            #check to see if items in curList are greater than maxIndex.
            checks = []
            if linkNums[i][0] > maxIndex:
                curList.append(linkNums[i][0])
                checks.append(linkNums[i][0])
            if linkNums[i][1] > maxIndex:
                curList.append(linkNums[i][1])
                checks.append(linkNums[i][1])
            
            #loop over checks looking for connecting values in the list that needs to be checked.
            if len(checks) > 0:
                for k in range(len(checks)):
                    #checks length, and then add all other values in this list to the current list
                    curCheck = linkageDir[checks[k]-(maxIndex+1)]
                    if len(curCheck) > 1:
                        for j in range(1,len(curCheck)):
                            #append the linkage numbers to the list of interest curList
                            curList.append(curCheck[j])

            linkageDir[i] = curList

    return linkageDir


def readAndPreProcess(file='',transform = 'None', scale ='None',func='else',first ='1'):
    '''
    readAndPreProcess

    This function is designed to remove the reading in and pre-processing of the data from the beginning of each function which needs to read-in and pre-process the data, saying large amounts of lines in this program.

    Input:
    transform
    scale
    func

    Output:
    Pre-processed data

    '''

    #check that the file the user selects is appropriate
    ###Should only be used when reading in excel files.
    metab_data = fileCheck(file=file)

    if metab_data is None:
        #log error message and return for soft exit.
        logging.error(': Error loading in the Excel sheet.')
        return

   
    if func =='CC':
        metab_dataCol = list(metab_data.columns)
        col_groups = metab_data.iloc[0]
        col_groups = col_groups.to_list()
        col_groups.pop(0)
        col_groups.pop(len(col_groups)-1)
        col_groups = pd.Series(col_groups,copy=False)
        metab_data = metab_data.drop(0,axis=0)
        #read in data
        # data = readInColumns(metab_data)
        data = readInColumns(metab_data)
        #put the data through the appropriate transformations. 
        data = transformations(data,transform=transform,scale=scale, first=first)

        return data, col_groups

    elif func == 'ANHM':
        #remove the first row of data.
        metab_data = metab_data.iloc[1:,:]
        # data = readInColumns(metab_data)
        data = readInColumns(metab_data)

        data = transformations(data,transform=transform,scale=scale,first=first)
        return data

    elif func =='Raw':
        metab_data = metab_data.drop(0,axis=0)
        metab_data.reset_index(inplace=True,drop=True)
        return metab_data

    else:
        #read in data
        # data = readInColumns(metab_data)
        data = readInColumns(metab_data)
        #transform the data
        data = transformations(data,transform=transform,scale=scale,first=first)
        return data


    
def createHeatmapFig(clMap):
    '''
    This function is responsible for creating the heatmap that is submitted by the user. 

    Input: 
    clMap: choosen color map scheme

    Output:
    editable pdf that with each of the selected clusters.
    '''

    #read in excel sheet of Heatmap.xlsx
    fileName = filedialog.askopenfilename()
    try:
        data = pd.read_excel(fileName)
    except:
        logging.error(': Likely that no file was selected, or there was an issue connecting to the drive!')
        messagebox.showerror(title="Error",message="Could not find file or no file was selected!")
        return
            
    #input the dataframe into the heatmap function
    g = sns.heatmap(data,yticklabels=False,cmap=clMap)

    #save the heatmap of the data
    plt.savefig('Heatmap.png',dpi=600)
    del(g)

    #create the pdf for publication
    pdf = FPDF('L','pt',(2880,3840))
    pdf.add_page()
    pdf.set_line_width(10)
    pdf.set_font('Arial','B',54)

    #add image to pdf
    pdf.image('Heatmap.png',0,0)


    #open the ClusterReference.txt file.
    with open('ClusterReference.txt') as f:
        lines = f.read()
    lineNew = lines.split("\n")
    
    numMetabs = int(lineNew[1])

    #determine the pixel ratio
    pR = 2220/numMetabs

    boxHeights = []
    for i in range(len(lineNew)):
        if i >= 2:
            #add cluster length in pixels to box height list
            boxHeights.append(pR*int(lineNew[i]))


    if len(boxHeights)%2 == 0:
        #alternate ceiling and floor for length of the list.
        for i in range(len(boxHeights)):
            if i%2 ==0:
                boxHeights[i] = math.ceil(boxHeights[i])
            else:
                boxHeights[i] = math.floor(boxHeights[i])

    else:
        #alternate ceiling and floor for length of the list until the last entry in the list.
        for i in range(len(boxHeights)):
            if i%2 == 0:
                #found using the ceiling no matter pixel fraction
                boxHeights[i] = math.ceil(boxHeights[i])

            else:
                #round using the floor no matter pixel fraction
                boxHeights[i] = math.floor(boxHeights[i])


    startPointX = 177
    startPointY = 346
    boxWidth = 2687

    curX = pdf.get_x()

    curY = pdf.get_y()
    
    diffX = startPointX - curX
    diffY = startPointY - curY

    pdf.ln(317.65)
    pdf.cell(149.65,h=100,ln=0)

    firstMetab = 1
    for i in range(len(boxHeights)):

        if i == 0:
            #make a box with the appropriate height and width, and the correct starting position.
            pdf.rect(startPointX,startPointY,boxWidth,boxHeights[i])
            f1 = str(i+1)
            pdf.cell(300,h=math.floor(boxHeights[i]/2),txt=f1,ln=2,align='C')
            eR = int((firstMetab-1) + int(lineNew[2]))
            rMetab = str(firstMetab) + ':' + str(eR)
            pdf.cell(300,h=math.floor(boxHeights[i]/2),txt=rMetab,ln=2,align='C')


        else:
            startPointY += boxHeights[i-1]
            pdf.rect(startPointX,startPointY,boxWidth,boxHeights[i])
            #put the i + 1 value in as the text
            f1 =str(i + 1)
            pdf.cell(300,h=math.floor(boxHeights[i]/2),txt =f1,ln=2,align='C')
            sR = eR + 1
            eR += int(lineNew[i+2])
            rMetab = str(sR) + ':' + str(eR)
            pdf.cell(300,h=math.floor(boxHeights[i]/2),txt=rMetab, ln=2,align='C')

    
    pdf.output('SelectedClusters.pdf','F')
    messagebox.showinfo(title='Success', message='Success, the pdf has been created!')
    return


def valPlotting(valIndex, mstOut, valMet = 'KMeansBased'):

    '''
    This function is responsible for plotting the validation outcome generated.

    Input:
    validation Index: containing the number of clusters and the validation metric values
    mstOut: containing the generated MST
    valMet: the validation metric run

    Output:
    Plot of the validation output
    csv of the validation output
    csv of the MST
    '''

    #put the number of clusters in K and the validation metric in y
    if valMet == 'KMeansBased':
        K = valIndex[:,1]
        y = valIndex[:,0]
        y[y.shape[0]-1] = y[y.shape[0]-2]*2

        y = 1/y

        for i in range(y.shape[0]):
            if K[i] >1:
                y[i] = (y[i]*K[i])/(K[i]-1)

    else:
        K=valIndex[:,0,1]
        y=valIndex[:,0,0]
        if valMet =='DBI':
            y[y.shape[0]-1] = y[y.shape[0]-2]*2
            y = 1/y

        else:
            y[y.shape[0]-1] = y[y.shape[0]-2]/2
            
        for i in range(y.shape[0]):
            if K[i] > 1:
                y[i] = (y[i]*K[i])/(K[i]-1)

    minValIndex = np.max(y)

    indMin = np.where(y==minValIndex)

    valOut = np.zeros((2,valIndex.shape[0]))
    for j in range(valIndex.shape[0]):
        #flip the array so that lower clusters start at lower positions
        if valMet == 'KMeansBased': 
            valOut[0,j] = valIndex[valIndex.shape[0]-j-1,0]
            valOut[1,j] = valIndex[valIndex.shape[0]-j-1,1]
        else:
            valOut[0,j] = valIndex[valIndex.shape[0]-j-1,0,0]
            valOut[1,j] = valIndex[valIndex.shape[0]-j-1,0,1]

    valIHeaders = list(valOut[1,:])
    valIndex = pd.DataFrame(valOut,columns=valIHeaders)
    valIndex = valIndex.drop(1,axis=0)
    rowLabels = ["Validation Index"]
    valIndex.insert(0,"Clusters",rowLabels)
    
    #save to a csv file
    mstOutFileName = 'MST_branches_' + valMet +'.csv'
    mstOut.to_csv(mstOutFileName, index=False)

    #save validation measure to csv file
    valIndexFileName = 'valIndex_' + valMet + '.csv'
    valIndex.to_csv(valIndexFileName, index=False)

    #logging the completion of the Minimum spanning tree
    logging.info(': Sucessfully completed clustering validation!')
    ax = plt.subplot(111)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_linewidth(4)
    ax.spines['bottom'].set_linewidth(4)
    ax.tick_params(width=4,length=10)

    ax.plot(K,y,linewidth=2.5)
    ax.plot(K,y,'k.')
    ax.plot(K[indMin[0]],minValIndex,'r.',markersize=20)
    font = {'family': 'serif','color':  'black','weight': 'bold','size': 20}
    ax.annotate(str(int(K[indMin[0]]))+' - Clusters!!',(K[indMin[0]]+0.1,minValIndex),xycoords='data',
    xytext=(K[indMin[0]]+0.2,minValIndex), textcoords='data',
                        horizontalalignment='left', verticalalignment='bottom',fontsize=18,fontname="Arial")
    plt.xlabel('Clusters',fontsize=28,fontname="Arial")
    plt.ylabel('Validation Index',fontsize=28,fontname="Arial")
    plt.xticks(fontsize=24,fontname="Arial")
    plt.yticks(fontsize=24,fontname="Arial")
    if valMet == "KMeansBased":
        title = "K-Means Based Validation"
    else:
        title = valMet + " Validation"
    plt.title(title,pad = 15,fontsize=36,fontname="Arial")
    pltFileName = valMet+"Validation.png"
    plt.savefig(pltFileName,bbox_inches='tight',dpi=600,transparent=True)
    plt.show()



def dataCheck(data):
    '''
    This function is responsible for checking for matching values within the input data, and reducing matching values down to a single value

    Input:
    original raw data

    Output:
    Corrected data.
    '''   
 
    #read in the data from fileCheck and look for matching values
    # data = readInColumns(data)

    dataMatches = {}
    toDelete = []
    dists = pdist(data)
    zerosL = dists==0

    indicies =list(itertools.combinations(range(data.shape[0]),2))
    #zeros checking
    zeros = ()

    zerosL = np.where(zerosL==True)[0].tolist()
    print('Starting...')
    # create reporting interval 
    report = np.linspace(0, len(zerosL),num=200,dtype=int)
    for i in range(len(zerosL)):
        if len(np.where(report==i)[0])>0:
            print((i/len(zerosL))*100)
        zeros += indicies[zerosL[i]]

    # create a numpy array of zeros
    zeros = np.array(zeros)
    # get the unique values of the array
    zeros = np.unique(zeros)

    #indicies of function, remove indicies if they do not match up with the unique values of the zeros 
    ind = list(range(data.shape[0]))
    
    for i in range(zeros.shape[0]):
        ind.remove(zeros[i])

    dataZeros = np.delete(data,ind,axis=0)

    for i in range(dataZeros.shape[0]):
        curMatch = []
        #skip the iteration if the value is already flagged to be deleted due to matching.
        if i in toDelete:
            continue

        for j in range(i+1,dataZeros.shape[0]):
            #skip the iteration if the value is already flagged to be deleted due to matching. 
            if j in toDelete:
                continue

            #check the current array v. the other arrays
            curCheck = dataZeros[i,:] == dataZeros[j,:]

            if np.all(curCheck):
                if len(curMatch) == 0:
                    #inputting the matched indicies to a list to be input to a dictionary.
                    curMatch.append(i)
                    curMatch.append(j)
                    toDelete.append(j)
                else:
                    curMatch.append(j)
                    toDelete.append(j)
        
        #if the length of list curMatch is non-zero add to the dictionary
        if len(curMatch)>0:
            #get dictionary of matching intensities and input to the dictionary using current length as key
            dictLen = len(dataMatches)

            #if this gives trouble automatically update the dictionary with each found match for each i
            dataMatches[dictLen] = curMatch
    
    #map to the appropriate indicies of the zeros then delete these indicies
    for i in range(len(toDelete)):
        toDelete[i] = zeros[toDelete[i]]


    print(toDelete)
    #delete the extraneous matching rows.
    data = np.delete(data,toDelete,axis=0)
    
    dataMessage = str(len(toDelete))
    removedIndicies = pd.DataFrame(toDelete)
    removedIndicies.to_excel('RemovedIndicies.xlsx',index=False)
    messagebox.showwarning(title="Matching Values removed",message=dataMessage + " matching values found and removed")

    return data, toDelete

        
def transformations(data, transform='None', scale='None',first='1'):
    '''
    This function is responsible for taking inputs from the broad range of functions needed data transformed or scaled and updating the data

    Input:
    data: raw data
    transform: selected transformation
    scale: selected data scaling

    Output:
    Pre-processed data

    '''

    ###-------------------------------------------------------------------------------------------------------------------------------------------------------------
    ###------------------------------------------------------------- Transforming and Scaling Data -----------------------------------------------------------------
    ###-------------------------------------------------------------------------------------------------------------------------------------------------------------

    ### UPDATED THE RANGE FOR THE FOR LOOP TO DATA INSTEAD OF METAB_DATA

    #Log transform no scaling
    if transform == 'Log transformation' and scale == 'None':
        for i in range(data.shape[0]):
            data[i,:] = logTrans(data[i,:])

    #Log transform, mean centering
    elif transform == 'Log transformation' and scale == 'Mean centering':
        #log transform
        for i in range(data.shape[0]):
            data[i,:] = logTrans(data[i,:])

        #mean center
        for i in range(data.shape[0]):
            data[i,:] = meanCentering(data[i,:])
    
    #log transform, auto scaling
    elif transform == 'Log transformation' and scale == 'Auto Scaling':
        #log transform
        for i in range(data.shape[0]):
            data[i,:] = logTrans(data[i,:])

        #Auto scale
        for i in range(data.shape[0]):
            data[i,:] = standardize(data[i,:])

    #log transform, pareto scaling
    elif transform == 'Log transformation' and scale == 'Pareto Scaling':
        #log transform
        for i in range(data.shape[0]):
            data[i,:] = logTrans(data[i,:])

        #Pareto scale
        for i in range(data.shape[0]):
            data[i,:] = paretoScaling(data[i,:])

    #log transform, range scaling
    elif transform == 'Log transformation' and scale == 'Range Scaling':
        #log transform
        for i in range(data.shape[0]):
            data[i,:] = logTrans(data[i,:])

        #Range scale
        for i in range(data.shape[0]):
            data[i,:] = rangeScaling(data[i,:])

    elif transform == 'Log transformation' and scale =='NormStand':
        #log transform
        for i in range(data.shape[0]):
            data[i,:] = logTrans(data[i,:])
        
        for i in range(data.shape[0]):
            data[i,:] = normStandardize(data[i,:],first)

    #Square root transform, no scaling
    elif transform == 'Square root transformation' and scale == 'None':
        #square root transform
        for i in range(data.shape[0]):
            data[i,:] = sqrtTrans(data[i,:])

    #square root transform, mean centering
    elif transform == 'Square root transformation' and scale == 'Mean centering':
        #square root transform
        for i in range(data.shape[0]):
            data[i,:] = sqrtTrans(data[i,:])

        #mean center
        for i in range(data.shape[0]):
            data[i,:] = meanCentering(data[i,:])

    #square root transform, auto scaling
    elif transform == 'Square root transformation' and scale == 'Auto Scaling':
        #square root transform
        for i in range(data.shape[0]):
            data[i,:] = sqrtTrans(data[i,:])

        #auto scale
        for i in range(data.shape[0]):
            data[i,:] = standardize(data[i,:])

    #square root transform, pareto scaling
    elif transform == 'Square root transformation' and scale == 'Pareto Scaling':
        #square root transform
        for i in range(data.shape[0]):
            data[i,:] = sqrtTrans(data[i,:])

        #pareto scale
        for i in range(data.shape[0]):
            data[i,:] = paretoScaling(data[i,:])

    #square root transform, range scaling
    elif transform == 'Square root transformation' and scale == 'Range Scaling':
        #square root transform
        for i in range(data.shape[0]):
            data[i,:] = sqrtTrans(data[i,:])

        #range scale
        for i in range(data.shape[0]):
            data[i,:] = rangeScaling(data[i,:])

    #cube root transform, no scaling
    elif transform == 'Cube root transformation' and scale == 'None':
        #cube root transform
        for i in range(data.shape[0]):
            data[i,:] = cubeRtTrans(data[i,:])

    #cube root transform, mean centering
    elif transform == 'Cube root transformation' and scale == 'Mean centering':
        #cube root transform
        for i in range(data.shape[0]):
            data[i,:] = cubeRtTrans(data[i,:])

        #mean centering
        for i in range(data.shape[0]):
            data[i,:] = meanCentering(data[i,:])
    
    #cube root transform, auto scale
    elif transform == 'Cube root transformation' and scale == 'Auto Scaling':
        #cube root transform
        for i in range(data.shape[0]):
            data[i,:] = cubeRtTrans(data[i,:])

        #auto scale
        for i in range(data.shape[0]):
            data[i,:] = standardize(data[i,:])

    elif transform == 'Cube root transformation' and scale == 'Pareto Scaling':
        #cube root transform
        for i in range(data.shape[0]):
            data[i,:] = cubeRtTrans(data[i,:])

        #pareto scale
        for i in range(data.shape[0]):
            data[i,:] = paretoScaling(data[i,:])

    elif transform == 'Cube root transformation' and scale == 'Range Scaling':
        #cube root transform
        for i in range(data.shape[0]):
            data[i,:] = cubeRtTrans(data[i,:])
        
        #range scale
        for i in range(data.shape[0]):
            data[i,:] = rangeScaling(data[i,:])

    elif transform == 'None' and scale == 'Mean centering':
        #mean centering
        for i in range(data.shape[0]):
            data[i,:] = meanCentering(data[i,:])

    elif transform == 'None' and scale == 'Auto Scaling':
        #auto scaling
        for i in range(data.shape[0]):
            data[i,:] = standardize(data[i,:])

    elif transform == 'None' and scale == 'Pareto Scaling':
        #pareto scale
        for i in range(data.shape[0]):
            data[i,:] = paretoScaling(data[i,:])

    elif transform == 'None' and scale == 'Range Scaling':
        #range scale 
        for i in range(data.shape[0]):
            data[i,:] = rangeScaling(data[i,:])

    return data



