#Creating a class containing functions that will be used in GUI
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram
from scipy.cluster.hierarchy import linkage
from scipy.spatial import distance_matrix
from scipy.spatial.distance import pdist,squareform
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import minimum_spanning_tree
import glob,sys,logging,time,getpass,fpdf,os
import statistics as stat
import GuiBackground as GB
from tkinter import filedialog
from scipy.signal import argrelextrema

class GUIUtils:
    def dataIntegrity(file):
        '''
        The data integrity function checks and corrects Volcano plot outputs from MetaboAnalyst for extra decimals. 

        Input:

        file - input full path to the file, use tkinter's filedialog for ease of getting file path. 

        Output: 
        
        This function outputs an excel file with corrected values, the original file name has  _corrected appended to the end. 

        '''
        #log that the user called the data integrity function
        logging.info(': User called the Data Integrity function.')

        #Read in Volcano Plot data
        try:
            volcano = pd.read_excel(file)
        except:
            logging.error(': Failed to read in the excel file. Please put error in the Github issues tab.')
            return
            
        #grab the first row the volcano data
        check = volcano['Unnamed: 0']

        #create array that can save the fixed data and the data that did not need to be fixed
        correctedArray = np.zeros(check.shape[0])

        #search each of the volcano data rows to determine if they have double decimals.
        for i in range(check.shape[0]):
            #grab the row corresponding to the current run through the loop
            curVal = check[i]

            #reset the number of decimals to 0 before checking the string for decimal points
            decimal = 0

            #creating a string that will contain the corrected string
            corrected = ''

            #Determine if the value is a float to allow for determination of whether or not we need to check it for the appropriate value
            if isinstance(curVal,float) != True:
                #Look through the strings to find data integrity issues. 
                for j in range(len(curVal)):
                    #Check for decimals otherwise add the value to a string
                    value = curVal[j]
                    if value == '.':
                        decimal += 1
                        if decimal == 1:
                            corrected += value
                    else:
                        corrected += value
                        if j == len(curVal)-1:
                            try:
                                correctedArray[i] = float(corrected)
                            except:
                                logging.error(': Unable to convert values to floats. Make sure all data values are only contain decimals or numberic values')
                                return
            else:
                #save the data that did not need to be corrected to the correctedArray
                correctedArray[i] = curVal

        #Replace the values in the dataframe with the appropriate values
        volcano['Unnamed: 0'] = correctedArray
        del(correctedArray,i,curVal,decimal,corrected)

        finalSlash = 0
        for i in range(len(file)):
            #determine the location of the final / in the name
            if file[i] == '/':
                finalSlash = i
        file = file[finalSlash+1:len(file)-5]
        #Replace the file name with the appropriate rename
        file += '_corrected.xlsx'
        #specify the file to write to
        output = pd.ExcelWriter(file)

        #write to excel file
        volcano.to_excel(output,index=False)
        del(volcano)
        #save the excel sheet
        output.save()
        logOut = 'Updated file saved as: ' + file
        logging.info(logOut)
        #log that the data integrity function has been sucessfully completed. 
        logging.info(': Data Integrity check sucessfully completed.')
        return

    def createClustergram(norm,linkFunc,distMet):
        '''
        The function is responsible for generating the clustergrams for multivariate data. This function is capable of
        using all the linkage functions and distance measures currently implemented to the scipy.hierarchy method. OF
        note the ward-euclidean distance is the only combination available from the scipy.hierarcy package. 
        
        Input:

        norm - input is binary, 0 gives non-normalized clustergrams, 1 gives a first column normalization (first column of metabolites used as normalizing values). 

        linkFunc - input a string for the linkage function you would like to use (i.e., 'ward')

        distMet - input a string for the distance measure you would like to use (i.e., 'euclidean')

        Output:

        This function outputs a .png of the generated clustergram. 
        '''

        #log that the user called the Create Clustergram function
        logging.info(': User called the Create Clustergram Function.')
        #check that the file the user selects is appropriate
        metab_data = GB.fileCheck()
        if metab_data is None:
            #log error message and return for soft exit.
            logging.error(': Error loading in the Excel sheet.')
            return

        #read in data
        data = GB.readInColumns(metab_data)

        #Standardize the data before clustering the results
        logging.info(': Standardizing the data.')
        for i in range(metab_data.shape[0]):
            data[i,:] = GB.standardize(data[i,:])
        del(metab_data)

        #create dendrogram and plot data
        GB.create_dendrogram(data,norm, link=linkFunc, dist=distMet)

        del(data,norm,linkFunc,distMet)

        logging.info(': Sucessfully created the wanted Clustergram')
        return

    def groupMedians():
        '''
        Determine the number of groups and then create a list or array of the appropriate
        beginning and ending of each group. This assumes that the groups are all of equal size which should be
        the goal for any and all analysis. Groups with out the same sizes should be considered
        inappropriate for analysis in this context, additionally it should be noted that statistics
        with out the same groups sizes can lead to incorrect analysis.
        
        Input:
        
        groupMedians does not accept any inputs, it will prompt you to select a file which you would like to have medians determined for. 

        Output:

        This function outputs a .csv file with _Medians appended to the end of the original file name. *** This will soon be updated to
        an excel file output for easy calling and input to the createClustergram, MST or Ensemble clustering functions. 
        '''
        #log that the user called the group medians function
        logging.info(': User called the Group Medians function.')
        file = filedialog.askopenfilename()
        medians = GB.fileCheck(file =file)
        if medians is None:
            #log error message and return the function for soft exit.
            logging.error(': Error reading in the excel sheet')
            return

        #get the group letters to know which groups to find the medians of
        groups = medians['mz']

        numObs = int(medians.shape[0]-1)
        

        num_groups = int(medians['Unnamed: 0'][numObs][0])

        factor = len(medians['Unnamed: 0'])/num_groups

        #create a numpy array containing 7 columns to allow for input of the m/z values and the groups
        mediansOut = np.zeros((medians.shape[1]-2,num_groups+1))
            
        #populate the first column of the array with the m/z values
        print(medians.columns[2:medians.shape[1]])
        mediansOut[:,0] = medians.columns[2:medians.shape[1]]
        print('Good')
        #Get the medians for each group and metabolite
        for i in range(num_groups):
            #calculate the start and end for each set of median calculations
            start = int(factor*i)
            end = int((factor*(i+1)))
            for j in range(mediansOut.shape[0]):
                #find the median for the first groups of values
                curMean = stat.median(medians[medians.columns[j+2]][start:end])
                #put medians into the appropriate table
                mediansOut[j,i+1] = curMean
        del(curMean,start,end,i,j,medians)
        #create list contains the headers for the files
        medianList = ['m/z']
        for i in range(num_groups):
            medianList.append('M' +str(i+1))
        #create dictionary that contains the data with there appropriate headers to input to a dataframe and
        #then be saved to a csv file
        medianDict = {}
        for i in range(num_groups+1):
            #input the appropriate data and key to the dictionary
            medianDict[medianList[i]] = mediansOut[:,i]
        del(mediansOut)
        #create dataframe that prepares the data to be input to a csv file
        mediansCSV = pd.DataFrame(data=medianDict)

        finalSlash = 0
        for i in range(len(file)):
            #determine the location of the final / in the name
            if file[i] == '/':
                finalSlash = i

        file = file[finalSlash+1:len(file)-5]
        file += '_Medians.csv'

        #specify the file that I want the program to write to.
        mediansCSV.to_csv(file,columns=medianList,index=False)
        del(file,medianList,medianDict,mediansCSV)
        #logging the completion of the group medians function
        logging.info(': Sucessfully grouped the Medians of each group!')
        return

    def linkageComparison(file,num_comps,linkList,distance):
        '''
        Compares 2-4 linkage functions on a given set of data. 
        
        linkageComparison requires a file, number of comparisons, and a list of linkage functions. 

        Input:

        file - include full file path, use the tkinter filedialog functionality for ease of obtaining file path

        num_of_comps - make sure to give an integer the same length as the link list. 
            
        linkList - list of linkage functions that you would like to have compared. 

        Output:

        linkageComparison saves a .png file of the output to the current working directory. 
        '''

        #set recursion limit above the common max for our data.
        sys.setrecursionlimit(10**8)
        #Log that user called linkage comparison function
        logging.info(': User called the Linkage Comparison function.')
        #check that the file is appropriate for our data set
        metab_data = GB.fileCheck(file)   

        if metab_data is None:
            #Logs error and returns function to ensure soft exit.
            logging.error(': Error loading in excel file check log file!')
            return

        #read in column data
        data = GB.readInColumns(metab_data)

        #Standardize the data before clustering the results
        logging.info('Standardizing the data.')
        for i in range(metab_data.shape[0]):
            data[i,:] = GB.standardize(data[i,:])
        del(i)
        #convert string to integer
        num_comps = int(num_comps)
        
        if num_comps == 2:
            #Create the linkage matrix
            linkageOne = linkage(data,linkList[0], metric=distance)
            distMeasure = pdist(data)
            distMeasure = squareform(distMeasure)
            linkageTwo = linkage(data,linkList[1], metric=distance)

            #Create the appropriate plt figure to allow for the comparison of linkage functions
            fig, axes = plt.subplots(1,2,figsize=(8,8))

            #create the dendrograms
            dend1 = dendrogram(linkageOne,ax=axes[0],above_threshold_color='y',orientation='left',no_labels=True)
            dend2 = dendrogram(linkageTwo,ax=axes[1],above_threshold_color='y',orientation='left',no_labels=True)
            del(linkageOne,linkageTwo,num_comps)
            print('Good')
        elif num_comps == 3:
            #Create the linkage matrix
            linkageOne = linkage(data,linkList[0],metric=distance)
            linkageTwo = linkage(data,linkList[1],metric=distance)
            linkageThree = linkage(data,linkList[2], metric=distance)

            #Create the appropriate plt figure to allow for the comparison of linkage functions
            fig, axes = plt.subplots(1,3,figsize=(8,8))

            #create the dendrograms
            dend1 = dendrogram(linkageOne,ax=axes[0],above_threshold_color='y',orientation='left',no_labels=True)
            dend2 = dendrogram(linkageTwo,ax=axes[1],above_threshold_color='y',orientation='left',no_labels=True)
            dend3 = dendrogram(linkageThree,ax=axes[2],above_threshold_color='y',orientation='left',no_labels=True)
            del(linkageOne,linkageTwo,linkageThree,num_comps)
            print('Good')
        elif num_comps == 4:
            #Create the linkage matrix
            linkageOne = linkage(data,linkList[0],metric=distance)
            linkageTwo = linkage(data,linkList[1],metric=distance)
            linkageThree = linkage(data,linkList[2],metric=distance)
            linkageFour = linkage(data, linkList[3],metric=distance)

            #Create the appropriate figure to allow for the comparison of linkage functions
            fig, axes = plt.subplots(2,2,figsize=(8,8))
            plt.title('Linkage Comparison')

            #create the dendrograms
            dend1 = dendrogram(linkageOne,ax=axes[0,0],above_threshold_color='y',orientation='left',no_labels=True)
            dend2 = dendrogram(linkageTwo,ax=axes[0,1],above_threshold_color='y',orientation='left',no_labels=True)
            dend3 = dendrogram(linkageThree,ax=axes[1,0],above_threshold_color='y',orientation='left',no_labels=True)
            dend4 = dendrogram(linkageFour,ax=axes[1,1],above_threshold_color='y',orientation='left',no_labels=True)
            del(linkageOne,linkageTwo,linkageThree,linkageFour,num_comps)
            print('Good')
        linkPre = 'LinkageComparison'
        linkSuf = '.png'
        sep = '_'
        firstCheck = linkPre+sep
        for i in range(len(linkList)):
            #create the first file check
            firstCheck += linkList[i] + sep

        firstCheck += '01' + linkSuf

        chkBuffer = glob.glob("*.png")
        count = 1
        if firstCheck in chkBuffer:
            checkVal = False
            firstCheck = firstCheck.strip(linkSuf)
            firstCheck = firstCheck.strip('01')
            while checkVal == False:
                count += 1
                #search the "buffer" for ensemble cluster
                if count < 10:
                    #determine if the file has already been made
                    curFileCheck = firstCheck + '0' + str(count) + linkSuf
                    if curFileCheck not in chkBuffer:
                        checkVal = True
                        linkFile = curFileCheck

                else:
                    curFileCheck = firstCheck + str(count) + linkSuf
                    if curFileCheck not in chkBuffer:
                        checkVal = True
                        linkFile = curFileCheck
            plt.savefig(linkFile,dpi=600)
        else:
            linkFile = firstCheck 
            plt.savefig(linkFile,dpi=600)
        print('Still good')
        plt.show()

        #log the completion of the linkage comparison
        logging.info(': Sucessfuly completed the comparison of the linkage functions!')
        return
            
    def compoundMatchUp():
        '''
        The compoundMatchUp function is responsible for matching the output compounds from mummichog to compounds from the KEGG data base spreadsheet. 

        Input:

        compoundMatchUp does not allow input. The file "mummichog_matched_compound_all.csv" needs to be updated prior to running this function. 

        *** Stay tuned this function will be updated soon. 
        '''
        logging.info(': Compound Match-Up function called!')
        # Reads in our Kegg Compound Dataset (Single Column)
        kegg_data = pd.read_excel("kegg_compound_IDs_3.xlsx")

        # Splits our single column into two more user friendly ones
        kegg_data[["ID","compound"]] = kegg_data["ID"].str.split(" ", 1, expand = True)

        # Pulls in our matched compound data
        #ask user for file input and read in the csv file. 
        file = filedialog.askopenfilename()
        my_data = pd.read_csv(file)

        # Makes an ID column
        my_data["ID"] = my_data["Matched.Compound"]

        # Deletes the unneeded columns
        gonecolumns = ["Query.Mass", "Matched.Form", "Mass.Diff", "Matched.Compound"]
        my_data = my_data.drop(axis = 1, labels = gonecolumns)

        # Filters the keggs data to only the ids that are
        # in the matched compound dataset
        my_final_data = kegg_data[kegg_data["ID"].isin(my_data["ID"])]

        # Writes this final dataset to a csv
        my_final_data.to_csv(path_or_buf = "CompoundMatchups.csv")
        return

    def ensembleClustering(optNum=2):
        '''
        The distance measures and linkage functions should be consistent but we could also develop
        a GUI that allows for the users to select various distance measures. The linkage functions 
        should be consistent for all ensemble clustering techniques

        Ensemble clustering generates the ensemble average of the clusterings from 13 different clusterings of the data. Please ask Brady Hislop
        if you have questions about how Ensemble clustering works. I will be happy to share notes and/or have a conversation with you about the 
        process that is entailed in generating these ensemble averages. 

        Input:

        optNum - Input the optimum number of clusters for these data based upon a minimum spanning tree optimization of the data. 

        Output:

        A figure output by these data will be saved as a .png to the current working directory. Additionally, the red-dashed lines around the 
        yellow portions of the graph represent the regions of metabolites which were clustered together 13 out of 13 times. Each of these boxed lines will
        be out as csv files, along with a csv file of the CoOccurence matrix. 

        '''
        #log that the user called ensemble clustering function
        logging.info(': User called Ensemble Clustering function.')

        #optimum number of clusters from validation index.
        sys.setrecursionlimit(10**8)

        #Make sure file can be read in. 
        metab_data = GB.fileCheck()
        if metab_data is None:
            logging.error(': File to not meet input requirements.')
            return

        #List for the use in creating and plotting the clustering results
        linkageList = ['single','complete','average']
        distList = ['euclidean','sqeuclidean','chebyshev','seuclidean']  #,'cosine']

        #calculate the number of clusterings based upon the size of the lists and an additional term for the ward-euclidean run. 
        numClusterings = (len(linkageList)*len(distList))+1

        #read in the data
        data = GB.readInColumns(metab_data)

        #Standardize the data before clustering the results
        logging.info(': Standardizing data.')
        for i in range(data.shape[0]):
            data[i,:] = GB.standardize(data[i,:])

        #determine the the number of clusters and the dictionary location that needs to be called. 
        numMetabs = data.shape[0]
        dictLoc = numMetabs-optNum-1

        #create co-occurrence matrix.
        coOcc = GB.cooccurrence(data)

        for i in range(len(linkageList)):
            for j in range(len(distList)):
                start = time.perf_counter()
                linkCur = linkage(data,linkageList[i],distList[j])
                valid = GB.clustConnectLink(linkCur)
                coOcc = GB.popCooccurrence(valid[dictLoc],coOcc,numClusterings)
                end = time.perf_counter()
                logging.info(str(linkageList[i])+'-'+str(distList[j]) +' done!')
                print(end-start)
        del(linkageList,distList)

        #add a ward euclidean clustering to the ensemble. 
        start = time.perf_counter()
        linkCur = linkage(data,'ward','euclidean')
        valid = GB.clustConnectLink(linkCur)
        coOcc = GB.popCooccurrence(valid[dictLoc],coOcc,numClusterings)
        end = time.perf_counter()
        logging.info('Ward-Euclidean done!')
        print(end-start)

        #create the ensemble dendrogram using ward-euclidean inputs. 
        GB.createEnsemDendrogram(coOcc,metab_data,norm=0,link='ward',dist='euclidean',func="ensemble")
        GB.plotting()

        #make the coOccurence matrix a dataframe.
        coOcc = pd.DataFrame(coOcc)
        try:
            #try to save the large .csv file of the CoOccurence matrix.
            coOcc.to_csv('EnsembleCoOcc.csv',index=False)
        except:
            logging.error(': Failed to save the Ensemble CoOccurence matrix!!')
            
        #Log the completgroupion of the ensemble clustering function
        logging.info(': Sucessfully completed Ensemble clustering!')
        return

    def MST(func = "base"):
        #*******Need for dataMST and mstoutNP??

        '''
        MST generates a minimum spanning tree of input data, and then validates the optimum number of clusters based upon a validation index of 
        the ***intra/inter*** cluster distances.

        Input:

        MST doesn't accept inputs, and will prompt you for an input file. 

        Output:
        
        MST will output two csv files, one with the indicies of the minimum spanning tree connections, the other will be contain a csv file with
        the the cluster number and the validation index value - local minimum is the optimum number of clusters. 
        '''

        #log that user called MST
        logging.info(': User called Minimum Spanning Tree function.')
        #check the stability of the file
        dataRaw = GB.fileCheck()
        if dataRaw is None:
            #log error and return function to ensure a soft closing of the class
            logging.error(': Error loading the Excel sheet.')
            return

        #read in raw data
        data = GB.readInColumns(dataRaw)
        #find the number of groups in that data
        num_groups = dataRaw.shape[1]-2
        
        #standardize the data before clustering
        logging.info(': Standardizing the Data.')
        for i in range(dataRaw.shape[0]):
            data[i,:] = GB.standardize(data[i,:])

        #find the distance matrix using the pairwise distance function, put into squareform, appropriate format for mst and submit. 
        pairWise = pdist(data)
        pairWise = squareform(pairWise)
        mstInput = csr_matrix(pairWise)
        mstOut = minimum_spanning_tree(mstInput)
        #get out the non-zero indicies. 
        mstOutInd = mstOut.nonzero()

        #create the matrix containing the various connections of the mst
        mstOutMat = mstOut.todense()
        dataMST = np.zeros([data.shape[0]-1,3])

        for i in range(data.shape[0]-1):
            #input the values of each matrix element to the numpy array for saving to csv file.
            dataMST[i,0] = mstOutInd[0][i]
            dataMST[i,1] = mstOutInd[1][i]
            dataMST[i,2] = mstOut[dataMST[i,0],dataMST[i,1]]

        #Input the dataMST into the dataframe to save the results of the MST for future use if needed
        mstOut = pd.DataFrame(dataMST, columns=['index1','index2','dist'])
        mstOut = mstOut.sort_values(by='dist')
        mstOutNp = mstOut.to_numpy()

        mstOutMat = pd.DataFrame(mstOutMat)
 
        #determine how the minimum spanning tree was created for validation of clusters
        validationClusters = GB.clustConnect(dataMST,mstOutNp)

        #Validate the number of clusters that should be used in the clustering solutions.
        valIndex = GB.Validate(validationClusters,data,num_groups)

        x=valIndex[1,:]
        y=valIndex[0,:]
        #find the local minimums
        minimums = argrelextrema(y,np.less)

        #determine the minimums y values to for plotting
        miniVals = np.zeros((len(minimums[0]),2))
        for i in range(len(minimums[0])):
            #input the number of clusters and the validation index value.
            miniVals[i,0] = x[minimums[0][i]]
            miniVals[i,1] = y[minimums[0][i]]
        ind = np.unravel_index(np.argmin(miniVals, axis=None), miniVals.shape)
        print(ind)
        minValIndex = miniVals[ind[0],:]

        valOut = np.zeros((2,valIndex.shape[1]))
        for j in range(valIndex.shape[1]):
            #flip the array so that lower clusters start at lower positions
            valOut[0,j] = valIndex[0,valIndex.shape[1]-j-1]
            valOut[1,j] = valIndex[1,valIndex.shape[1]-j-1]

        valIHeaders = list(valOut[1,:])
        valIndex = pd.DataFrame(valOut,columns=valIHeaders)
        valIndex = valIndex.drop(1,axis=0)
        rowLabels = ["Validation Index"]
        valIndex.insert(0,"Clusters",rowLabels)

        #save to a csv file
        mstOut.to_csv('MST_branches.csv',index=False)

        #save validation measure to csv file
        valIndex.to_csv('valIndex.csv',index=False)
        if func == "base":
            #logging the completion of the Minimum spanning tree
            logging.info(': Sucessfully completed MST and clustering validation!')
            plt.plot(valOut[1,:],valOut[0,:])
            plt.plot(minValIndex[0],minValIndex[1],'r*')
            font = {'family': 'serif','color':  'black','weight': 'bold','size': 20}
            plt.text(valIndex.shape[1]/2, 0.75, str(int(minValIndex[0]))+' - Clusters!!', fontdict=font)
            plt.xlabel('Clusters')
            plt.ylabel('Validation Index')
            plt.title('Cluster Validation')
            plt.show()

        return minValIndex[0]

    def peaksToPathways():
        '''
        Create input files for the mummichog algorithm, using files output from the ensemble clustering and the in the future from the clustergram function. 

        Input:

        peaksToPathways does not accept any inputs but will prompt the user for two inputs. First, the user will need to select the original files for there ensemble clustered data. 
        Next, the user will need to select the directory containing the files generated by the ensemble clustering function. 

        Output:

        This function will output csv files containing the m/z value and p-values fo the matched metabolites (p-values =0.04), and the remaining the metabolites with p-values equal to 1. 
        '''
        logging.info(': Entering the Peaks to Pathways generator!')
        #ask user to input the file name of the original data
        filename = filedialog.askopenfilename()

        dataRaw = GB.fileCheck(file = filename)
        if dataRaw is None:
            #log error and return function to ensure a soft closing of the class
            logging.error(': Error loading the reference Excel sheet.')
            return

        #ask user to select the directory containing the csv files from ensemble clustering output (currently only method available)
        direct = filedialog.askdirectory()

        #create string to check for csv files
        globCheck = direct + "/*.csv"
        files = glob.glob(globCheck)

        nameCheckStriper = direct +'\\'
        ensemFiles = []
        for i in range(len(files)):
            #strip the beginning of the strip off and then check the first 5 characters of the stripped string
            curCheck = files[i].strip(nameCheckStriper)
            if curCheck[0:5] == 'Ensem':
                #append to ensemFiles
                ensemFiles.append(direct + '/' + curCheck)

        columnsHead = list(dataRaw.columns)
        columnFirst = columnsHead[0]
        for i in range(len(ensemFiles)):
            dataClust = np.ones((dataRaw.shape[0],3))
            dataClust[:,0] = dataRaw[columnFirst]
            dataClust[:,1] = dataRaw["rtmed"]
            #start the process of reading in and creating the ensemble output files. 
            try:
                dataCur = pd.read_csv(ensemFiles[i])
            except:
                logging.error(': Failed to read in the excel sheet check the sheet name is Medians')
                return

            dataMzRt = np.zeros((dataCur.shape[0],2))
            dataMzRt[:,0] = dataCur["Identities"]
            dataMzRt[:,1] = dataCur["rt_med"]

            for j in range(dataMzRt.shape[0]):
                #determine the location of the m/z values in the dataClust
                locRT =  np.where(abs(dataClust[:,1]-dataMzRt[j,1]) < 0.0001)
                if len(locRT[0]) > 1:
                    for k in range(len(locRT[0])):
                        if abs(dataClust[locRT[0][k],0] - dataMzRt[j,0]) < 0.0001:
                            dataClust[locRT[0][k],2] = 0.04

                elif len(locRT[0]) == 1:
                    dataClust[locRT[0][0],2] = 0.04

                else:
                    print(j+1)
                    logging.warning(': Creation of peaks to pathway files halted due to non-matching values, please make sure you have selected appropriate reference file.')
                    return
            #create the files that can be submitted to the csv saving file. 
            dataOut = np.zeros((dataRaw.shape[0],2))

            dataOut[:,0] = dataClust[:,0]
            dataOut[:,1] = dataClust[:,2]

            dataOut = pd.DataFrame(dataOut,columns=["m.z","p.value"])
            dataOut = dataOut.sort_values(by=["p.value"])

            p2pPre = 'PeaksToPathways'
            p2pSuf = '.csv'
            firstCheck = p2pPre + '01' + p2pSuf

            chkBuffer = glob.glob("*.csv")
            count = 1
            if firstCheck in chkBuffer:
                checkVal = False
                while checkVal == False:
                    count += 1
                    #search the "buffer" for ensemble cluster
                    if count < 10:
                        #determine if the file has already been made
                        curFileCheck = p2pPre + '0' + str(count) + p2pSuf
                        if curFileCheck not in chkBuffer:
                            checkVal = True
                            p2pFile = curFileCheck

                    else:
                        curFileCheck = p2pPre + str(count) + p2pSuf
                        if curFileCheck not in chkBuffer:
                            checkVal = True
                            p2pFile = curFileCheck
                dataOut.to_csv(p2pFile, index=False)
            else:
                p2pFile = p2pPre + '0'+ str(count) + p2pSuf 
                dataOut.to_csv(p2pFile, index=False)
            logging.info(':Success!')
        logging.info(': Leaving the Peaks to Pathways Function!')
        return

    def PDFGenerator():
        '''
        Generates PDF of the MetaboAnalystR results. 

        Input:

        PDFGenerator does not accept input, but asks for the directory containing the output images. 

        Output:

        PDF report of the results from a MetaboAnalystR run. 

        '''

        #log that the function has been called
        logging.info(': Entering PDF Generator function.')

        #create the pdf and title for each page.
        pdf = fpdf.FPDF('P','mm','Letter')

        #Create the title and set the default font
        directory = filedialog.askdirectory()
        os.chdir(directory)
        #determine the current user
        curUser = getpass.getuser()
        curUser = GB.who(curUser)

        #Create the first page
        title = 'Metabolanalyst Results' + '-' + curUser
        pdf.add_page()
        pdf.set_font('Arial','B',24)
        pdf.cell(197,10,title,0,0,'C')
        pdf.set_font('Arial','B',14)
        pdf.set_font('')
        pdf.ln(10)
        #*******************************************
        #create the variability in the pdf's here.
        #*******************************************
        files = GB.files(directory)
        #first page should always contain the normalized and sample normalizations
        #for the first iteration look for the Normalization and sample normalization
        norm = 'norm_0_dpi300.png' in files
        snorm = 'snorm_0_dpi300.png' in files

        if norm is True and snorm is True:
            #input the sample and full data normalization onto the first page of pdf document
            pdf.cell(197,10,'Normalization',0,0,'L')
            pdf.ln(10)
            imRatio = GB.imageSize('norm_0_dpi300.png') 
            pdf.image('norm_0_dpi300.png',55,30,100*imRatio,100)
            
            pdf.ln(120)
            pdf.cell(197,10,'Sample Normalization',0,0,'L')
            pdf.ln(10)
            imRatio = GB.imageSize('snorm_0_dpi300.png')
            pdf.image('snorm_0_dpi300.png',55,160,100*imRatio,100)

            #remove the first two files used above from the list
            files.remove('norm_0_dpi300.png')
            files.remove('snorm_0_dpi300.png')

        #determine the number of pages needed in the pdf report to generate the appropriate for loop.
        pages = int((len(files))/2)
        for i in range(pages):
            # iterating through the images to get create a pdf
            pdf.add_page()

            # grab the first file that is available and send it to the image size function and then the naming function
            fileOne = files[0]
            imRatio = GB.imageSize(fileOne)
            headerOne = GB.pdfHeader(fileOne)

            #add the first image to the page.
            pdf.cell(197,10,headerOne,0,0,'L')
            pdf.ln(10)
            pdf.image(fileOne,55,30,100*imRatio,100)

            # grab the second file that is available and send it to the image size function and then the naming function
            fileTwo = files[1]
            imRatio = GB.imageSize(fileTwo)
            headerTwo = GB.pdfHeader(fileTwo)

            #add the second image to the page
            pdf.ln(120)
            pdf.cell(197,10,headerTwo,0,0,'L')
            pdf.ln(10)
            pdf.image(fileTwo,55,160,100*imRatio,100)

            #delete the first two files from the list of files
            files.remove(fileOne)
            files.remove(fileTwo)

        #create the pdf of the results
        ending = '.pdf'
        fileName = ''
        curTime = time.strftime("%a_%b_%d_%Y_%H_%M")
        fileName += curUser + '_' + curTime + ending
        pdf.output(fileName,'F')
        #log the sucessful creation of the pdf
        logging.info(': Sucessfully created a pdf of the results!')
        logging.info(': Leaving the pdf PDF Generator Function!')
        return