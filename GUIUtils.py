#Creating a class containing functions that will be used in GUI
import re
from numpy.lib.arraysetops import isin
import pandas as pd
import numpy as np
from matplotlib import cm, pyplot as plt
from scipy.cluster.hierarchy import dendrogram
from scipy.cluster.hierarchy import linkage
from scipy.spatial import distance_matrix
from scipy.spatial.distance import pdist,squareform
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import minimum_spanning_tree
from scipy.stats import t
from scipy.stats import normaltest
import scipy.stats as stats
import glob,sys,logging,time,getpass,fpdf,os
import statistics as stat
from multiprocessing import Pool
import seaborn as sns
import config
import random
import math

from sklearn import cluster
import GuiBackground as GB
from tkinter import filedialog, messagebox
from scipy.signal import argrelextrema,argrelmin, argrelmax
import mplcursors
from LocallyWeighted import LocallyWeighted as LW 
import ValidationMetric

from Bio.KEGG import REST
from Bio.KEGG import Compound

from ValidationMetric import ValidationMetric as VM

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

        try:
            #Read in Volcano Plot data
            if file[len(file)-1] == 'x':
                volcano = pd.read_excel(file)
            elif file[len(file)-1 == 'v']:
                volcano = pd.read_csv(file)
        except:
            logging.error(': Failed to read in the excel file. Please put error in the Github issues tab.')
            messagebox.showerror(title='Error',message='Failed to read in the excel file. Please let Brady know!!')
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
                    if decimal == 2:
                        correctedArray[i] = corrected
                        continue
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
        messagebox.showinfo(title="Success",message="Removed data integrity issues!!")
        return

    def createClustergram(norm,linkFunc,distMet,cmap, transform = 'None', scale ='None'):
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
        logging.info(':-------------------------------------------------------------')
        logging.info(': User called the Create Clustergram Function.')
        logMessage = ': Linkage Function:' + linkFunc
        logging.info(logMessage)
        logMessage = ': Distance Metric:' + distMet
        logging.info(logMessage)
        logMessage = ': Data Transform: ' + transform +'; Data Scaling: ' + scale
        logging.info(logMessage)

        try:
            data, col_groups = GB.readAndPreProcess(file='',transform=transform,scale=scale,func="CC")
        except TypeError:
            logging.error(': No file selected!')
            messagebox.showerror(title='Error',message='No file selected, returning to GUI. If you wish to continue with creating a clustergram, click continue and submit!')
            return

        #create dendrogram and plot data        
        GB.create_dendrogram(data,col_groups, norm, link=linkFunc, dist=distMet,color = cmap)

        del(data,norm,linkFunc,distMet)

        logging.info(': Sucessfully created the wanted Clustergram')
        return

    def groupMedians(rmZeros=0):
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
        medians = GB.fileCheck(file=file)
        if medians is None:
            #log error message and return the function for soft exit.
            logging.error(': Error reading in the excel sheet')
            return

        #get the group letters to know which groups to find the medians of
        groups = medians['mz']

        metaboliteIdentities = list(medians.columns)
        metaboliteIdentities = metaboliteIdentities[2:len(metaboliteIdentities)]

        #put groups into a list so we know how many groups we need to find the median for.
        groups = set(groups)
        groups = sorted(groups)

        #get the groupID for each sample
        samples = list(medians['mz'])

        #create numpy array that contains the medians...
        mediansOut = np.zeros((len(metaboliteIdentities),len(groups)+2))


        for i in range(len(groups)):
            #create list of indicies for the program to obtain
            curInd = []
            numFound = 0
            #get theh current number of times a particular group is present.
            numLocs = samples.count(groups[i])
            for j in range(len(samples)):
                #find the numLocs
                curEntry = samples[j]
                if curEntry == groups[i]:
                    #add one to the numFound
                    numFound += 1
                    curInd.append(j)
                    if numFound == numLocs:
                        break

            #determine the median for each metabolite across the samples for each group
            for j in range(len(metaboliteIdentities)):
                #for each metabolite grab the current list of data
                curMetabMed = []
                for k in range(len(curInd)):
                    #add values too list that are from the list of values
                    curMetabMed.append(medians[metaboliteIdentities[j]][curInd[k]])

                curMedian = stat.median(curMetabMed)
                mediansOut[j,i+1] = curMedian

        for i in range(len(metaboliteIdentities)):
            #check if theh current value is a string, if it is then fix string. 
            if isinstance(metaboliteIdentities[i], str):
                #correct the string to a float.
                #find the first instance of the decimal, then search for values after
                curString = metaboliteIdentities[i]
                first = curString.find('.')
                startRemove = curString.find('.',first+1,len(curString))
                updatedString = curString[0:startRemove]
                curVal = float(updatedString)
                mediansOut[i,0] = curVal

            else:
                mediansOut[i,0] = metaboliteIdentities[i]

        if rmZeros == 1:
            #create a numpy array that cna continually be added to. 
            mediansNZ = np.zeros((1,mediansOut.shape[1]))

            #remove the zeros from the data sheet, allowing the final data sheet to only contain detected metabolites.
            #number of groups to determine whether or not to keep the column. 
            numGroups = len(groups)
            for i in range(len(metaboliteIdentities)):
                #get the current row from the numpy array
                curRow = mediansOut[i,1:mediansOut.shape[1]-1]
                curZero = np.where(curRow == 0)
                numZero = len(curZero[0])
                if numZero/numGroups < 0.5:
                    #add the row to the output of medians NZ
                    newRow = np.zeros((1,mediansOut.shape[1]))
                    newRow[0,0] = mediansOut[i,0]
                    newRow[0,newRow.shape[1]-1] = mediansOut[i,mediansOut.shape[1]-1]
                    newRow[0,1:newRow.shape[1]-1] = curRow
                    mediansNZ = np.vstack([mediansNZ, newRow])

            mediansNZ = np.delete(mediansNZ,0,0)
            mediansOut = mediansNZ
            logging.info(': Zeros removed from the data set!')
        else:
            logging.info(': Zeros are not being removed!')


        columns =['m/z']
        columns.extend(groups)
        columns.append('rt_med')

        mediansExcel = pd.DataFrame(data=mediansOut,columns=columns)
        file = file[0:len(file)-5]
        file = file + 'Medians.xlsx'
        mediansExcel.to_excel('MediansOutput.xlsx',index=False,sheet_name="Medians")

        #logging the completion of the group medians function
        logging.info(': Successfully grouped the Medians of each group!')
        messagebox.showinfo(title="Success",message="Successfully created MediansOutput.xlsx file!!")
        return

    def linkageComparison(file,num_comps,linkList,distance, transform,scale):
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

        data, col_groups = GB.readAndPreProcess(file =file,transform=transform,scale=scale,func="CC")
        del(col_groups)
        #input the arguments to the log file so user has record of what was input.
        logging.info(':-------------------------------------------------------------')
        logMessage = file
        logging.info(logMessage)
        logMessage = ': Number of comparisons: ' + str(num_comps)
        logging.info(logMessage)
        logMessage = ': Linkage functions: ' + str(linkList)
        logging.info(logMessage)
        logMessage = ': Distance metric: ' + distance
        logging.info(logMessage)
        logMessage = ': Data Transform: ' + transform +'; Data Scaling: ' + scale
        logging.info(logMessage)

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
            axes[0].set_title(linkList[0])
            axes[1].set_title(linkList[1])
            #grab the last entry of the linkage list
            maxList = np.zeros((1,2))
            maxList[0,0] = linkageOne[len(linkageOne)-1][0]
            maxList[0,1] = linkageOne[len(linkageOne)-1][1]
            maxLinkNum = int(np.amax(maxList))
            sameColor = []
            for i in range(maxLinkNum+2):
                sameColor.append('k')
            #create the dendrograms
            dend1 = dendrogram(linkageOne,ax=axes[0],above_threshold_color='y',orientation='left',no_labels=True, link_color_func= lambda x: sameColor[x])
            dend2 = dendrogram(linkageTwo,ax=axes[1],above_threshold_color='y',orientation='left',no_labels=True, link_color_func= lambda x: sameColor[x])
                
            del(linkageOne,linkageTwo,num_comps)
        elif num_comps == 3:
            #Create the linkage matrix
            linkageOne = linkage(data,linkList[0],metric=distance)
            linkageTwo = linkage(data,linkList[1],metric=distance)
            linkageThree = linkage(data,linkList[2], metric=distance)

            #Create the appropriate plt figure to allow for the comparison of linkage functions
            fig, axes = plt.subplots(1,3,figsize=(8,8))

            axes[0].set_title(linkList[0])
            axes[1].set_title(linkList[1])
            axes[2].set_title(linkList[2])

            #grab the last entry of the linkage list
            maxList = np.zeros((1,2))
            maxList[0,0] = linkageOne[len(linkageOne)-1][0]
            maxList[0,1] = linkageOne[len(linkageOne)-1][1]
            maxLinkNum = int(np.amax(maxList))
            sameColor = []
            for i in range(maxLinkNum+2):
                sameColor.append('k')
            #create the dendrograms
            dend1 = dendrogram(linkageOne,ax=axes[0],above_threshold_color='y',orientation='left',no_labels=True, link_color_func= lambda x: sameColor[x])
            dend2 = dendrogram(linkageTwo,ax=axes[1],above_threshold_color='y',orientation='left',no_labels=True, link_color_func= lambda x: sameColor[x])
            dend3 = dendrogram(linkageThree,ax=axes[2],above_threshold_color='y',orientation='left',no_labels=True, link_color_func= lambda x: sameColor[x])
            del(linkageOne,linkageTwo,linkageThree,num_comps)
        elif num_comps == 4:
            print(data)

            #Create the linkage matrix
            linkageOne = linkage(data,linkList[0],metric=distance)
            linkageTwo = linkage(data,linkList[1],metric=distance)
            linkageThree = linkage(data,linkList[2],metric=distance)
            linkageFour = linkage(data, linkList[3],metric=distance)

            #Create the appropriate figure to allow for the comparison of linkage functions
            fig, axes = plt.subplots(2,2,figsize=(8,8))

            axes[0,0].set_title(linkList[0],fontsize=24)
            axes[0,1].set_title(linkList[1],fontsize=24)
            axes[1,0].set_title(linkList[2],fontsize=24)
            axes[1,1].set_title(linkList[3],fontsize=24)

            #grab the last entry of the linkage list
            maxList = np.zeros((1,2))
            maxList[0,0] = linkageOne[len(linkageOne)-1][0]
            maxList[0,1] = linkageOne[len(linkageOne)-1][1]
            maxLinkNum = int(np.amax(maxList))
            sameColor = []
            for i in range(maxLinkNum+2):
                sameColor.append('k')

            #create the dendrograms
            dend1 = dendrogram(linkageOne,ax=axes[0,0],above_threshold_color='y',orientation='left',no_labels=True, link_color_func= lambda x: sameColor[x])
            dend2 = dendrogram(linkageTwo,ax=axes[0,1],above_threshold_color='y',orientation='left',no_labels=True, link_color_func= lambda x: sameColor[x])
            dend3 = dendrogram(linkageThree,ax=axes[1,0],above_threshold_color='y',orientation='left',no_labels=True, link_color_func= lambda x: sameColor[x])
            dend4 = dendrogram(linkageFour,ax=axes[1,1],above_threshold_color='y',orientation='left',no_labels=True, link_color_func= lambda x: sameColor[x])
            del(linkageOne,linkageTwo,linkageThree,linkageFour,num_comps)
        elif num_comps == 1:
            #Create the linkage matrix
            linkageOne = linkage(data,linkList[0], metric=distance)
            distMeasure = pdist(data)
            distMeasure = squareform(distMeasure)

            #Create the appropriate plt figure to allow for the comparison of linkage functions
            fig, axes = plt.subplots(1,1,figsize=(8,8))
            axes.set_title(linkList[0],fontsize=24)
            #grab the last entry of the linkage list
            maxList = np.zeros((1,2))
            maxList[0,0] = linkageOne[len(linkageOne)-1][0]
            maxList[0,1] = linkageOne[len(linkageOne)-1][1]
            maxLinkNum = int(np.amax(maxList))
            sameColor = []
            for i in range(maxLinkNum+2):
                sameColor.append('k')

            dend1 = dendrogram(linkageOne,ax=axes,above_threshold_color='y',orientation='left',no_labels=True, link_color_func= lambda x: sameColor[x])
                

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
            plt.savefig(linkFile,dpi=600,transparent=True)
        else:
            linkFile = firstCheck 
            plt.savefig(linkFile,dpi=600,transparent=True)

        plt.show()

        #log the completion of the linkage comparison
        logging.info(': Sucessfuly completed the comparison of the linkage functions!')
        return
            
    def compoundMatchUp(typeFile = 'all'):
        '''
        The compoundMatchUp function is responsible for matching the output compounds from mummichog to compounds from the KEGG data base spreadsheet. 

        Input:

        compoundMatchUp does not allow input. The file "mummichog_matched_compound_all.csv" needs to be updated prior to running this function. 

        *** Stay tuned this function will be updated soon. 
        '''
        

        logging.info(': Compound Match-Up function called!')

        # Pulls in our matched compound data
        #ask user for file input and read in the csv file. 
        file = filedialog.askopenfilename()
        my_data = pd.read_csv(file)

        
        if typeFile == 'all':

            my_final_data = np.zeros((len(my_data["Matched.Compound"]),2))
            my_final_data = pd.DataFrame(my_final_data,columns=['ID', 'Compound Name'])
            #grab the compound ID of interest
            lenCompounds = len(my_data['Matched.Compound'])
            for i in range(len(my_data["Matched.Compound"])):
                if (i+1)%100 ==0:
                    x = ((i+1)/lenCompounds)*100
                    x = float("{0:.2f}".format(x))
                    logging.info(': ' + str(x)+'%' + ' completed!')

                #input the values into a request from KEGG API
                if my_data['Matched.Compound'][i][0] == 'G':
                    my_final_data['ID'][i] = my_data['Matched.Compound'][i]
                    my_final_data['Compound Name'][i] = 'GAG subunits'

                elif my_data['Matched.Compound'][i][0] == 'C':
                    if my_data['Matched.Compound'][i][1] == 'E':
                        my_final_data['ID'][i] = my_data['Matched.Compound'][i]
                        my_final_data['Compound Name'][i] = 'Not in KEGG, will update soon!'

                    else:
                        try:
                            request = REST.kegg_get(my_data['Matched.Compound'][i])
                            txtFCur = my_data['Matched.Compound'][i] + '.txt'
                            open(txtFCur,'w').write(request.read())
                            records = Compound.parse(open(txtFCur))
                            record = list(records)[0]
                            os.remove(txtFCur)
                            my_final_data['ID'][i] = my_data['Matched.Compound'][i]
                            my_final_data['Compound Name'][i] = record.name

                        except:
                            logging.error(": No KEGG match found! Check KEGG Website!")
                            #messagebox.showerror(title="Error",message="ID not found in KEGG, as of 04.17.22, this is likely because it is a Glycan!!")

                            logString = my_data['Matched.Compound'][i]
                            logString = ': Failed to find-' + logString
                            logging.info(logString)
                            my_final_data['ID'][i] = my_data['Matched.Compound'][i]
                            my_final_data['Compound Name'][i] = 'No match in KEGG, investigate this compound further.'
                            continue

                else:
                    my_final_data['ID'][i] = my_data['Matched.Compound'][i]
                    my_final_data['Compound Name'][i] = 'Unknown'
            #save data frame as csv file for users
            my_final_data.to_csv(path_or_buf="CompoundMatchUps.csv", index=False)
            messagebox.showinfo(title="Success", message="Successfully generaterd CompoundMatchUps file!!")
            return

        elif typeFile == 'enrich':
            
            for i in range(len(my_data["Cpd.Hits"])):
                curString = my_data["Cpd.Hits"][i]
                #print(type(curString))
                curCpds = curString.split(';')
            
                #loop through each curCpds list to find matching compounds
                for j in range(len(curCpds)):
                    if curCpds[j][0] == 'G':
                        try:
                            curCpds[j] = 'GAG subunits'
                            my_data["Unnamed: 0"][i] = "GAG Metabolism"
                        except:
                            logging.error(': Failed to update the row Cpd.Hits!')
                            messagebox.showerror(title="Error",message="Unable to update excel sheet let Brady know and send him input spreadsheet")
                            return
                        
                    elif curCpds[j][0] =='C':

                        if curCpds[j][1] != 'E':
                            try:
                                curHit = REST.kegg_get(curCpds[j])
                            except:
                                logString = curCpds[j]
                                logString = ': Failed to find-' + logString
                                logging.info(logString)
                                curCpds[j] = curCpds
                                continue

                            try:
                                open('Compound.txt','w').write(curHit.read())

                            except:
                                logging.error(': Failed to open text! Let Brady know, this should rarely if ever occur!')
                                messagebox.showerror(title='Error',message='Failed to open text! Let Brady know, this should rarely if ever occur!')
                                return

                            records = Compound.parse(open('Compound.txt'))
                            record = list(records)[0]
                            curCpds[j] = record.name


                try:
                    my_data["Cpd.Hits"][i] = curCpds
                except:
                    logging.error(': Failed to update the row Cpd.Hits!')
                    messagebox.showerror(title="Error",message="Unable to update excel sheet let Brady know and send him input spreadsheet")
                    return

            #save the updated DataFrame as a EnrichmentIdentifications.csv
            my_data.to_csv(path_or_buf="EnrichmentIdentifications.csv",index=False)
            messagebox.showinfo(title="Sucess",message="EnrichmentIdentifications.csv has been been successfully created!!")

        return

    def compoundList(tol):
        '''
        Input a list of exact monoisotopic masses and determine the compounds associated with the masses

        Input:
        Excel sheet with exact masses. 

        Output:
        Updated excel sheet with Compound matches

        '''
        filename = filedialog.askopenfilename()
        glyList = pd.read_excel('C:/Users/Public/Documents/ClusteringGUI-develop/Glycans.xlsx')
        glyList = glyList.to_numpy()
        glyListMasses = glyList[:,1]
        glyListGID = glyList[:,0]



        #read in the data as a dataframe and convert to 
        data = GB.fileCheck(file=filename)
        compoundLookUp = data.to_numpy()

        tol = int(tol)/(10**6)
        

        compoundMatches = []
        for i in range(compoundLookUp.shape[0]):
            low = compoundLookUp[i] - (compoundLookUp[i]*tol)
            high = compoundLookUp[i] + (compoundLookUp[i]*tol)
            inputEM = str(low[0]) +'-'+ str(high[0])
            

            #input exact masses to kegg_find to determine the compound IDs on KEGG
            try:
                request = REST.kegg_find('compound',inputEM,"exact_mass")

            except:
                logging.error(': No matches found!')
                request = ['No matches']

            if isinstance(request,list):
                compoundMatches.append(request)
            else:
                
                try:
                    open('CompoundMatches.txt','w').write(request.read())

                except:
                    logging.error(': Failed to open text file! Let Brady know, this should rarely if ever happen!!')
                    messagebox.showerror(title='Error',message='Failed to open text file! Let Brady know, this should rarely if ever happen!')
                    return
                

                #open found compoundMatches file and put compound matches into a list. 
                lines = []
                with open('CompoundMatches.txt') as f:
                    line = f.readline()
                    while line:
                        line = f.readline()
                        lines.append(line)

                if len(lines) > 1:
                    for j in range(len(lines)):
                        if len(lines[j])>0:
                            #reformat string
                            curLine = lines[j].strip()
                            curLine = curLine.lstrip('cpd:')
                            curLine = curLine.split("\t")
                            lines[j] = curLine[0]

                        lookUpMasses = glyListGID[(np.where((glyListMasses >= low) & (glyListMasses <= high)))]
                    
                    if len(lookUpMasses) > 0:
                        lines.append(lookUpMasses)
                    compoundMatches.append(lines)

                else:
                    lookUpMasses = glyListGID[(np.where((glyListMasses >= low) & (glyListMasses <= high)))]
                    if len(lookUpMasses) == 0:
                        compoundMatches.append('No Matches')
                    else:
                        lookUpMasses = lookUpMasses.tolist()
                        compoundMatches.append(lookUpMasses)
                


    def ensembleClustering(optNum=2, minMetabs = 0, colorMap='viridis',linkParams=[],transform = 'None',scale='None', type='base'):
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
        file = filedialog.askopenfilename()
        data, col_groups = GB.readAndPreProcess(file=file,transform=transform,scale=scale,func='CC')
        del(col_groups)

        #determine whether data read in or not.
        if data is None:
            messagebox.showerror(title='Error',message='No file selected, returning to GUI. If you wish to continue with ensemble clustering, click continue and then select file!')
            return
        
        #read in data as dataframe for ease of use in recClusters, and ensembleClustersOut
        metab_data = GB.readAndPreProcess(file=file, transform='None', scale='None', func='Raw')
        #List for the use in creating and plotting the clustering results
        # linkParams = [['ward','euclidean'],['single','euclidean'],['single','sqeuclidean'],['single','seuclidean'],['single','chebyshev'],['complete','euclidean'],['complete','sqeuclidean'],['complete','seuclidean'],['complete','chebyshev'],['average','euclidean'],['average','sqeuclidean'],['average','seuclidean'],['average','chebyshev']]

        #calculate the number of clusterings based upon the size of the lists and an additional term for the ward-euclidean run. 
        numClusterings = (len(linkParams))

        #determine the the number of clusters and the dictionary location that needs to be called. 
        numMetabs = data.shape[0]
        dictLoc = numMetabs-optNum-1

        #create co-occurrence matrix.
        coOcc = GB.cooccurrence(data)

        for i in range(len(linkParams)):
            start = time.perf_counter()
            linkCur = linkage(data,linkParams[i][0],linkParams[i][1])
            valid = GB.clustConnectLink(linkCur)
            coOcc = GB.popCooccurrence(valid[dictLoc],coOcc,numClusterings)
            end = time.perf_counter()
            logging.info(': ' +str(linkParams[i][0])+'-'+str(linkParams[i][1]) +' done!')
            logging.info(str(end-start))
        del(linkParams)

        #make the coOccurence matrix a dataframe.
        coOcc1 = pd.DataFrame(coOcc)
        try:
            #try to save the large .csv file of the CoOccurence matrix.
            coOcc1.to_csv('EnsembleCoOcc.csv',index=False)
        except:
            logging.error(': Failed to save the Ensemble CoOccurence matrix!!')
            messagebox.showerror(title='Error',message='Unable to save Ensemble CoOccurent matrix, please inform Brady!')

        #create the ensemble dendrogram using ward-euclidean inputs. 
        GB.createEnsemDendrogram(coOcc,metab_data,norm=0,minMetabs=minMetabs,numClusts=numClusterings,link='ward',dist='euclidean',func="ensemble",colMap=colorMap)

        #Log the completgroupion of the ensemble clustering function
        logging.info(': Sucessfully completed Ensemble clustering!')
        
        return

    def MST(self,transform ='None',scale ='None', func = 'k-means based'):
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
        logging.info(': User called Cluster validation function.')

        filename = filedialog.askopenfilename()
        try:
            data, col_groups = GB.readAndPreProcess(file=filename, transform = transform, scale =scale, func='CC')
        except BaseException:
            logging.error(': Unable to proceed, due to file error!')
            messagebox.showerror(title='Error',message='Unable to proceed, try again or return to homepage!')
            return

        num_groups=data.shape[1]

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

        if func == 'k-means based':
            #Validate the number of clusters that should be used in the clustering solutions.
            #create a list of tuples containing the single cluster set, with the data and the num_groups
            argsMulti = []
            
            logging.info(": Starting k-means based cluster validation!")
            start = time.perf_counter()
            if len(validationClusters) < 100:
                for i in range(len(validationClusters)):
                    if i >= len(validationClusters)-(int(len(validationClusters)/2)):
                        argsMulti.append(({0:validationClusters[i]},data,num_groups))
            else:
                for i in range(len(validationClusters)):
                    if i >= len(validationClusters)-101:
                        argsMulti.append(({0:validationClusters[i]},data,num_groups))
            
            if __name__ == 'GUIUtils':
                with Pool(config.numThreads) as p:
                    valIndex = p.starmap(GB.Validate,argsMulti)
            
            end = time.perf_counter()
            logging.info(':'+str(end-start))
            valIndex = np.asarray(valIndex)
            GB.valPlotting(valIndex,mstOut)

        elif func=='DBI':
            logging.info(": Starting Davies-Bouldin cluster validation!")

            #start tracking the performance of non-threaded DBI validation
            start = time.perf_counter()
            valIndex = np.zeros((len(validationClusters),2))

            argsMulti = []
            if len(validationClusters) < 100:
                for i in range(len(validationClusters)):
                    if i >= len(validationClusters)-(int(len(validationClusters)/2)):
                        argsMulti.append(({0:validationClusters[i]},data,num_groups))
            else:
                for i in range(len(validationClusters)):
                    if i >= len(validationClusters)-101:
                        argsMulti.append(({0:validationClusters[i]},data,num_groups))

            if __name__ == 'GUIUtils':
                with Pool(config.numThreads) as p:
                    valIndex = p.starmap(VM.daviesBouldin,argsMulti)

            
            end = time.perf_counter()
            valIndex = np.asarray(valIndex)

            GB.valPlotting(valIndex,mstOut,valMet = func)

        elif func == 'Dunn':
            logging.info(": Starting Dunn cluster validation!")

            start = time.perf_counter()
            valIndex = np.zeros((len(validationClusters),2))

            argsMulti = []
            if len(validationClusters) < 100:
                for i in range(len(validationClusters)):
                    if i >= len(validationClusters)-(int(len(validationClusters)/2)):
                        argsMulti.append(({0:validationClusters[i]},data,num_groups))
            else:
                for i in range(len(validationClusters)):
                    if i >= len(validationClusters)-101:
                        argsMulti.append(({0:validationClusters[i]},data,num_groups))

            if __name__ == 'GUIUtils':
                with Pool(config.numThreads) as p:
                    valIndex = p.starmap(VM.dunnIndex,argsMulti)

            
            end = time.perf_counter()
            valIndex = np.asarray(valIndex)

            GB.valPlotting(valIndex,mstOut,valMet='Dunn')

        elif func == 'PBM':
            logging.info(": Starting PBM cluster validation!")

            #find the center of all the data.
            dataPatCenter = np.mean(data, axis=0)

            #put the center at the top of the numpy array to find the sum of the distances
            dataC = np.vstack([dataPatCenter,data])
            
            #find the distances and sum them
            distances = pdist(dataC)
            distances = squareform(distances)
            Eo = np.sum(distances[0,:])

            #start tracking the performance of the threaded PBM valdidation metric
            start = time.perf_counter()
            valIndex = np.zeros((len(validationClusters),2))

            argsMulti = []

            if len(validationClusters) < 100:
                for i in range(len(validationClusters)):
                    if i >= len(validationClusters)-(int(len(validationClusters)/2)):
                        argsMulti.append(({0:validationClusters[i]},data,num_groups,Eo))
            else:
                for i in range(len(validationClusters)):
                    if i >= len(validationClusters)-101:
                        argsMulti.append(({0:validationClusters[i]},data,num_groups,Eo))

            if __name__ == 'GUIUtils':
                with Pool(config.numThreads) as p:
                    valIndex = p.starmap(VM.PBM,argsMulti)

            
            end = time.perf_counter()
            valIndex = np.asarray(valIndex)
            GB.valPlotting(valIndex,mstOut,valMet='PBM')
            
        elif func == 'Silhouette':
            logging.info(": Starting Silhouette cluster validation!")

            #start tracking the performance of non-threaded DBI validation
            start = time.perf_counter()
            valIndex = np.zeros((len(validationClusters),2))

            argsMulti = []
            if len(validationClusters) < 100:
                for i in range(len(validationClusters)):
                    if i >= len(validationClusters)-(int(len(validationClusters)/2)):
                        argsMulti.append(({0:validationClusters[i]},data,num_groups))
            else:
                for i in range(len(validationClusters)):
                    if i >= len(validationClusters)-101:
                        argsMulti.append(({0:validationClusters[i]},data,num_groups))

            if __name__ == 'GUIUtils':
                with Pool(config.numThreads) as p:
                    valIndex = p.starmap(VM.Silhouette,argsMulti)

            end = time.perf_counter()
            valIndex = np.asarray(valIndex)
            GB.valPlotting(valIndex,mstOut,valMet="Silhouette")


    def peaksToPathways():
        '''
        Create input files for the mummichog algorithm, using files output from the ensemble clustering and the in the future from the clustergram function. 

        Input:

        peaksToPathways does not accept any inputs but will prompt the user for two inputs. First, the user will need to select the original files for there ensemble clustered data. 
        Next, the user will need to select the directory containing the files generated by the ensemble clustering function. 

        Output:

        This function will output csv files containing the m/z value and p-values fo the matched metabolites (p-values =0.04), and the remaining the metabolites with p-values equal to 1. 
        '''
        print(os.getcwd())
        logging.info(': Entering the Peaks to Pathways generator!')
        #ask user to input the file name of the original data
        messagebox.showinfo(title='File selection', message="Please select the original data file submitted for clustering!!")
        filename = filedialog.askopenfilename()

        dataRaw = GB.fileCheck(file = filename)
        if dataRaw is None:
            #log error and return function to ensure a soft closing of the class
            logging.error(': Error loading the reference Excel sheet.')
            return

        #ask user to select the directory containing the csv files from ensemble clustering output (currently only method available)
        messagebox.showinfo(title="Directory Selection",message="Please select the directory containing the ensemble clustering output files!")
        direct = filedialog.askdirectory()
        curDir = os.getcwd()

        #change the current working directory to 
        os.chdir(direct)



        files = glob.glob('*.xlsx')

        ensemFiles = []
        dirLog = os.getcwd()
        for i in range(len(files)):
            #strip the beginning of the strip off and then check the first 5 characters of the stripped string
            #curCheck = files[i].strip(nameCheckStriper)
            curCheck =files[i]
            if curCheck[0:5] == 'Ensem':
                #append to ensemFiles
                ensemFiles.append(direct + '/' + curCheck)
            elif curCheck[0:5] == 'Clust':
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
                dataCur = pd.read_excel(ensemFiles[i])
            except:
                logging.error(': Failed to read in the excel sheet, it is recommend to upload an excel workbook with a single sheet!')
                messagebox.showerror(title='Error', message="Failed to read the excel sheet, it is recommended to upload an excel workbook with a single sheet!")
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
                    logging.warning(': Creation of peaks to pathway files halted due to non-matching values, please make sure you have selected appropriate reference file.')
                    return
            #create the files that can be submitted to the csv saving file. 
            dataOut = np.zeros((dataRaw.shape[0],3))

            dataOut[:,0] = dataClust[:,0]
            dataOut[:,1] = dataClust[:,2]
            dataOut[:,2] = dataClust[:,1]

            dataOut = pd.DataFrame(dataOut,columns=["m.z","p.value",'r.t'])
            dataOut = dataOut.sort_values(by=["p.value"])

            p2pPre = 'PeaksToPathways'
            p2pSuf = '.csv'
            firstCheck = p2pPre + '01' + p2pSuf

            #create and/or navigate to P2PFiles folder to contain peaksToPathways output
            if os.path.isdir('P2PFiles'):
                os.chdir('P2PFiles')
            else:
                os.mkdir('P2PFiles')
                os.chdir('P2PFiles')

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

            os.chdir('..')
        logging.info(': Leaving the Peaks to Pathways Function!')
        os.chdir(curDir)
        messagebox.showinfo(title="Success",message="Success Peaks to Pathway files have been generated!!")
        return

    def PDFGenerator():
        '''
        Generates image of the selected clusters from the Cluster selection tool. 

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

    def heatmapAnalysis(linkFunc,distMet,cmap, transform = 'None', scale ='None'):
        '''
        Allows users to input a subset of the original clutergram from heatmap analysis. 

        Input:
        linkage function
        distance metric
        color map choice
        transform
        scale

        Output:

        Heatmap of the subset of metabolites given as input.

        '''
        try:
            file = filedialog.askopenfilename()
            data, col_groups = GB.readAndPreProcess(file=file,transform=transform,scale=scale, func="CC")
            if data is None:
                raise ValueError
        except:
            logging.error(': No file selected or issue with connecting to drive!')
            messagebox.showerror(title='Error loading Data',message='File was not selected or trouble connecting to drive')
            return

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
        groupCluster = np.transpose(data)

        g = sns.clustermap(data, method=linkFunc,metric=distMet, figsize=(7, 5), col_cluster=True,col_colors=col_groups,cmap=cmap,yticklabels=False,xticklabels=True)
        #g = sns.clustermap(data, figsize=(7, 5), yticklabels=False, xticklabels=True, row_cluster=False, col_linkage=groupLink, col_colors=col_groups, cmap=cmap, cbar_pos=(0.01, 0.8, 0.025, 0.175))
        plt.savefig('HeatMap.png',dpi=600,transparent=True)
        plt.show()
        


    def selectClusters(link,dist,transform = 'None', scale = 'None',cmap = 'viridis'):
        '''
        Function that pulls out the information from the plot and saves it until the user is ready to submit the clusters to the peaks to pathways function. 
        
        Input:
        linkage function
        distance metric
        transform
        scale
        color map

        Output:
        dendrogram allowing users to select clusters of interest

        '''

        #log that the user called the Create Clustergram function
        logging.info(': User called the Create Clustergram Function.')
        #check that the file the user selects is appropriate
        global metab_data
        file = filedialog.askopenfilename()
        metab_data = GB.readAndPreProcess(file=file,transform='None',scale='None',func='Raw')
        
        if metab_data is None:
            #log error message and return for soft exit.
            logging.error(': Error loading in the Excel sheet.')
            return  

        #get columns for the usage later in the creation of a dataframe to save for    
        columns = list(metab_data.columns)
        #read in data
        data = GB.readInColumns(metab_data)
        data_orig = metab_data.to_numpy()
        #Standardize the data before clustering the results
        logging.info(': Pre-processing the data.')

        #send the data off to the readAndPreProcess function for analysis. 
        data, col_groups = GB.readAndPreProcess(file=file,transform=transform,scale=scale,func="CC")
        del(col_groups)
        #create messagebox explaining to users how they need to select clusters.
        messagebox.showinfo(title='Cluster Selection Info.', message='Select clusters of interest, cluster and peak to pathway files will be automatically generated!')

        #Create the appropriate plt figure to allow for the comparison of linkage functions
        fig, axes = plt.subplots(1,1,figsize=(8,8))

        #find the linkages
        linkageOne = linkage(data,link,metric=dist)

        if len(linkageOne[:,2]) == len(np.unique(linkageOne[:,2])):
            logging.info('No need to jitter data!')

        else:
            logging.info(': Matching distance need to jitter distances')
            values, counts = np.unique(linkageOne[:,2],return_counts=True)

            #get the locations where the counts are greater than 1 (i.e., the distances are matching)
            matchingDists = np.where(counts>1)

            for j in range(len(matchingDists[0])):
                #find the location of the values which have matching distances
                curLinkListLoc = np.where(linkageOne[:,2]==values[matchingDists[0][j]])
                curLinkListLoc = curLinkListLoc[0]
                for k in range(len(curLinkListLoc)):
                    if k > 0:
                        linkageOne[curLinkListLoc[k],2] += (k*0.000001)+0.000001

        groupCluster = np.transpose(data)
        linkageG = linkage(groupCluster,link,metric=dist)
        #create the dendrogram
        dend = dendrogram(linkageOne,ax=axes,above_threshold_color='y',orientation='left',no_labels=True)
        dendG = dendrogram(linkageG,ax=axes,above_threshold_color='y',orientation='left',no_labels=True)
        #Rework the data to create the clustergram
        metaboliteDendLeaves = dend['leaves']
        #find the maximum leaf to know what the index must be larger than for filling in the color
        maxLeaf = np.array(metaboliteDendLeaves)
        maxLeaf = np.amax(maxLeaf)
        groupDendLeaves = dendG['leaves']
        plt.close()
        fig, axes = plt.subplots(1,1,figsize=(8,8))
        dataFinal = np.zeros((data.shape[0],data.shape[1]))

        for i in range(data.shape[1]):
            #rearranging the data for heatmap
            for j in range(data.shape[0]):
                #going through the metabolites
                dataFinal[j,i] = data[metaboliteDendLeaves[j],groupDendLeaves[i]]

        columns.pop(0)
        columns.pop(len(columns)-1)
        columnsNew = []

        for i in range(len(groupDendLeaves)):
            columnsNew.append(columns[groupDendLeaves[i]])

        dataFinalDF = pd.DataFrame(dataFinal,columns=columnsNew)
        dataFinalDF.to_excel('Heatmap.xlsx',index=False)
        #create the axes in which the heatmap will be mapped upon
        plt.cla()
        heatmapAxes = [0.3, 0, 0.68, 1]
        heatmapAxes = fig.add_axes(heatmapAxes)
        heatmapAxes.matshow(dataFinal,aspect ='auto',origin='upper',cmap= cmap)
        
        maxList = np.zeros((1,2))
        maxList[0,0] = linkageOne[len(linkageOne)-1][0]
        maxList[0,1] = linkageOne[len(linkageOne)-1][1]
        maxLinkNum = int(np.amax(maxList)+2)

        metabDendAxes =[0,0, 0.3, 1]
        metabAxes = fig.add_axes(metabDendAxes)
        plt.cla()

        for i in range(len(dend['icoord'])):
            #plot each of the linkages one by one. 
            x = np.array(dend['icoord'][i])
            y = np.array(dend['dcoord'][i])
        
            plt.plot(-y, -x,'k')
            plt.draw()
        right, left = plt.xlim()
        plt.xlim(right,0)
        #find the length of data rows, to adjust the axes to fit the current heatmap 
        lenRows = data.shape[0]
        right = -(10*lenRows)
        plt.ylim(right,0)

        #sending linkageOne to the function which will point the selection to the appropriate
        #number of linkages to color.
        linkDir = GB.linkDir(linkageOne,maxLeaf)
        linkageClusters = GB.clustConnectLink(linkageOne)

        colSel = 0
        open('ClustColor.txt','w').write(str(colSel))
        open('ClusterReference.txt','w').write(str(time.strftime("%a_%b_%d_%Y_%H_%M_%S")))
        open('ClusterReference.txt','a').write("\n"+str(len(dataFinalDF[columnsNew[0]])))
        #create an interactive cursor
        cursor = mplcursors.cursor(multiple=True)
        cursor.visible =False
        cursor.connect("add", lambda sel: GB.select(sel.target,dend,linkageOne,linkDir,linkageClusters,data_orig))
        plt.show()


    def localWeighted():
        '''
        This function performs locally-weighted ensemble clustering

        Input: TBD

        Output:
        Recommended clusters
        Ensemble Clustergram

        '''

        #optimum number of clusters from validation index.
        sys.setrecursionlimit(10**8) 
        metab_data = GB.fileCheck()

        if metab_data is None:
            logging.error(': File does not meet input requirements.')   
            return
       
        #List for the use in creating and plotting the clustering results
        linkageList = ['single','complete','average']
        distList = ['euclidean','sqeuclidean','chebyshev','seuclidean'] 
        
        #calculate the number of clusterings based upon the size of the lists and an additional term for the ward-euclidean run. 
        numClusterings = (len(linkageList)*len(distList))+1

        #read in the data
        data = GB.readInColumns(metab_data)
        
        #Standardize the data before clustering the results
        logging.info(': Pre-processing data.')
        for i in range(data.shape[0]):
            data[i,:] = GB.standardize(data[i,:])

        #creates empty dictionary for clusterings
        clusters = {}

        #performs first 12 base clusterings and populates clusters dictionary
        for i in range(len(linkageList)):
            for j in range(len(distList)):
                linkCur = linkage(data,linkageList[i],distList[j])
                valid = GB.clustConnectLink(linkCur)
                index = str(linkageList[i] + '_' + distList[j])
                clusters.update({index:valid})
                logging.info(str(linkageList[i])+'-'+str(distList[j]) +' done!')

        #performs 13th base clustering and populates clusters dictionary
        linkCur = linkage(data, 'ward', 'euclidean')
        valid = GB.clustConnectLink(linkCur)
        logging.info(str('ward-euclidean done!'))
        clusters.update({'ward_euclidean':valid})


        optNum = 2
        
        refClust = LW.clustFinder(data=data,optNum=optNum,clusters=clusters)
        
        
        ECI = LW.clustCompare(refClust)
       
        consensusMat = LW.consensus(ECI,refClust,data)
        regionsOut = LW.regions(consensusMat)

    def enzymeLookUp(numSheets):
        '''
        '''

        #have the user select the file they would like to have read in.
        filename = filedialog.askopenfilename()
        
        #heatmapEnzyme Outputs
        outFile = 'HeatmapEnzyme.xlsx'
        writer = pd.ExcelWriter(outFile, engine='xlsxwriter')
        for i in range(numSheets):
            #read in each sheet
            dataCur = pd.read_excel(filename,sheet_name=i)

            # except:
            #     messagebox.showerror(title="Error opening file", message="Select file to continue!")
            #     logging.error(': No file selected sending back to GUI!')
            #     return

            dFDict = {}
            print("Starting sheet number: " + str(i+1))
            #get the compounds and determine how many pathways hits there are.
            for j in range(len(dataCur['Cpd.Hits'])):
                #for the current compound hits find the length
                CpdList = dataCur['Cpd.Hits'][j]

                #split the current list by ;
                CpdList = CpdList.split(';')
                
                for k in range(len(CpdList)):
                    #get the enzyme numbers from KEGG
                    try:
                        request = REST.kegg_get(CpdList[k])

                    except:
                        messagebox.showerror(title="Error",message="Cannot find compound, this should not happen")
                        logging.error(': Compound not found this should not happen!')

                    txtFCur = CpdList[k] + '.txt'
                    open(txtFCur,'w').write(request.read())
                    
                    records = Compound.parse(open(txtFCur))
                    #os.remove(txtFCur)

                    #get the record of the compound currently being looked up.
                    try:
                        record = list(records)[0]
                    except:
                        logging.error(': Almost for sure a glycan was found.')

                    if k == 0 and j == 0:
                        dict = {'Cluster #':[i+1],
                                'Pathway':dataCur['Pathway'][0],
                                'Pathway Total':dataCur['Pathway total'][0],
                                'Hits.total':dataCur['Hits.total'][0],
                                'Hits.sig':dataCur['Hits.sig'][0],
                                'Gamma':dataCur['Gamma'][0],
                                'Cpd.Hits':CpdList[0],
                                'Enzyme #s':[0]
                                }
                        #create a spreadsheet for the current hits
                        dFDict[i] = pd.DataFrame(dict)
                        dFDict[i]['Enzyme #s'][0] = record.enzyme
                    
                    elif k == 0 and j !=0:
                        #create a spreadsheet for the current hits
                        dFDict[i].loc[len(dFDict[i].index)] = [i+1,dataCur['Pathway'][j],dataCur['Pathway total'][j],dataCur['Hits.total'][j],dataCur['Hits.sig'][j],dataCur['Gamma'][j],CpdList[0], record.enzyme]

                    else:
                        dFDict[i].loc[len(dFDict[i].index)] = [None,None,None,None,None,None,CpdList[k],record.enzyme]
        
            dFDict[i].to_excel(writer,sheet_name=str(i+1))
        try:
            writer.save()

        except:
            messagebox.showerror(title='No worky',message='Need to investigate further')


        messagebox.showinfo(title="Success", message="Successfully completed getting enzyme IDs for each compound!")

    def anovaHM(transform ='Log transformation',scale='Auto Scaling',cMap = 'viridis'):
        '''
        Allows users to plot the top ### of data objects from ANOVA analysis

        Input:
        transform
        scale
        color map

        Output:
        heatmap figure
        Raw data for top features
        '''
        func = 'CC'
        #ask the user for the input excel workbook needs to contain two sheets
        file = filedialog.askopenfilename()

        #open sheet 0 - containing the original data
        #open sheet 1 - containing all ANOVA outcomes or the pre-truncated ANOVA results
        try:
            anovaRes = pd.read_excel(file,sheet_name=1)
        except:
            logging.error(': Unable to open file, may due to no input file')
            messagebox.showerror(title="Error", message='Unable to open file!')
            return



        dataOrig = pd.read_excel(file,sheet_name=0)
        dataOrig = dataOrig.iloc[1:,:]
        colHeaders = list(dataOrig.columns)
        colHeaders.pop(0)
        colHeaders.pop(len(colHeaders)-1)
        anovaResC = list(anovaRes.columns)


        #get column names out
        metabNames = list(dataOrig.columns)
        metabNames = dataOrig[metabNames[0]]
        del(dataOrig)
        
        #reads in all column headers and trims off first and last columns
        dataOrig = GB.readAndPreProcess(file=file,transform=transform,scale=scale,func='ANHM')

        dataUpdated = dataOrig[metabNames.isin(anovaRes[anovaResC[0]])]
        del(dataOrig)

        #send data to be transformed

        dataUpdated = pd.DataFrame(dataUpdated,index = anovaRes[anovaResC[0]],columns=colHeaders)
        dataUpdated.to_excel('RawTopANOVA.xlsx')
        print(cMap)
        g=sns.clustermap(dataUpdated,cMap=cMap)
        plt.show()


        return

    def confidenceIntervals(sampSize, confidenceLevel = 95):
        '''
        '''
        messagebox.showinfo(title="Input Order",message="Please first select a column oriented metabolomics files with the intensities between m/z (do not label this column) and rtmed columns. Next, select t_test.csv file from Metaboanalyst output.")
        #get file name and read in the original data file
        dataOrigLoc = filedialog.askopenfilename()
        dataOrig = pd.read_excel(dataOrigLoc)

        #put the original mz values into a list
        mzOrig = dataOrig['Unnamed: 0'].tolist()

        #get file name and read in the tTests for p-values <0.1
        tTestsLoc = filedialog.askopenfilename()
        tTests = pd.read_csv(tTestsLoc)


        #put the tTestMzs into a list for searching
        tTestMzs = tTests['Unnamed: 0'].tolist()

        #have the user input the sample size and confidence level
        sampleSize = sampSize
        confLevel = confidenceLevel
        sampleSize = float(sampleSize)
        confLevel = float(confLevel)

        #calculate the df and quantile, given the user inputs
        df = (sampleSize*2)-2
        q = 1- (1- ((confLevel)/100))/2

        tCritical = t.ppf(q,df)

        #search the values in the list against the original values to find the values of interest.
        CIs = []
        for i in range(len(tTestMzs)):
            #determine if the current list item has one or two decimals.
            if tTestMzs[i].count('.') > 1:
                #enumerate the original mz values each time
                mzNums = enumerate(mzOrig)
                #find the locations in  the list containing the wanted m/z value

                curList = [k for k, j in mzNums if j -float(tTestMzs[i][:tTestMzs[i].rfind('.')])<=0.0000001]

                #convert to numpy array... and remove m/z and rtmed. 
                ser = dataOrig.iloc[curList[int(tTestMzs[i][tTestMzs[i].rfind('.')+1:])]]
                ser = ser.to_numpy()
                ser = ser[1:-1]
                #add small "jitter" so that log transform does not fail
                ser = ser + 0.0000001
                
                #for now log-transform
                ser = np.log10(ser)

                #from the sample size pick out the two groups
                g1 = ser[0:int(sampleSize)]
                g2 = ser[int(sampleSize):]

                ciUpper = (stat.mean(g1)-stat.mean(g2)) + (tCritical* (((stat.variance(g1)/sampleSize) + (stat.variance(g2)/sampleSize))**.5))
                ciLower = (stat.mean(g1)-stat.mean(g2)) - (tCritical* (((stat.variance(g1)/sampleSize) + (stat.variance(g2)/sampleSize))**.5))

                ciUpper = 10**ciUpper
                ciLower = 10**ciLower
                CIs.append((round(ciLower,2),round(ciUpper,2)))

            else:
                #enumerate the original mz values each time
                mzNums = enumerate(mzOrig)
                curList = [k for k, j in mzNums if float(tTestMzs[i])-j <= 0.0000001]

                ser = dataOrig.iloc[curList[0]]
                ser = ser.to_numpy()
                ser = ser[1:-1]
                #add small "jitter" so that log transform does not fail
                ser = ser + 0.0000001
                
                #for now log-transform
                ser = np.log10(ser)

                #from the sample size pick out the two groups
                g1 = ser[0:int(sampleSize)]
                g2 = ser[int(sampleSize):]

                ciUpper = (stat.mean(g1)-stat.mean(g2)) + (tCritical* (((stat.variance(g1)/sampleSize) + (stat.variance(g2)/sampleSize))**.5))
                ciLower = (stat.mean(g1)-stat.mean(g2)) - (tCritical* (((stat.variance(g1)/sampleSize) + (stat.variance(g2)/sampleSize))**.5))
                ciUpper = 10**ciUpper
                ciLower = 10**ciLower
                CIs.append((round(ciLower,2),round(ciUpper,2)))

        tTests['CIs'] = CIs

        tTests.to_excel('t_testWCIs.xlsx',index=False)
        messagebox.showinfo(title="Success",message="A t_testWCIs.xlsx file has successfully been created!")
        return

    def bootstrapping(numReSamp,numPerSamp):
        '''
        Hello, there!!
        '''

        #log that user called MST
        logging.info(': User called the bootstrapping function.')

        #get the file of interest
        filename = filedialog.askopenfilename()

        try:
            data = GB.readAndPreProcess(file=filename, func='else')
        
        except BaseException:
            logging.error(': Unable to proceed, due to file error!')
            messagebox.showerror(title='Error',message='Unable to proceed, try again or return to homepage!')
            return


        #determine the number of groups in the sample set.
        num_groups=data.shape[1]
        num_metabs = data.shape[0]


        if num_groups >= int(numPerSamp):
            bootCur = []
            for i in range(int(numReSamp)):
                bootCur.append(stat.mean(random.choices(data[1,:], k=int(numReSamp))))

        #convert the current list to a numpy array for ease of using the sort function
        bootCur = np.array(bootCur)

        #sort the bootCur numpy array in descending order
        bootCur = np.sort(bootCur)

        #calculate LB and UB indicies for analysis
        indLB = int(numReSamp)*0.025
        indUB = int(numReSamp)*0.975

        #linear interpolation of indLB and indUB to land at 95%CI
        if math.floor(indLB) != math.ceil(indLB):
            #calculate linear interpolated value
            bootCILB = bootCur[int(math.floor(indLB))] + ((indLB-math.floor(indLB))(bootCur[int(math.ceil(indLB))] - bootCur[int(math.floor(indLB))]))
        else:
            bootCILB = bootCur[int(math.floor(indLB))]

        if math.floor(indUB) != math.ceil(indUB):
            #calculate linear interpolated value
            bootCIUB = bootCur[int(math.floor(indUB))] + ((indUB-math.floor(indUB))(bootCur[int(math.ceil(indUB))] - bootCur[int(math.floor(indUB))]))

        else:
            bootCIUB = bootCur[int(math.ceil(indUB))]

        bootCI = (bootCILB,bootCIUB)
        

        g = sns.kdeplot(bootCur,color='r',fill=True)
        #plot vertical lines of upper and lower bounds
        g.vlines(bootCI,0,g.get_ylim()[1],colors=['k'])

        plt.show()



    def normalityCheck(transform=config.curTrans,scale=config.curScale):
        '''
        '''

        #have user input the wanted file
        file =filedialog.askopenfilename()

        #read in the metabolites file
        print(transform,scale)
        data = GB.readAndPreProcess(file=file,transform=transform,scale=scale)
        dataOrig = GB.readAndPreProcess(file)

        #reformat both datasets to incorporate the appropriate data format for kdeplots
        data = data.reshape(int(data.shape[0]*data.shape[1]),1)
        data = pd.DataFrame(data,columns=['mz'])
        dataOrig = dataOrig.reshape(int(dataOrig.shape[0]*dataOrig.shape[1]),1)
        dataOrig = pd.DataFrame(dataOrig,columns=['mz'])

        fig, axes = plt.subplots(1, 2, figsize=(10,5))
        go = sns.kdeplot(data=dataOrig,x='mz',ax=axes[0])
        axes[0].set_title("Raw Data")
        g = sns.kdeplot(data=data,x='mz',ax=axes[1])
        axes[1].set_title("Transform and/or Scaled")
        plt.savefig("Normalized.png",dpi=600,transparent=True)
        plt.show()


