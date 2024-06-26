#Creating a class containing functions that will be used in GUI
from numpy.lib.arraysetops import isin
import pandas as pd
import numpy as np
import time
import platform
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.chrome.options import Options
from matplotlib import pyplot as plt
import matplotlib as mpl
from scipy.cluster.hierarchy import dendrogram
from scipy.cluster.hierarchy import linkage, ward, fcluster
from scipy.spatial.distance import pdist, squareform
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import minimum_spanning_tree
from scipy.stats import t
from sklearn.cluster import AgglomerativeClustering as AC
from sklearn import metrics

import glob,sys,logging,time,getpass,fpdf,os
import statistics as stat
from multiprocessing import Pool
import seaborn as sns
import config
import random
import math
from itertools import combinations
import zipfile
import os

import GuiBackground as GB
from tkinter import filedialog, messagebox
import mplcursors

from Bio.KEGG import REST
from Bio.KEGG import Compound

from ValidationMetric import ValidationMetric as VM

class GUIUtils:
    def dataIntegrity(file):
        '''
        The data integrity allows users to determine whether the first column of data in a dataframe contains double decimals. This error mostly arises when working with MetaboAnalyst results. 

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
        output.close()
        logOut = 'Updated file saved as: ' + file
        logging.info(logOut)
        #log that the data integrity function has been sucessfully completed. 
        logging.info(': Data Integrity check sucessfully completed.')
        messagebox.showinfo(title="Success",message="Removed data integrity issues!!")
        return

    def createClustergram(norm,linkFunc,distMet,cmap,colOrder=[], transform = 'None', scale ='None'):
        '''
        The function is responsible for generating the clustergrams for multivariate data. This function is capable of
        using all the linkage functions and distance measures currently implemented to the scipy.hierarchy method. Note the ward-euclidean distance is the only combination 
        available from the scipy.hierarcy package. 
        
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

        fileI = filedialog.askopenfilename()
        
        try:
            if norm == 0 or norm == 2:
                
                data, col_groups = GB.readAndPreProcess(file=fileI,transform=transform,scale=scale,func="CC")

            elif norm == 1:
                
                data, col_groups = GB.readAndPreProcess(file=fileI,transform=transform,scale=scale,func="CC",first=colOrder[0]) 
        except TypeError:
            logging.error(': No file selected!')
            messagebox.showerror(title='Error',message='Error loading in the data, normilizing it or something else. Contact Brady if cannot figure out.')
            return

        #create dendrogram and plot data        
        GB.create_dendrogram(data,col_groups, norm, link=linkFunc, dist=distMet,color = cmap,colOrder=colOrder)

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

        #getting the needed data. 
        df = pd.read_excel(file)
        rt = df[list(df.columns)[-1]]
        mz = df[list(df.columns)[0]]

        #transpose the matrix to get what I want.
        df = df.T
        out = df.groupby([0]).median()
        out = out.T

        out.insert(0,"mz",mz)
        out.insert(out.shape[1],"rtmed",rt)

        #set up a place holder with first and last places as None
        ph = [i for i in range(0,out.shape[1])]
        ph[0] = None;ph[-1]=None
        out.loc[0,list(out.columns)] = ph

        #reorder the indicies and save the file with _medians attached
        out = out.sort_index(ascending=True)
        outFile = os.path.basename(file)
        outFile = outFile[:-5]
        outFile += "_medians.xlsx"
        out.to_excel(outFile,index=False)

        #logging the completion of the group medians function
        logging.info(': Successfully calculated the medians of each metabolite in each group!')
        messagebox.showinfo(title="Success",message="Successfully created medians file!!")
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

        data, col_groups = GB.readAndPreProcess(file=file, transform=transform, scale=scale, func="CC")
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
            messagebox.showinfo(title="Success",message="EnrichmentIdentifications.csv has been been successfully created!!")

        return

    def compoundList(tol):
        '''
        Input a list of exact monoisotopic masses and determine the compounds associated with the masses

        Input:
        Excel sheet with exact masses. 

        Output:
        Updated excel sheet with Compound matches

        '''
        logging.info(': Starting exact mass comparison against KEGG compounds, and glycans!')
        filename = filedialog.askopenfilename()
        data = pd.read_excel(filename)
        lookUpList = pd.read_excel(config.keggGlycanLoc)

        #create a numpy array of just the masses
        lookUpMasses = lookUpList['Masses']
        lookUpMasses=lookUpMasses.to_numpy()
        lookUpMasses = lookUpMasses.reshape(np.shape(lookUpMasses)[0],1)

        #convert the look-up data to a numpy array
        cmdLU = data.to_numpy()

        #allowable tolerance as an vector of tolerances
        tol = (10/(10**6))*np.ones((np.shape(cmdLU)[0],1))

        #calculate the tolerances for each compound
        cmdTols = np.matmul(np.diagflat(cmdLU), tol)

        #creating matricies of the appropriate size for comparison with lookUpMasses List
        M_plus1     = np.matmul(np.ones((np.shape(lookUpMasses)[0],1)),np.transpose(cmdLU))
        M_plusNa    = np.matmul(np.ones((np.shape(lookUpMasses)[0],1)),np.transpose(cmdLU)) - (22.989769*np.ones((np.shape(lookUpMasses)[0],np.shape(cmdLU)[0])))
        M_plusH     = np.matmul(np.ones((np.shape(lookUpMasses)[0],1)),np.transpose(cmdLU)) - (1.00784*np.ones((np.shape(lookUpMasses)[0],np.shape(cmdLU)[0])))
        M_plusNaH   = np.matmul(np.ones((np.shape(lookUpMasses)[0],1)),np.transpose(cmdLU)) - ((22.989769+1.00784)*np.ones((np.shape(lookUpMasses)[0],np.shape(cmdLU)[0])))

        #Create a look up table for given number of compounds, and glycans
        lookUpMat = np.matmul(lookUpMasses,np.ones((1,np.shape(cmdLU)[0])))

        #compare all four compound Look Ups
        M_plus1_c = M_plus1 - lookUpMat
        M_plusNa_c = M_plusNa - lookUpMat
        M_plusH_c = M_plusH - lookUpMat
        M_plusNaH_c = M_plusNaH - lookUpMat

        #loop over the matricies to determine the most likely compound matches
        cID = []
        cIDM = []
        adduct = []
        for i in range(np.shape(M_plus1_c)[1]):
            #get current checking value for each array 
            curCheck_mp1  = np.min(abs(M_plus1_c[:,i]))
            curCheck_mpNa = np.min(abs(M_plusNa_c[:,i]))
            curCheck_mpH = np.min(abs(M_plusH_c[:,i]))
            curCheck_mpNaH = np.min(abs(M_plusNaH_c[:,i]))
            
            #compare the minimums of each to the appropriate tolerance
            curCheck_mp1_tf = curCheck_mp1 < cmdTols[i]
            curCheck_mpNa_tf = curCheck_mpNa < cmdTols[i]
            curCheck_mpH_tf = curCheck_mpH < cmdTols[i]
            curCheck_mpNaH_tf = curCheck_mpNaH < cmdTols[i]
            
            #number of potential matches given ppm
            tol_check = [curCheck_mp1_tf, curCheck_mpNa_tf, curCheck_mpH_tf, curCheck_mpNaH_tf]
            minDiffCheck = [curCheck_mp1, curCheck_mpNa, curCheck_mpH, curCheck_mpNaH]
            possibleCount = tol_check.count(True)
            
            compoundList = []
            if possibleCount > 0:
                #append all matches
                for j in range(len(tol_check)):
                    if tol_check[j] == True:
                        compoundList.append(minDiffCheck[j])
                    else:
                        compoundList.append(1000)
                        
                #calculate and location minimum of the list
                indexOfInterest = compoundList.index(min(compoundList))
                
                if indexOfInterest == 0:
                    cID.append(lookUpList['KEGG ID'][np.where(abs(M_plus1_c[:,i])==compoundList[indexOfInterest])[0][0]])
                    cIDM.append(lookUpList['Masses'][np.where(abs(M_plus1_c[:,i])==compoundList[indexOfInterest])[0][0]])
                    adduct.append('M[1+]')
                elif indexOfInterest == 1:
                    cID.append(lookUpList['KEGG ID'][np.where(abs(M_plusNa_c[:,i])==compoundList[indexOfInterest])[0][0]])
                    cIDM.append(lookUpList['Masses'][np.where(abs(M_plusNa_c[:,i])==compoundList[indexOfInterest])[0][0]])
                    adduct.append('M+Na[1+]')
                elif indexOfInterest == 2:
                    cID.append(lookUpList['KEGG ID'][np.where(abs(M_plusH_c[:,i])==compoundList[indexOfInterest])[0][0]])
                    cIDM.append(lookUpList['Masses'][np.where(abs(M_plusH_c[:,i])==compoundList[indexOfInterest])[0][0]])
                    adduct.append('M+H[1+]')
                elif indexOfInterest == 3:
                    cID.append(lookUpList['KEGG ID'][np.where(abs(M_plusNaH_c[:,i])==compoundList[indexOfInterest])[0][0]])
                    cIDM.append(lookUpList['Masses'][np.where(abs(M_plusNaH_c[:,i])==compoundList[indexOfInterest])[0][0]])
                    adduct.append('M+Na+H[1+]')
                
            else:
                cID.append('No matches')
                cIDM.append(0)
                adduct.append('No adduct found')

        dataCID = data.assign(CompoundID = cID)
        dataFinal = dataCID.assign(CompoundMass = cIDM)
        dataFinal = dataFinal.assign(Adduct = adduct)

        dataFinal.to_excel('CompoundMatches.xlsx',index=False)
        logging.info(': Completed!')
        messagebox.showinfo(title="Success", message="Compound Matches has been generated!")
        return
    
    def ensembleClusteringFullOpt(parameters, numClusts = 10, transform='None',scale='None'):
        '''
        '''

        #log that the user called ensemble clustering function
        logging.info(': User called Ensemble Clustering function.')

        #optimum number of clusters from validation index.
        sys.setrecursionlimit(10**8)
        file = filedialog.askopenfilename()
        data, col_groups = GB.readAndPreProcess(file=file,transform=transform,scale=scale,func='CC')
        metab_data = GB.readAndPreProcess(file=file,transform='None',scale='None',func='Raw')
        del(col_groups)

        #determine whether data read in or not.
        if data is None:
            messagebox.showerror(title='Error',message='No file selected, returning to GUI. If you wish to continue with ensemble clustering, click continue and then select file!')
            return
        
        # create labels
        best_labels = [None]*parameters.shape[0]

        #setting the functions into a validation, and distance metrics. 
        valIndex = {
            'CH':GB.calinskiHarabasz_correlation,
            'CH_':metrics.calinski_harabasz_score,
            'SIL':metrics.silhouette_score,
            'DBI':GB.daviesBouldinScore_correlation,
            'DBI_':metrics.davies_bouldin_score
        }

        distance = {
            'CNS': GB.correlationNosqrt,
            'CS': GB.correlationSqrt,
            'PW': GB.pairWise
        }


        #create co-occurrence matrix.
        coOcc = GB.cooccurrence(data)
        for i in range(parameters.shape[0]):
            bestScore = 0
            optClust = [None]*2
            #calculate the distance matrix
            if parameters['Linkage'][i] == 'ward':
                dist = distance[parameters['Distance'][i]](data,metric=parameters['Correlation'][i])
                dist = squareform(dist)
                link_mat = ward(dist)
                for j in range(numClusts):
                    labels_ = fcluster(link_mat,j+2,criterion='maxclust')
                    #update the best clustering solutions
                    score = valIndex[parameters['Optimizer'][i]](data,labels_,parameters['Distance'][i])
                    if score > bestScore:
                        optClust[0],optClust[1] = j+2, labels_-1
                        bestScore = score
                best_labels[i]= optClust[1]    
                optClusters = dict.fromkeys(list(range(0,optClust[0])),[])
                for k in optClusters:
                    optClusters.update({k:np.where(optClust[1]==k)[0].tolist()})
                #update the co-occurrence matrix
                coOcc = GB.popCooccurrence(optClusters,coOcc,parameters.shape[0])

            elif parameters['Optimizer'][i] =="DBI":
                #update such that best score is lowest
                bestScore =10**10

                #get the distance metric out
                dist = distance[parameters['Distance'][i]](data,metric=parameters['Correlation'][i])
                for j in range(numClusts):
                    #calculate the clustering solutions
                    agglo = AC(n_clusters=j+2,linkage=parameters['Linkage'][i],
                            metric='precomputed').fit(dist)
                    
                    #update the best clustering solutions
                    score = valIndex[parameters['Optimizer'][i]](data,agglo.labels_+1,parameters['Distance'][i])
                    if score < bestScore:
                        optClust[0],optClust[1] = j+2, agglo.labels_
                        bestScore = score
                best_labels[i]= optClust[1]    
                optClusters = dict.fromkeys(list(range(0,optClust[0])),[])
                for k in optClusters:
                    optClusters.update({k:np.where(optClust[1]==k)[0].tolist()})
                #update the co-occurrence matrix
                coOcc = GB.popCooccurrence(optClusters,coOcc,parameters.shape[0])

            else:
                #get the distance metric for clustering
                dist = distance[parameters['Distance'][i]](data,metric=parameters['Correlation'][i])      
                for j in range(numClusts):
                    #calculate the clustering solutions
                    agglo = AC(n_clusters=j+2,linkage=parameters['Linkage'][i],
                            metric='precomputed').fit(dist)
                    
                    #update the best clustering solutions
                    score = valIndex[parameters['Optimizer'][i]](dist,agglo.labels_,metric='precomputed')
                    if score > bestScore:
                        optClust[0],optClust[1] = j+2, agglo.labels_
                        bestScore = score
                best_labels[i]= optClust[1]    
                optClusters = dict.fromkeys(list(range(0,optClust[0])),[])
                for k in optClusters:
                    optClusters.update({k:np.where(optClust[1]==k)[0].tolist()})
                
                coOcc = GB.popCooccurrence(optClusters,coOcc,parameters.shape[0])

        #make the coOccurence matrix a dataframe.
        coOcc1 = pd.DataFrame(coOcc)
        try:
            #try to save the large .csv file of the CoOccurence matrix.
            coOcc1.to_excel('EnsembleCoOcc.xlsx',index=False)
        except:
            logging.info('Issue saving the co-occurence matrix as .csv. This error should not occur, should have been caught earlier')

        #generate linkage function
        dissim = 1 - np.around(coOcc,decimals=3)
        dissim = squareform(dissim)
        linkageMetabOut = linkage(dissim,'average')

        #Get the labels for the optimal solution and plot the comparisons 
        labelsCoOcc =fcluster(linkageMetabOut,2,'maxclust')


        #get the labels in a tuple to send to the function    
        ensemble_s = parameters[['Distance','Linkage']]
        distLink = tuple(ensemble_s.itertuples(index=False,name=None))
        GB.randComp(best_labels,distLink)
        GB.adjRandComp(best_labels,distLink)
        GB.coOccMonoComp(best_labels,labelsCoOcc,distLink)

        #create the ensemble dendrogram using ward-euclidean inputs. 
        GB.createEnsemDendrogramNew(coOcc,metab_data,data,norm=0,minMetabs=0,numClusts=len(parameters),
                                link='average',dist='euclidean',func="ensemble",colMap='viridis');
    
        return      

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
    
    def validateMono(self, upperLimClust=10, transform ='None',scale ='None',linkage='average',dissimilarity='euclidean', func = 'SIL'):
        '''
        This functionality is designed to perform clustering optimization on a set of input data. Specifically this is designed to 
        generate agglomerative hierarchical clustering solutions to be optimized using a optimization metric.

        Input:
        self - the object for the UI
        
        upperLimClust - Upper limit of the number of clusters the user would like to consider for optimization purposes

        transform - which transform would the user like to use on their data. No transform is the default

        scale - which data scaling would the user like to use on their data. No scaling is the default.

        linkage - which agglomerative hierarchical clustering linkage function would the user like to use. Average linkage function is the default

        dissimilarity - which dissimilarity measure would the user like to use. Euclidean dissimilarity is the default. 

        func - which optimization metric would the user like to use. The default is 'SIL' (i.e., Silhouette)

        Output: 
        Matplotlib plot of the optimization

        .csv of the optimization metric scores for the considered metric

        Recommended number of clusters given the optimization metric. 

        '''
        


        valIndex = {
        'Calinski-Harabasz':metrics.calinski_harabasz_score,
        'Silhouette':metrics.silhouette_score,
        'Davies-Bouldin':metrics.davies_bouldin_score
        }

        #message until this is fully implemented as I want. 
        messagebox.showinfo(title='INFORMATION', message='Working to create mono-clustering validation functionality')

        #log that user called MST
        logging.info(': User called Cluster validation function.')

        filename = filedialog.askopenfilename()
        try:
            data, col_groups = GB.readAndPreProcess(file=filename, transform = transform, scale =scale, func='CC')
        except BaseException:
            logging.error(': Unable to proceed, due to file error!')
            messagebox.showerror(title='Error',message='Unable to proceed, try again or return to homepage!')
            return
        
        #define the optimal cluster list for saving
        optClust = [None]*2
        optScores = np.ones((len(range(2,upperLimClust+1)),2))
        if linkage == 'ward':
            #calculate the distance matrix, which in this from should give me the sparse form. So, I should not need to use squareform
            dist = pdist(data, dissimilarity)

            # dist = squareform(dist)
            link_mat = ward(dist)

            #define the best score as 0 if Silhouette or Calinski-Harabasz
            if func == 'Davies-Bouldin':
                bestScore = 10**10
                for j in range(2, upperLimClust+1):
                    labels_ = fcluster(link_mat,j,criterion='maxclust')
                                
                    #update the best clustering solutions
                    optScores[j-2,1] = valIndex[func](data,labels_-1)
                    optScores[j-2,0] = j
                    if optScores[j-2,1] < bestScore:
                        optClust[0],optClust[1] = j, labels_-1
                        bestScore = optScores[j-2,1]
            else:
                bestScore = -1
                for j in range(2, upperLimClust+1):
                    labels_ = fcluster(link_mat,j,criterion='maxclust')
                                
                    #update the best clustering solutions
                    optScores[j-2,1] = valIndex[func](data,labels_-1)
                    optScores[j-2,0] = j
                    if optScores[j-2,1] > bestScore:
                        optClust[0],optClust[1] = j, labels_-1
                        bestScore = optScores[j-2,1]
        else:
            #calculate the distance matrix, which in this from should give me the sparse form. So, I should not need to use squareform
            dist = pdist(data, dissimilarity)
            # dist = distance[ensemble['Distance'][i]](data)

            for j in range(2, upperLimClust+1):
                #calculate the clustering solutions
                agglo = AC(n_clusters=j,linkage=linkage,
                        metric='precomputed').fit(dist)
                
                if func == 'Davies-Bouldin':
                    bestScore = 10**10
                    #update the best clustering solutions
                    optScores[j-2,1] = valIndex[func](data,labels_-1)
                    optScores[j-2,0] = j
                    if optScores[j-2,1] < bestScore:
                        optClust[0],optClust[1] = j, agglo.labels_
                        bestScore = optScores[j-2,1]
                else:
                    bestScore = -1
                    #update the best clustering solutions
                    optScores[j-2,1] = valIndex[func](data,labels_-1)
                    optScores[j-2,0] = j
                    if optScores[j-2,1] < bestScore:
                        optClust[0],optClust[1] = j, agglo.labels_
                        bestScore = optScores[j-2,1]
        plt.plot(optScores[:,0],optScores[:,1],'r-')
        plt.show()
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
            dataMST[i,2] = mstOut[int(dataMST[i,0]),int(dataMST[i,1])]

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

    def heatmapAnalysis(linkFunc,distMet,cmap,norm, colOrder = [], transform = 'None', scale ='None'):
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

        #read in the file
        messagebox.showinfo(message='Input file you would like to have heatmap of.')
        file = filedialog.askopenfilename()
        #send the data off to the readAndPreProcess function for analysis. 
        data, col_groups = GB.readAndPreProcess(file=file,transform=transform,scale=scale,func="CC")

        del(col_groups)
        #create messagebox explaining to users how they need to select clusters.
        print('Select clusters of interest, cluster and peak to pathway files will be automatically generated!')

        #Create the appropriate plt figure to allow for the comparison of linkage functions
        fig, axes = plt.subplots(1,1,figsize=(8,8))

        #find the linkages
        linkageOne = linkage(data,linkFunc,metric=distMet)

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
        linkageG = linkage(groupCluster,linkFunc,metric=distMet)
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

        #build
        sns.heatmap(dataFinal,fmt='0.2f', cmap=cmap, yticklabels=False)
        plt.show()

        return
        


    def selectClusters(link,dist,norm=0, colOrder=[], transform = 'None', scale = 'None',cmap = 'viridis'):
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
        if norm == 0:
            data, col_groups = GB.readAndPreProcess(file=file,transform=transform,scale=scale,func="CC")
        elif norm == 1:
            data, col_groups = GB.readAndPreProcess(file=file,transform=transform,scale=scale,func="CC",first =colOrder[0])

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

        if norm == 0:
            for i in range(data.shape[1]):
                #rearranging the data for heatmap
                for j in range(data.shape[0]):
                    #going through the metabolites
                    dataFinal[j,i] = data[metaboliteDendLeaves[j],groupDendLeaves[i]]

        else:
            colOrder = [int(i) for i in colOrder]
            colOrder = [i-1 for i in colOrder]
            for i in range(data.shape[1]):
                #rearranging the data for heatmap
                for j in range(data.shape[0]):
                    #going through the metabolites
                    dataFinal[j,i] = data[metaboliteDendLeaves[j],colOrder[i]]


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

    def enzymeLookUp(numSheets):
        '''
        '''

        #have the user select the file they would like to have read in.
        filename = filedialog.askopenfilename()
        
        #heatmapEnzyme Outputs
        outFile = 'HeatmapEnzyme.xlsx'
        writer = pd.ExcelWriter(outFile, engine='xlsxwriter')

        messagebox.showinfo(title="Starting Enzyme Look-Up", message="Please click ok to start function, this will take some time please be patient.")
        for i in range(numSheets):
            #read in each sheet
            dataCur = pd.read_excel(filename,sheet_name=i)
            dFDict = {}
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
            writer.close()
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


    def normalityCheck(transform='None',scale='None'):
        '''
        '''

        #have user input the wanted file
        file =filedialog.askopenfilename()

        #read in the metabolites file, and reshape for the analysis
        data, col_groups = GB.readAndPreProcess(file=file,transform=transform,scale=scale,func='CC')
        dataOrig,cols =  GB.readAndPreProcess(file=file,transform='None',scale='None',func='CC')
        data = np.reshape(data,(data.size,1))
        dataOrig = np.reshape(dataOrig,(dataOrig.size,1))

        fig, axes = plt.subplots(1, 2, figsize=(10,5))
        go = sns.kdeplot(data=dataOrig,ax=axes[0])
        axes[0].set_title("Raw Data")
        axes[0].set_xlabel('Intensities')
        g = sns.kdeplot(data=data,ax=axes[1])
        axes[1].set_title("Transform and/or Scaled")
        axes[1].set_xlabel('Intensities')
        plt.savefig("Normalized.png",dpi=600,transparent=True)
        plt.show()

    def progenesis():
        '''
        '''

        #have the user select the progenesis output spreadsheet. 
        messagebox.showinfo(title='File selection',message="Select progenesis output file, with two sheets.")
        fileName = filedialog.askopenfilename()
        d1 = pd.read_excel(fileName, sheet_name=0)
        d2 = pd.read_excel(fileName, sheet_name=1)
        index = list(d2["Identifications"].where(d2["Identifications"]!=0).dropna().index)
        d2 = d2.iloc[index]

        #columns 
        d2 = d2[["m/z","Neutral mass (Da)","Retention time (min)","Identifications","Accepted Compound ID","Accepted Description"]]

        #mz from datasheet two
        mz_d = d2[["m/z"]].dropna().to_numpy()

        # read in search values (SV) were looking for
        messagebox.showinfo(title='File selection',message="Select excel workbook with mz, and comparisons.")
        SV_file = filedialog.askopenfilename()
        SV_pd_ = pd.read_excel(SV_file)
        SV_pd = SV_pd_["m.z"].to_numpy()

        LL = mz_d * np.ones(len(SV_pd))

        mz_SV = SV_pd[:]
        QL = mz_SV * np.ones([len(mz_d),1])

        # get tolerance values for each of the SV's 
        tol = mz_SV * (20/1000000)
        A = LL - QL
        # finds where the two are very close
        temp = np.where(abs(A) <= tol)
        looking_at = mz_d[temp[:][0]]

        #checking our matches against the input data file... typically a vip matches
        compDict = {}
        for i in range(0,len(SV_pd)):
            # for the chosen searching value
            
            # look at where it has matches
            temp_indexes = np.where(temp[1]==i)[0]
            d2_indexes =[]
            for j in range(0,len(temp_indexes)):
                # for each of the matches
                # look at where in d2 it matches
                temp_d2_indexes = temp[0][temp_indexes[j]]
                d2_indexes.append(temp_d2_indexes)
            
            # add all d2 indexes for each
            compDict[i] = d2_indexes

        matches = 0
        for i in range(0,len(compDict)):
            # for each entry in the dictionary
            d1_match = []
            for j in range(0,len(compDict[i])):
                # for each entry in the dictionary's array find its mz value in d2
                mz_d2 = mz_d[compDict[i][j]]
                # then look for that in d1
                d1_match.extend(d1["m/z"].where(d1["m/z"].to_numpy() == mz_d2).dropna().index)

            # get all the data from those matches
            match_DF = d1.loc[d1_match]

            #check for the number of matches present... if there aren't any matches move on the next iteration
            if len(match_DF)==0:
                continue
            else: 
                matches +=1


            match_DF=match_DF.sort_values("Score",ascending=False)
            match_DF=match_DF.drop(columns=["Compound","Not fragmented","Not identified"])
            sheetName = "Query List " + str(mz_SV[i])
            Query = len(d1_match)*[None]
            comp = len(d1_match)*[None]

            #set the first value in the Query to the mz_SV
            Query[0] = mz_SV[i]
            comp[0]  = SV_pd_['comparison'][i]
            match_DF.insert(0, "Query", Query, True)
            match_DF.insert(match_DF.shape[1],"Comparisons",comp,True)

            #place the matches into a matrix of interest
            if matches == 1:
                final_df = match_DF
            else:
                final_df=pd.concat([final_df,match_DF])

        if matches == 0:
            messagebox.showinfo(message='No Matches found!')
        else:
            final_df.to_excel('VIPwComps.xlsx',index=False)
            messagebox.showinfo(message='Saved matching queries to VIPwComps.xlsx')

        return
    
    
    def mfgUtil(fileName,variant):
        '''
        The goal of this function is to streamline the generation of files for analysis in MetaboAnalyst. 

        Input: 

        fileName - (path or name of file in current directory) name of the file that you want to split into MetaboAnalyst Files
        variant - "All", "Multi", or "Uni". Allows the user to specify the comparisons they would like the function to create files for. Either all possible, only multivariate or only univariate. 

        Output:
        Comma separated values for the wanted comparisons. These files can be directly submitted to the metaboBot for analysis. 

        - Connor Boone, 2023
        '''
        logging.info(": called MetaboAnalyst File generation functions,")
        logging.info(": User selected {variant}")
        dir = os.getcwd()

        excelFile = fileName
        data = pd.read_excel(excelFile)
        df = pd.DataFrame(data)
        df = df.drop("rtmed",axis=1)
        dfOld = df

        numCols = sum(data["rtmed"].isna())

        # DO THIS FOR EACH ROW TO GET ALL UNIQUE IDENTIFIERS
        identifierDict = {}

        for i in range(0,numCols):
            uni=pd.unique(df.loc[i,:])
            title = uni[0]
            levels = uni[1:len(uni)]
            
            identifierDict[title] = list(levels)

        comparisons = list(identifierDict.values())
        logging.info("Comparing: " + str(comparisons))
        uni_F1 = []*len(comparisons)
        for k in range(0,len(comparisons)):
            temp = []
            # get all unique combinations
            for i in range(2,len(comparisons[k])+1):
                    temp.extend(combinations(comparisons[k],i))

            uni_F1.append(temp)

        # combine everything we need to compare
        comparitors = sum(comparisons,[])

        tempDict = {}
        myDict = {}
        lookupTable = {}
        iter = 0
        # for all comparitors
        for i in range(0,len(comparisons)):
            # for each level
            for j in range(0,len(comparisons[i])):                
                # for each level, look at the comparison list
                for m in range(0,len(comparisons)):
                    
                    # skip the comparison list that includes the selected level
                    if m==i:
                        pass
                    else:
                        tempDict = {}
                        currentIndex = 0

                        # run through each comparison and level 
                        for k in range(0,len(uni_F1[m])):
                            
                            #determine the types of files the user wanted out. 
                            if variant == "All":
                                tempDict[currentIndex] = uni_F1[m][k]
                                currentIndex += 1
                            elif variant == "Multi" and len(uni_F1[m][k]) > 2:                                
                                tempDict[currentIndex] = uni_F1[m][k]
                                currentIndex += 1
                            elif variant == "Uni" and len(uni_F1[m][k]) == 2:                                
                                tempDict[currentIndex] = uni_F1[m][k]
                                currentIndex += 1
                            pass
                # look where the level were lookin at shows up in the comparitors list
                tempIndex = comparitors.index(comparisons[i][j])
                # then save the temp dictionary
                myDict[tempIndex] = tempDict

                # save this index in the lookup table
                lookupTable[tempIndex] = comparisons[i][j]

        # for each comparitor
        for j in range(0,len(myDict)):
            # for each comparitor's variate sets
            if not len(myDict[j]) == 0:
                # if the current index is not empty 
                for k in range(0,len(myDict[j])):
                    os.chdir(dir)
            
                    data = pd.read_excel(excelFile)
                    df = pd.DataFrame(data)
                    df = df.drop("rtmed",axis=1)
                    # the variate combination we're looking at
                    variateCombo = myDict[j][k]

                    variateList = []
                    fileName = lookupTable[j] + "_"
                    # creates the csv file name using the variates in the file
                    for i in range(0,len(variateCombo)):
                        fileName += variateCombo[i]
                        if i+1 != len(variateCombo):
                            fileName += "_"
                    
                    variateList.append(lookupTable[j])
                    variateList.extend(list(variateCombo))
                    
                    tempArray = variateList
                    # Injury or Age
                    selected = tempArray[0]
                    
                    # get the number of rows with factors
                    numCols = sum(data["rtmed"].isna())

                    data_a = data
                    
                    myKeys = list(identifierDict.keys())
                    index = -1
                    
                    # finds the index where the selected variant is
                    for i in range(0,len(myKeys)):
                        for strings in identifierDict[myKeys[i]]:
                            if strings == selected:
                                index = i
                                break
                    
                    myIndex = myKeys.index(myKeys[index]);
                    #remove NA's
                    df = df[list(df.iloc[myIndex].dropna().index)] # removes any NA
                    rowNumber = df.loc[df['mz'] == myKeys[index]].index # gets the row number of the selected header                  
                    
                    notNeeded = []
                    for k in range(0,len(myKeys)):
                        rowVals = df.iloc[k] # gets the values in row k
                        
                        # checks to see if we care about this data based on the varaint selection above
                        temp = []
                        for i in range(0,len(rowVals)):
                            if rowVals[i] in tempArray:
                               
                                temp.extend([True])
                            else:
                                temp.extend([False])
                        # temp
                        temp = np.array(temp)
                        a = np.where(temp != True)[0]
                        
                        notNeeded.extend(list(a))

                    df = df.drop([index]) # drops the row of data we dont care about
                    data_a = df.drop(columns=data_a.columns[notNeeded]) # removes the columns with data we dont care about
                    
                    fileName = ""
                    # creates the csv file name using the variates in the file
                    for i in range(0,len(tempArray)):
                        fileName += tempArray[i]
                        if i+1 != len(tempArray):
                            fileName += "_"
            
                    # create a new directory if the name below isn't present. 
                    newpath = 'CSV_Files' # name of new folder
                    if not os.path.exists(newpath): # check to see if folder exists
                        os.makedirs(newpath) # if it does not, create a new folder
                    os.chdir(newpath) # move to the new folder
            
                    # determines whether the chosen comparison is uni or multi variate 
                    if len(tempArray) == 3:
                        newFolderName = "Uni"
                    else:
                        newFolderName = "Multi"
            
                    # creates the directory
                    if not os.path.exists(newFolderName):
                        os.makedirs(newFolderName)
                    os.chdir(newFolderName)
                    # saves the file to that directory
                    data_a = data_a.T
                    data_a.to_csv(fileName + ".csv",header=False)
                    
                    # goes back out to the "CSV_Files" directory
                    os.chdir("..")
        os.chdir(dir)
        logging.info(":Completed")
        return 

    def metaboBot(analysis='Uni',varianceFilter='SD',sampleNorm='None',trans='Log transformation',scale='Auto Scaling'):
        '''

        '''
        start =time.time()
        #set up the browser runner, update the directory, get all the csv files within the folder of interest
        directory = filedialog.askdirectory()

        curDir = os.getcwd()
        os.chdir(directory)
        files = glob.glob( '*.csv')
        driver = webdriver.Chrome()#options=chrome_options

        for i in range(len(files)):


            ################################################################################################################################
            ################################################################################################################################
            ################################################################################################################################
            ################################################################################################################################
            ################################################################################################################################
            ################################### Pre-processing before univariate or multi-variate analysis #################################
            #get file of interest and 
            file = files[i]
            # loc, name = os.path.split(file)
            name = file.strip('.csv')
            file = directory+'/'+file

            driver.get('https://www.metaboanalyst.ca/MetaboAnalyst/upload/StatUploadView.xhtml')


            #wait the homepage loads, if it doesn't load in 10 seconds quit, the expected file input should be a row (unpaired file) 
            try:
                #wait for the website to load and then select the Peak Intensities radio button. 
                element = WebDriverWait(driver,10).until(
                    EC.presence_of_element_located((By.XPATH,'//*[@id="j_idt12:j_idt20"]/div[3]/div'))
                )
                #select the radio button peak intensities
                driver.find_element(By.XPATH,'//*[@id="j_idt12:j_idt20"]/div[3]/div').click()
                #sent the selected file to the website
                driver.find_element(By.XPATH,'//*[@id="j_idt12:j_idt26_input"]').send_keys(file)
                #click the submit button. 
                driver.find_element(By.XPATH,'//*[@id="j_idt12:j_idt27"]').click()
            except:
                #quit the website, should the bot not work. 
                driver.quit()
                messagebox.showerror(message='I quit on homepage')
                return

            #wait for the data check to complete  
            try:
                #wait for element to show up
                element = WebDriverWait(driver,10).until(
                    EC.presence_of_element_located((By.XPATH,'//*[@id="form1:j_idt18"]'))
                )
                #click proceed, it is the onus of the user to ensure that the files are correct for metaboanalyst.
                driver.find_element(By.XPATH,'//*[@id="form1:j_idt18"]').click()
            except:
                #quit the website, should the bot not work.
                driver.quit()
                messagebox.showerror(message='I quit on data processing page')
                return

            #check for the correct button, and the appropriate filtering then move on. 
            try:
                #look for standard deviation filter... if it is present the rest will be as well. 
                element = WebDriverWait(driver,10).until(
                    EC.presence_of_element_located((By.XPATH,'//*[@id="j_idt14:j_idt24"]/div[2]/div/div/div[2]/span'))
                )
                if varianceFilter == 'SD':
                    #standard deviation filtering into pre-processing
                    driver.find_element(By.XPATH,'//*[@id="j_idt14:j_idt24"]/div[2]/div/div/div[2]/span').click()

                elif varianceFilter == 'MAD':
                    #select the median absolute deviation (MAD)
                    driver.find_element(By.XPATH,'//*[@id="j_idt14:j_idt24"]/div[3]/div/div/div[2]/span').click()

                elif varianceFilter == 'RSD':
                    #Relative standard deviation
                    driver.find_element(By.XPATH,'//*[@id="j_idt14:j_idt24"]/div[4]/div/div/div[2]/span').click()

                elif varianceFilter == 'MAD_m':
                    #Relative standard deviation
                    driver.find_element(By.XPATH,'//*[@id="j_idt14:j_idt24"]/div[5]/div/div/div[2]/span').click()
                
                #submit
                driver.find_element(By.XPATH,'//*[@id="j_idt14:j_idt41"]').click()
                #proceed
                driver.find_element(By.XPATH,'//*[@id="j_idt14:j_idt42"]').click()

            except:
                #quit the website, should the bot not work.
                driver.quit()
                #add message box here if you end up here.
                messagebox.showerror(message='I quit on the pre-processing page')

            #get the data normalized
            try:
                #check for the log-transformation button
                element = WebDriverWait(driver,10).until(
                    EC.presence_of_element_located((By.XPATH,'//*[@id="form1:j_idt74"]/div[2]/span'))
                )

                #select the correct button should the user select a sample normalization
                if sampleNorm == 'Sum':
                    #click on the normalization by sum 
                    driver.find_element(By.XPATH,'//*[@id="form1:j_idt41"]/div[2]/span').click()

                elif sampleNorm == 'Median':
                    #click on the normalization by median
                    driver.find_element(By.XPATH,'//*[@id="form1:j_idt44"]/div[2]/span').click()
                
                elif sampleNorm =='Quantile':
                    #click on the Quantile Normalization 
                    driver.find_element(By.XPATH,'//*[@id="form1:j_idt56"]/div[2]/span').click()


                #select the correct data transformation.
                if trans == 'Log transformation':
                    #click on log transformation
                    driver.find_element(By.XPATH,'//*[@id="form1:j_idt74"]/div[2]/span').click()

                elif trans == 'Square root transformation':
                    #click on square-root tranformation
                    driver.find_element(By.XPATH,'//*[@id="form1:j_idt78"]/div[2]/span').click()

                elif trans == 'Cube root transformation':
                    #click on the cube-root transformation
                    driver.find_element(By.XPATH,'//*[@id="form1:j_idt82"]/div[2]/span').click()

                #select the dataScale = ('None','Mean-Center','Auto-scale','Pareto-scale','Range-scale')correct data scaling
                if scale == 'Mean centering':
                    #click on mean centering
                    driver.find_element(By.XPATH,'//*[@id="form1:j_idt91"]/div[2]/span').click()

                elif scale == 'Auto Scaling':
                    #click on auto-scaling
                    driver.find_element(By.XPATH,'//*[@id="form1:j_idt94"]/div[2]/span').click()

                elif scale == 'Pareto Scaling':
                    #click on paret-scaling
                    driver.find_element(By.XPATH,'//*[@id="form1:j_idt97"]/div[2]/span').click()

                elif scale =='Range Scaling':
                    #click on range-scaling
                    driver.find_element(By.XPATH,'//*[@id="form1:j_idt100"]/div[2]').click()

                #make sure that the normalize button is present and selected each time. 
                element = WebDriverWait(driver,10).until(
                    EC.presence_of_element_located((By.XPATH,'//*[@id="form1:j_idt104"]/span'))
                )
                driver.find_element(By.XPATH,'//*[@id="form1:j_idt104"]/span').click()
            except:
                driver.quit()
                #add message here.
                messagebox.showerror(message='I quit on the data normalization page') 
                return

            #proceed to next step after pre-processing
            try:
                #the element is present but not clickable. 
                time.sleep(2)
                #wait for element to show up
                element = WebDriverWait(driver,15).until(
                    EC.presence_of_element_located((By.XPATH,'//*[@id="form1:nextBn"]/span'))
                )
                driver.find_element(By.XPATH,'//*[@id="form1:nextBn"]/span').click()  

            except:
                driver.quit()
                #add appropriate message
                messagebox.showerror(message='I quit on the preprocessing page')
                return

            ################################################################################################################################
            ################################################################################################################################
            ################################################################################################################################
            ################################################################################################################################
            ################################################################################################################################
            ###########################################Univariate or multi-variate analysis ################################################
            if analysis == 'Uni':
                #perform fold-change
                try:
                #   wait for element to show up
                    element = WebDriverWait(driver,10).until(
                        EC.presence_of_element_located((By.XPATH,'//*[@id="j_idt12"]/table/tbody/tr[2]/td/table/tbody/tr[2]/td/table/tbody/tr[1]/td/table/tbody/tr/td[1]/a'))
                    )
                    driver.find_element(By.XPATH,'//*[@id="j_idt12"]/table/tbody/tr[2]/td/table/tbody/tr[2]/td/table/tbody/tr[1]/td/table/tbody/tr/td[1]/a').click()
                except:
                    driver.quit()
                    #add appropriate message
                    messagebox.showerror(message='I quit on the fold-change selection page')
                    return

                #peform a t-test
                try:
                    element = WebDriverWait(driver,10).until(
                        EC.presence_of_element_located((By.XPATH,'//*[@id="treeForm:j_idt78:3_1"]/div'))
                    )
                    driver.find_element(By.XPATH,'//*[@id="treeForm:j_idt78:3_1"]/div').click()
                except:
                    driver.quit()
                    messagebox.showerror(message='I quit on the t-test page')
                    return


                #perform volcano plot analysis
                try:
                    element = WebDriverWait(driver,10).until(
                        EC.presence_of_element_located((By.XPATH,'//*[@id="treeForm:j_idt93:3_2"]/div'))
                    )
                    driver.find_element(By.XPATH,'//*[@id="treeForm:j_idt93:3_2"]/div').click()

                    element = WebDriverWait(driver,10).until(
                        EC.presence_of_element_located((By.XPATH,'//*[@id="form3:j_idt51"]'))
                    )
                    driver.find_element(By.XPATH,'//*[@id="form3:j_idt51"]').clear()
                    driver.find_element(By.XPATH,'//*[@id="form3:j_idt51"]').send_keys(0.05)
                    driver.find_element(By.XPATH,'//*[@id="form3:j_idt52"]/div[2]/div/div[2]/span').click()
                    driver.find_element(By.XPATH,'//*[@id="form3:j_idt59"]').click()

                except:
                    driver.quit()
                    messagebox.showerror(message='I quit on the volcano plot analyses page')
                    return

                #principal components analysis.
                try:
                    time.sleep(4)
                    #wait for element
                    element = WebDriverWait(driver,10).until(
                        EC.presence_of_element_located((By.XPATH,'//*[@id="treeForm:j_idt112:3_7"]/div'))
                    )
                    driver.find_element(By.XPATH,'//*[@id="treeForm:j_idt112:3_7"]/div').click()
                    #wait for element
                    element = WebDriverWait(driver,10).until(
                        EC.presence_of_element_located((By.XPATH,'//*[@id="ac"]/ul/li[3]/a'))
                    )
                    driver.find_element(By.XPATH,'//*[@id="ac"]/ul/li[3]/a').click()
                except:
                    driver.quit()
                    messagebox.showerror(message='I quit on the principal components analysis page')
                    return

                #PLS-DA
                try:
                    
                    #wait for element
                    element = WebDriverWait(driver,10).until(
                        EC.presence_of_element_located((By.XPATH,'//*[@id="treeForm:j_idt167:3_8"]/div'))
                    )
                    driver.find_element(By.XPATH,'//*[@id="treeForm:j_idt167:3_8"]/div').click()

                except:
                    driver.quit()
                    #add needed message.
                    messagebox.showerror(message='I quit on initial PLS-DA')
                    return

                #PLS-DA extra's
                try:
                    time.sleep(2)
                    #wait for element
                    element = WebDriverWait(driver,10).until(
                        EC.presence_of_element_located((By.XPATH,'//*[@id="ac"]/ul/li[2]/a'))
                    )
                    
                    driver.find_element(By.XPATH,'//*[@id="ac"]/ul/li[2]/a').click()
                    #wait for element
                    element = WebDriverWait(driver,10).until(
                        EC.presence_of_element_located((By.XPATH,'//*[@id="ac"]/ul/li[4]/a'))
                    )
                    driver.find_element(By.XPATH,'//*[@id="ac"]/ul/li[4]/a').click()

                except:
                    driver.quit()
                    #add appropriate message
                    messagebox.showerror(message='I quit on PLS-DA extras')
                    return

                #create a dendrogram 
                try:
                    element = WebDriverWait(driver,10).until(
                        EC.presence_of_element_located((By.XPATH,'//*[@id="treeForm:j_idt220:3_13"]/div'))
                    )
                    driver.find_element(By.XPATH,'//*[@id="treeForm:j_idt220:3_13"]/div').click()

                except:
                    driver.quit()
                    #add appropriate message
                    messagebox.showerror(message='I quit at the dendrogram')
                    return

                #create a heatmap
                try:
                    element = WebDriverWait(driver,10).until(
                        EC.presence_of_element_located((By.XPATH,'//*[@id="treeForm:j_idt35:3_14"]/div'))
                    )
                    driver.find_element(By.XPATH,'//*[@id="treeForm:j_idt35:3_14"]/div').click()

                except:
                    driver.quit()
                    #add appropriate message
                    messagebox.showerror('I quit at the heatmap')
                    return
                
                # go to the download page
                try:
                    element = WebDriverWait(driver,10).until(
                        EC.presence_of_element_located((By.XPATH,'//*[@id="treeForm:j_idt127:4"]/div'))
                    )
                    time.sleep(3)
                    driver.find_element(By.XPATH,'//*[@id="treeForm:j_idt127:4"]/div').click()

                except:
                    driver.quit()
                    #add appropriate message
                    messagebox.showerror(message='I didn''t make it to the download page' )
                    return


            elif analysis == 'Multi':
                #ANOVA
                try:
                    #the element is present but not clickable. 
                    #wait for element to show up
                    element = WebDriverWait(driver,15).until(
                        EC.presence_of_element_located((By.XPATH,'//*[@id="j_idt12"]/table/tbody/tr[2]/td/table/tbody/tr[2]/td/table/tbody/tr[2]/td/a'))
                    )
                    driver.find_element(By.XPATH,'//*[@id="j_idt12"]/table/tbody/tr[2]/td/table/tbody/tr[2]/td/table/tbody/tr[2]/td/a').click()  

                except:
                    driver.quit()
                    #add appropriate message
                    messagebox.showerror(message=' I quit on ANOVA page')
                    return


                #principal components analysis.
                try:
                    #wait for element
                    element = WebDriverWait(driver,10).until(
                        EC.presence_of_element_located((By.XPATH,'//*[@id="treeForm:j_idt92:3_7"]/div'))
                    )
                    driver.find_element(By.XPATH,'//*[@id="treeForm:j_idt92:3_7"]/div').click()
                    #wait for element
                    element = WebDriverWait(driver,10).until(
                        EC.presence_of_element_located((By.XPATH,'//*[@id="ac"]/ul/li[3]/a'))
                    )
                    driver.find_element(By.XPATH,'//*[@id="ac"]/ul/li[3]/a').click()
                except:
                    driver.quit()
                    messagebox.showerror(message='I quit on the principal components analysis page')
                    return


                #PLS-DA
                try:
                    
                    #wait for element
                    element = WebDriverWait(driver,10).until(
                        EC.presence_of_element_located((By.XPATH,'//*[@id="treeForm:j_idt167:3_8"]/div'))
                    )
                    driver.find_element(By.XPATH,'//*[@id="treeForm:j_idt167:3_8"]/div').click()

                except:
                    driver.quit()
                    #add needed message.
                    messagebox.showerror(message='I quit on initial PLS-DA')
                    return

                #PLS-DA extra's
                try:
                    time.sleep(2)
                    #wait for element
                    element = WebDriverWait(driver,10).until(
                        EC.presence_of_element_located((By.XPATH,'//*[@id="ac"]/ul/li[2]/a'))
                    )
                    
                    driver.find_element(By.XPATH,'//*[@id="ac"]/ul/li[2]/a').click()
                    #wait for element
                    element = WebDriverWait(driver,10).until(
                        EC.presence_of_element_located((By.XPATH,'//*[@id="ac"]/ul/li[4]/a'))
                    )
                    driver.find_element(By.XPATH,'//*[@id="ac"]/ul/li[4]/a').click()

                except:
                    driver.quit()
                    #add appropriate message
                    messagebox.showerror(message='I quit on PLS-DA extras')
                    return

                #create a dendrogram 
                try:
                    element = WebDriverWait(driver,10).until(
                        EC.presence_of_element_located((By.XPATH,'//*[@id="treeForm:j_idt220:3_13"]/div'))
                    )
                    driver.find_element(By.XPATH,'//*[@id="treeForm:j_idt220:3_13"]/div').click()

                except:
                    driver.quit()
                    #add appropriate message
                    messagebox.showerror(message='I quit at the dendrogram')
                    return

                #create a heatmap
                try:
                    element = WebDriverWait(driver,10).until(
                        EC.presence_of_element_located((By.XPATH,'//*[@id="treeForm:j_idt35:3_14"]/div'))
                    )
                    driver.find_element(By.XPATH,'//*[@id="treeForm:j_idt35:3_14"]/div').click()

                except:
                    driver.quit()
                    #add appropriate message
                    messagebox.showerror(message='I quit at the heatmap')
                    return

 
                # go to the download page
                try:
                    element = WebDriverWait(driver,10).until(
                        EC.presence_of_element_located((By.XPATH,'//*[@id="treeForm:j_idt127:4"]/div'))
                    )
                    driver.find_element(By.XPATH,'//*[@id="treeForm:j_idt127:4"]/div').click()

                except:
                    driver.quit()
                    #add appropriate message
                    messagebox.showerror(message='I didn''t make it to the download page' )
                    return
                
            ################################################################################################################################
            ################################################################################################################################
            ################################################################################################################################
            ################################################################################################################################
            ################################################################################################################################
            ########################################### Downloading and saving results. ####################################################

            ##### After analysis trying to download the data.  
            #download results
            try:
                element = WebDriverWait(driver,10).until(
                    EC.presence_of_element_located((By.XPATH,'//*[@id="ac:form1:j_idt20_data"]/tr[1]/td[1]/a'))
                )
                driver.find_element(By.XPATH,'//*[@id="ac:form1:j_idt20_data"]/tr[1]/td[1]/a').click()

            except:
                driver.quit()
                #Add appropriate message
                messagebox.showerror(message='I didn''t download the data. ')
                return


            #rename zip file, but first check that it has downloaded.
            rename = name+'.zip'
            basepath = os.path.expanduser('~')
            basepath +='/Downloads/Download.zip'
            #get the start time to give plenty of time for the download to occur. 
            downloadTime = 0
            downloaded = False
            start = time.time()
            while downloadTime < 60:
                #checking for download.
                downloadTime = time.time() - start
                if os.path.exists(basepath):
                    #file has downloaded
                    downloaded = True
                    break

            if downloaded:
                time.sleep(2)
                #rename the file and unzip into the appropriate folder
                os.rename(basepath,rename)
                directory = os.getcwd()
                zip_dir = directory + '/'+ name
                with zipfile.ZipFile(rename,'r') as zip_ref:
                    zip_ref.extractall(zip_dir)
                #remove the downloaded and renamed zip file. 
                os.remove(rename)

        logging.info(' Completed the analysis for the selected inputs. ')
        messagebox.showinfo(message='Successfully completed runs for all files in the files in the selected directory.')

    def mummiBot(pval = 0.05,lcMode='Positive',db = 'Mouse (KEGG)'):
        '''
        
        '''
        #set up the browser runner, update the directory, get all the csv files within the folder of interest
        directory = filedialog.askdirectory()

        curDir = os.getcwd()

        os.chdir(directory)
        files = glob.glob( '*.csv')

        driver = webdriver.Chrome()

        for i in range(len(files)):
            #get file of interest and 
            file = files[i]
            name = file.strip('.csv')
            file = directory+'/'+file

            ## Navigate to the MetaboAnalyst Website for the one-factor stats analysis
            driver.get('https://www.metaboanalyst.ca/MetaboAnalyst/upload/PeakUploadView.xhtml')
            
            
            #determine the machine mode used for analysis. 
            if lcMode == 'Positive':
                try:
                    #wait for the website to load and then select the Peak Intensities radio button. 
                    element = WebDriverWait(driver,10).until(
                        EC.presence_of_element_located((By.XPATH,'//*[@id="ac:form2:j_idt19"]/div[2]'))
                    )
                    #click on the drop down menu
                    driver.find_element(By.XPATH,'//*[@id="ac:form2:j_idt19"]/div[2]').click()
                    driver.find_element(By.XPATH,'//*[@id="ac:form2:j_idt19_0"]').click() #click on the wanted mode, currently this is positive.

                except:
                    messagebox.showinfo(message='Failed to select the correct mode')
                    return


            #upload the files you would like to have analyzed
            driver.find_element(By.XPATH,'//*[@id="ac:form2:j_idt46_input"]').send_keys(file)
            #submit the file for analysis
            driver.find_element(By.XPATH,'//*[@id="ac:form2:j_idt48"]').click()


            #proceed after data sanity check
            try:
                element = WebDriverWait(driver,10).until(
                    EC.presence_of_all_elements_located((By.XPATH,'//*[@id="form1:j_idt18"]'))
                )
                driver.find_element(By.XPATH,'//*[@id="form1:j_idt18"]').click()

            except:
                messagebox.showerror(message='Failed to get past the data sanity check, check input sheet.')
                return

            
            
            try:
                #verify that the p-value cutoff is present and set it to 0.05
                element = WebDriverWait(driver,10).until(
                    EC.presence_of_element_located((By.XPATH,'//*[@id="j_idt13:j_idt25"]/tbody/tr[1]/td[2]/table/tbody/tr[1]/td[3]/table/tbody/tr/td/span/input'))
                )
                #clear out the p-value and verify that it is 0.05
                driver.find_element(By.XPATH,'//*[@id="j_idt13:j_idt25"]/tbody/tr[1]/td[2]/table/tbody/tr[1]/td[3]/table/tbody/tr/td/span/input').clear()
                driver.find_element(By.XPATH,'//*[@id="j_idt13:j_idt25"]/tbody/tr[1]/td[2]/table/tbody/tr[1]/td[3]/table/tbody/tr/td/span/input').send_keys(pval)
                
                #determine the database and click on the correct one. 
                if db == 'Mouse (KEGG)':
                    driver.find_element(By.XPATH,'//*[@id="j_idt13:j_idt110"]/div[2]/span').click()
                elif db == 'Human (BioCyc)':
                    driver.find_element(By.XPATH,'//*[@id="j_idt13:j_idt98"]/div[2]/span').click()
                elif db == 'Human (KEGG)':
                    driver.find_element(By.XPATH,'//*[@id="j_idt13:j_idt100"]/div[2]/span').click()
                elif db == 'Mouse (BioCyc)':
                    driver.find_element(By.XPATH,'//*[@id="j_idt13:j_idt108"]/div[2]/span').click()
                elif db == 'Rat (KEGG)':
                    driver.find_element(By.XPATH,'//*[@id="j_idt13:j_idt112"]/div[2]/span').click()
                elif db == 'Cow (KEGG)':
                    driver.find_element(By.XPATH,'//*[@id="j_idt13:j_idt120"]/div[2]/span').click()

                driver.find_element(By.XPATH,'//*[@id="j_idt13:j_idt439"]').click()
            except:
                messagebox.showerror(message='Failed to input the wanted database, and p-value.')
                return


            try:
                #making sure the download tab is available
                element = WebDriverWait(driver,10).until(
                    EC.presence_of_element_located((By.XPATH,'//*[@id="treeForm:j_idt77:4"]/div'))
                )
                driver.find_element(By.XPATH,'//*[@id="treeForm:j_idt77:4"]/div').click()

            except:
                messagebox.showerror(message='No download tab')
                return

            try:
                time.sleep(2)
                element = WebDriverWait(driver,10).until(
                    EC.presence_of_element_located((By.XPATH,'//*[@id="ac:form1:j_idt20_data"]/tr[1]/td[1]/a'))
                )
                driver.find_element(By.XPATH,'//*[@id="ac:form1:j_idt20_data"]/tr[1]/td[1]/a').click()

            except:
                messagebox.showerror(message='File did not start to download')
                return

            #rename zip file, but first check that it has downloaded.
            rename = name+'.zip'
            basepath = os.path.expanduser('~')
            basepath +='/Downloads/Download.zip'
            #get the start time to give plenty of time for the download to occur. 
            downloadTime = 0
            downloaded = False
            start = time.time()
            while downloadTime < 60:
                #checking for download.
                downloadTime = time.time() - start
                if os.path.exists(basepath):
                    #file has downloaded
                    downloaded = True
                    break

            if downloaded:
                time.sleep(2)
                #rename the file and unzip into the appropriate folder
                os.rename(basepath,rename)
                directory = os.getcwd()
                zip_dir = directory + '/'+ name
                with zipfile.ZipFile(rename,'r') as zip_ref:
                    zip_ref.extractall(zip_dir)
                #remove the downloaded and renamed zip file. 
                os.remove(rename)

        #close out the current browser window
        driver.close()


    def externalCriteria(comp='rand',trans ='None',scale='None'):
        '''

        '''

        #set up the dictionary for the distance measure.
        distance = {
        'CNS': GB.correlationNosqrt,
        'CS': GB.correlationSqrt,
        'PW': GB.pairWise
        }

        #setting the functions into a validation, and distance metrics. 
        valIndex = {
            'CH':GB.calinskiHarabasz_correlation,
            'SIL':metrics.silhouette_score,
            'DBI':GB.daviesBouldinScore_correlation
        }

        #In this functionality I need to read in the data, and ask for a .csv of the clustering parameters of interest similar to ensemble. 

        #read in the data file of interest.
        messagebox.showinfo(message='Select the input file of interest.')
        dataFile = filedialog.askopenfilename()
        data, col_groups = GB.readAndPreProcess(dataFile,transform=trans,scale=scale,func="CC")

        #ask the user for the clustering parameters they would like to use.
        #read in the csv, for ensemble creation
        messagebox.showinfo(message='Select clustering parameters file.')
        clustsFile = filedialog.askopenfilename()
        clusts = pd.read_csv(clustsFile)
        clusts_s = clusts[['Distance','Linkage']]
        distLink = tuple(clusts_s.itertuples(index=False,name=None))

        #a list for saving the best labels from each clustering
        best_labels = [None]*clusts.shape[0]
        
        #go through the clustering parameters of interest and generate the 
        for i in range(clusts.shape[0]):
            bestScore = 0
            optClust = [None]*2
            #calculate the distance matrix
            if clusts['Linkage'][i] == 'ward':
                #calculate the distance metric for the clustering of interest. 
                dist = distance[clusts['Distance'][i]](data,metric=clusts['Correlation'][i])
                dist = squareform(dist)
                link_mat = ward(dist)
                for j in range(10):
                    labels_ = fcluster(link_mat,j+2,criterion='maxclust')
                                
                    #update the best clustering solutions
                    score = valIndex[clusts['Optimizer'][i]](data,labels_,clusts['Distance'][i])
                    if score > bestScore:
                        optClust[0],optClust[1] = j+2, labels_-1
                        bestScore = score
                #save the labels that were best.
                best_labels[i]= optClust[1]   
                
            elif clusts['Optimizer'][i] =="DBI":
                #update such that best score is lowest
                bestScore =10**10

                #get the distance metric out
                dist = distance[clusts['Distance'][i]](data,metric=clusts['Correlation'][i])
                for j in range(10):
                    #calculate the clustering solutions
                    agglo = AC(n_clusters=j+2,linkage=clusts['Linkage'][i],
                            metric='precomputed').fit(dist)
                    
                    #update the best clustering solutions
                    score = valIndex[clusts['Optimizer'][i]](data,agglo.labels_+1,clusts['Distance'][i])
                    if score < bestScore:
                        optClust[0],optClust[1] = j+2, agglo.labels_
                        bestScore = score
                #save the labels that were best
                best_labels[i]= optClust[1]    

            else:
                #get the distance metric for clustering
                dist = distance[clusts['Distance'][i]](data,metric=clusts['Correlation'][i])      
                for j in range(10):
                    #calculate the clustering solutions
                    agglo = AC(n_clusters=j+2,linkage=clusts['Linkage'][i],
                            metric='precomputed').fit(dist)
                    
                    #update the best clustering solutions
                    score = valIndex[clusts['Optimizer'][i]](dist,agglo.labels_,metric='precomputed')
                    if score > bestScore:
                        optClust[0],optClust[1] = j+2, agglo.labels_
                        bestScore = score
                #save the labels that were best
                best_labels[i]= optClust[1]    

        #check for what is what and send labels. 
        if comp =='Rand-index':
            #send to the rand index function
            GB.randComp(best_labels,distLink)
            messagebox.showinfo(title='Completed',message='Successfully saved the output of the Rand-Index comparison.')

        elif comp =='Adjusted Rand-index':
            #send to the adjusted rand function
            GB.adjRandComp(best_labels,distLink)
            messagebox.showinfo(title='Completed',message='Successfully saved the output of the Adjusted Rand-Index comparison.')

        elif comp =='Normalized Mutual Info.':
            #send to the mutual info function
            GB.mutualInfo(best_labels,distLink,'norm')
        
        elif comp =='Adjusted Mutual Info.':
            #send to the mutual info function
            GB.mutualInfo(best_labels,distLink,'adj')

        else:
            #eventually add here so that I can support this functionality. 
            messagebox.showinfo(message='Currently not supporting the comparison between mono-clustering and ensemble clustering solutions.')


    def geneToPathway():
        '''
        '''


        #give the user the input that tells them what to submit
        messagebox.showinfo(message="Select a csv file of genes, Fold-change, and p-values")

        file = filedialog.askopenfilename()

        #read in the input parameters
        genes = pd.read_csv(file)

        #add 10 columns to genes data frame
        pathIndex = [f"P{i}" for i in range(1, 11)]
        newCols = pd.DataFrame(np.zeros((genes.shape[0],10)),columns=pathIndex)

        #add columns for 
        genes = pd.concat([genes, newCols],axis=1)
        genes = genes.replace(0, np.nan)

        # Specify the file path and name
        file_path = 'Results.xlsx'

        # Use ExcelWriter with the openpyxl engine
        with pd.ExcelWriter(file_path, engine='openpyxl') as writer:
            genes.to_excel(writer, index=False, sheet_name='Sheet1')
            
            # Load the openpyxl workbook object
            workbook = writer.book
            sheet = writer.sheets['Sheet1']
            columnN = 3
            count = 0
            for j in range(genes.shape[0]):
                request = REST.kegg_find("mmu",genes['Gene'][j])

                #write the file to the wanted location
                txtFCur = 'Testing'+ '.txt'
                open(txtFCur,'w').write(request.read())

                #open the file and read, parse out to search names
                geneSearch = open(txtFCur,"r")
                geneMatches = geneSearch.read()

                #start parsing
                geneList = geneMatches.split("\n")
                for i in range(len(geneList)):
                    #check the current gene_i for the name, removing tab
                    curCheck = geneList[i].split("\t")
                    if len(curCheck) == 1:
                        continue
                    #hello
                    curCheck1 = curCheck[1].split(",")

                    if genes['Gene'][j] in curCheck1:
                        count +=1  
                        requestGene = REST.kegg_get(curCheck[0])
                        #write the file to the wanted location
                        txtFCur = 'TestingGet'+ '.txt'
                        open(txtFCur,'w').write(requestGene.read())


                        #open the txtFCur,read split then look for start and end
                        check = open(txtFCur,"r")
                        check1 = check.read()
                        pathwayStart = check1.split("\n")
                        #get the starting and ending positions for pathway then get pathways. 
                        pathway = False
                        pathwayStartI = None
                        pathwayEndI = None
                        #search for start and end positions of pathway
                        for k in range(len(pathwayStart)):
                            if pathway == False:
                                #figure out where the start of the pathway is 
                                if 'PATHWAY' in pathwayStart[k]:
                                    pathway = True
                                    pathwayStartI = k 
                                else:
                                    continue

                            else:
                                if pathwayStart[k][0].isalpha():
                                    #keep as the wrong index for easy input to the range function
                                    pathwayEndI = k
                                    break

                        if pathwayEndI == None:
                            continue


                        #get out the pathways of interest
                        pathList = 10*[None]
                        count2 = -1

                        for m in range(pathwayStartI,pathwayEndI):
                            count2 +=1
                            if count2 == 10:
                                break
                            
                            if m == pathwayStartI:
                                #removing pathway and getting the list
                                pathList[count2] = pathwayStart[m].strip('PATHWAY').strip()
                            
                            else:
                                pathList[count2] = pathwayStart[m].strip()
            
                #save the hyperlinks to the appropriate location.
                columnN = 3
                for l in range(len(pathList)):
                    if pathList[l] == None:
                        continue
                    columnN += 1
                #go through each pathway and add to the spreadsheet.
                    # Use openpyxl to add a hyperlink
                    # Note: Excel uses 1-based indexing, adjust cell accordingly
                    sheet.cell(row=j+2, column=columnN).hyperlink = "https://www.genome.jp/pathway/" + pathList[l][:8]
                    sheet.cell(row=j+2, column=columnN).value = pathList[l]
                    sheet.cell(row=j+2, column=columnN).style = "Hyperlink"
