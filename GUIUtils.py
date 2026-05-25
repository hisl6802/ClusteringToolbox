#Creating a class containing functions that will be used in GUI
import time, glob, sys, logging, warnings, config, mplcursors
import pandas as pd
import numpy as np
import statistics as stat
from multiprocessing import Pool
import seaborn as sns
import GuiBackground as GB
import seaborn as sns
import matplotlib.pyplot as plt
import networkx as nx

from numpy import isin
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage, ward, fcluster
from scipy.spatial.distance import pdist, squareform
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import minimum_spanning_tree
from scipy.stats import t
from sklearn.cluster import AgglomerativeClustering as AC, KMeans
from sklearn import metrics
from itertools import combinations
from tkinter import filedialog, messagebox, simpledialog
from Bio.KEGG import REST, Compound
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

try:
    from pydeseq2.default_inference import DefaultInference
except ImportError:
    DefaultInference = None

# Suppress a noisy OpenMP runtime warning emitted via threadpoolctl during
# some sklearn calls on certain mixed runtime environments.
warnings.filterwarnings(
    "ignore",
    message=r"(?s).*Found Intel OpenMP \('libiomp'\).*LLVM OpenMP \('libomp'\).*",
    category=RuntimeWarning,
    module=r"threadpoolctl",
)

class GUIUtils:
    def dataIntegrity(file):
        '''
        The data integrity allows users to determine whether the first column 
        of data in a dataframe contains double decimals. This error mostly 
        arises when working with MetaboAnalyst results. 

        Input:

        file - input full path to the file, use tkinter's filedialog for 
        ease of getting file path. 

        Output: 
        
        This function outputs an excel file with corrected values, the original
        file name has  _corrected appended to the end. 

        '''
        #log that the user called the data integrity function
        logging.info(': User called the Data Integrity function.')

        try:
            #Read in Volcano Plot data
            if file[len(file)-1] == 'x':
                volcano = pd.read_excel(file)
                
            elif file[len(file)-1] == 'v':
                volcano = pd.read_csv(file)
        except:
            logging.error(': Failed to read in the excel file. Please put error in the Github issues tab.')
            messagebox.showerror(title='Error',message='Failed to read in the excel file. Please let Brady know!!')
            return
            
        # grab the first row of the volcano data (treat as a plain sequence,
        # ignoring the index labels to avoid KeyError on non-0-based indices)
        check = volcano.iloc[0]

        #create array that can save the fixed data and the data that did not need to be fixed
        correctedArray = np.zeros(check.shape[0])

        #search each of the volcano data rows to determine if they have double decimals.
        for i in range(check.shape[0]):
            # grab the value by position, not by label
            curVal = check.iloc[i]

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
        volcano.iloc[0] = correctedArray
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

    def createClustergram(norm, linkFunc, distMet, cmap, colOrder=[], transform='None', scale='None', file=None):
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

        fileI = file or filedialog.askopenfilename()
        
        try:
            if norm == 0 or norm == 2:
                data, col_groups = GB.readAndPreProcess(file=fileI,transform=transform,scale=scale,func="CC")

            elif norm == 1:
                data, col_groups = GB.readAndPreProcess(file=fileI,transform=transform,scale=scale,func="CC",first=colOrder[0]) 
        except TypeError:
            logging.error(': No file selected!')
            messagebox.showerror(title='Error',message='Error loading in the data, normilizing it or something else. Contact Brady if cannot figure out.')
            return

        #create dendrogram and plot data, remove the unneed data, log what the user has called and completed and return to the UI.        
        GB.create_dendrogram(data,col_groups, norm, link=linkFunc, dist=distMet,color = cmap,colOrder=colOrder)
        del(data,norm,linkFunc,distMet)
        logging.info(': Sucessfully created the wanted Clustergram')
        return

    def groupMedians(rmZeros=0, file=None):
        '''
        Determine the number of groups and then create a list or array of the appropriate
        beginning and ending of each group. This assumes that the groups are all of equal size which should be
        the goal for any and all analysis. Groups with out the same sizes should be considered
        inappropriate for analysis in this context, additionally it should be noted that statistics
        with out the same groups sizes can lead to incorrect analysis.
        
        Input:
        groupMedians will prompt you to select a file which you would like to have medians determined for. 

        Output:
        This function outputs a .csv file with _Medians appended to the end of the original file name. *** This will soon be updated to
        an excel file output for easy calling and input to the createClustergram, MST or Ensemble clustering functions. 
        '''

        #log that the user called the group medians function
        logging.info(': User called the Group Medians function.')
        file = file or filedialog.askopenfilename()

        #getting the needed data. 
        df = pd.read_excel(file)
        rt = df[list(df.columns)[-1]]
        mz = df[list(df.columns)[0]]

        #transpose the matrix to get what I want.
        df = df.T
        out = df.groupby([0]).median()
        out = out.T
        #create a mz, and rtmed column
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
        #create the file prefixes and suffixes for saving the resulting comparisons.
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
            
    def compoundMatchUp(typeFile='all', file=None):
        '''
        The compoundMatchUp function is responsible for matching the output compounds from mummichog to compounds from the KEGG data base spreadsheet. 

        Input:
        compoundMatchUp expects the user to select the directory or the file they would like to run a compound match-up on. 
        The file "mummichog_matched_compound_all.csv" needs to be updated prior to running this function. 
        '''
        
        logging.info(': Compound Match-Up function called!')

        # Pull in matched-compound data from either:
        #  - a single CSV file, or
        #  - a directory containing many matched-compound CSV files (from MummiBot)
        if file is None:
            if typeFile == 'all':
                # Ask for folder containing Mummichog output (from MummiBot), not P2P files
                file = filedialog.askdirectory(
                    title="Select folder containing Mummichog output (mummichog_matched_compound_*.csv)\nRun MummiBot on your P2P files first if you don't have these."
                )
            else:
                file = filedialog.askopenfilename(
                    title="Select Mummichog matched-compound CSV",
                    filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
                )
            if not file:
                return

        if typeFile == 'all':

            # ------------------------------------------------------------------
            # CASE 1: user supplied a DIRECTORY -> batch mode
            # ------------------------------------------------------------------
            if os.path.isdir(file):
                paths = [
                    os.path.join(file, f)
                    for f in os.listdir(file)
                    if f.lower().endswith(".csv")
                ]
                if not paths:
                    messagebox.showerror(
                        title="Error",
                        message="No CSV files were found in the selected folder."
                    )
                    return

                # Only process files that have Matched.Compound (Mummichog output)
                processed = 0
                # loop over each matched-compound CSV and write a corresponding
                # CompoundMatchUps_<basename>.csv next to it
                for path in paths:
                    try:
                        my_data = pd.read_csv(path)
                    except Exception:
                        logging.error(f": Failed to read Mummichog CSV: {path}")
                        continue

                    if "Matched.Compound" not in my_data.columns:
                        logging.warning(f": Skipping file without 'Matched.Compound' column: {path}")
                        continue

                    processed += 1
                    my_final_data = np.zeros((len(my_data["Matched.Compound"]), 2))
                    my_final_data = pd.DataFrame(my_final_data, columns=['ID', 'Compound Name'])
                    #grab the compound ID of interest
                    lenCompounds = len(my_data['Matched.Compound'])
                    for i in range(len(my_data["Matched.Compound"])):
                        if (i+1)%100 ==0:
                            x = ((i+1)/lenCompounds)*100
                            x = float("{0:.2f}".format(x))
                            logging.info(': ' + str(x)+'%' + ' completed!')

                        cur_id = str(my_data['Matched.Compound'][i])

                        #input the values into a request from KEGG API
                        if cur_id and cur_id[0] == 'G':
                            my_final_data.at[i, 'ID'] = cur_id
                            my_final_data.at[i, 'Compound Name'] = 'GAG subunits'

                        elif cur_id and cur_id[0] == 'C':
                            if len(cur_id) > 1 and cur_id[1] == 'E':
                                my_final_data.at[i, 'ID'] = cur_id
                                my_final_data.at[i, 'Compound Name'] = 'Not in KEGG, will update soon!'

                            else:
                                try:
                                    request = REST.kegg_get(cur_id)
                                    txtFCur = cur_id + '.txt'
                                    open(txtFCur,'w').write(request.read())
                                    records = Compound.parse(open(txtFCur))
                                    record = list(records)[0]
                                    os.remove(txtFCur)
                                    my_final_data.at[i, 'ID'] = cur_id
                                    my_final_data.at[i, 'Compound Name'] = record.name

                                except:
                                    logging.error(": No KEGG match found! Check KEGG Website!")
                                    logString = cur_id
                                    logString = ': Failed to find-' + logString
                                    logging.info(logString)
                                    my_final_data.at[i, 'ID'] = cur_id
                                    my_final_data.at[i, 'Compound Name'] = 'No match in KEGG, investigate this compound further.'
                                    continue

                        else:
                            my_final_data.at[i, 'ID'] = cur_id
                            my_final_data.at[i, 'Compound Name'] = 'Unknown'

                    base, fname = os.path.split(path)
                    root, _ = os.path.splitext(fname)
                    out_name = f"CompoundMatchUps_{root}.csv"
                    out_path = os.path.join(base, out_name)
                    my_final_data.to_csv(path_or_buf=out_path, index=False)

                if processed == 0:
                    messagebox.showerror(
                        title="No Mummichog output found",
                        message="No CSV files in this folder contain a 'Matched.Compound' column.\n\n"
                                "Compound ID requires Mummichog output (e.g. mummichog_matched_compound_all.csv), "
                                "not Peaks-to-Pathways (P2P) files.\n\n"
                                "Run MummiBot on your P2P files first, then select the folder where MummiBot saved the downloaded CSVs."
                    )
                    return
                messagebox.showinfo(
                    title="Success",
                    message=f"Successfully generated CompoundMatchUps for {processed} file(s) in the folder."
                )
                return

            # ------------------------------------------------------------------
            # CASE 2: user supplied a SINGLE CSV -> keep legacy filename
            # ------------------------------------------------------------------
            else:
                path = file
                try:
                    my_data = pd.read_csv(path)
                except Exception:
                    logging.error(f": Failed to read Mummichog CSV: {path}")
                    messagebox.showerror(
                        title="Error",
                        message=f"Failed to read CSV file:\n{path}"
                    )
                    return

                if "Matched.Compound" not in my_data.columns:
                    messagebox.showerror(
                        title="Wrong file type",
                        message="This file is not a Mummichog matched-compound output.\n\n"
                                "Compound ID needs the output from MummiBot (e.g. mummichog_matched_compound_all.csv), "
                                "which has a 'Matched.Compound' column.\n\n"
                                "Run MummiBot on your P2P (Peaks-to-Pathways) files first, then select the Mummichog CSV or its folder."
                    )
                    return

                my_final_data = np.zeros((len(my_data["Matched.Compound"]), 2))
                my_final_data = pd.DataFrame(my_final_data, columns=['ID', 'Compound Name'])
                #grab the compound ID of interest
                lenCompounds = len(my_data['Matched.Compound'])
                for i in range(len(my_data["Matched.Compound"])):
                    if (i+1)%100 ==0:
                        x = ((i+1)/lenCompounds)*100
                        x = float("{0:.2f}".format(x))
                        logging.info(': ' + str(x)+'%' + ' completed!')

                    cur_id = str(my_data['Matched.Compound'][i])

                    #input the values into a request from KEGG API
                    if cur_id and cur_id[0] == 'G':
                        my_final_data.at[i, 'ID'] = cur_id
                        my_final_data.at[i, 'Compound Name'] = 'GAG subunits'

                    elif cur_id and cur_id[0] == 'C':
                        if len(cur_id) > 1 and cur_id[1] == 'E':
                            my_final_data.at[i, 'ID'] = cur_id
                            my_final_data.at[i, 'Compound Name'] = 'Not in KEGG, will update soon!'

                        else:
                            try:
                                request = REST.kegg_get(cur_id)
                                txtFCur = cur_id + '.txt'
                                open(txtFCur,'w').write(request.read())
                                records = Compound.parse(open(txtFCur))
                                record = list(records)[0]
                                os.remove(txtFCur)
                                my_final_data.at[i, 'ID'] = cur_id
                                my_final_data.at[i, 'Compound Name'] = record.name

                            except:
                                logging.error(": No KEGG match found! Check KEGG Website!")
                                logString = cur_id
                                logString = ': Failed to find-' + logString
                                logging.info(logString)
                                my_final_data.at[i, 'ID'] = cur_id
                                my_final_data.at[i, 'Compound Name'] = 'No match in KEGG, investigate this compound further.'
                                continue

                    else:
                        my_final_data.at[i, 'ID'] = cur_id
                        my_final_data.at[i, 'Compound Name'] = 'Unknown'

                # For single-file usage, match the original behavior and
                # write CompoundMatchUps.csv to the current working directory.
                out_path = os.path.join(os.getcwd(), "CompoundMatchUps.csv")
                my_final_data.to_csv(path_or_buf=out_path, index=False)

                messagebox.showinfo(
                    title="Success",
                    message=f"Successfully generated {os.path.basename(out_path)}."
                )
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
        Excel sheet with exact monoisotopic masses. 

        Output:
        Updated excel sheet with KEGG Compound matches

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
    
    def ensembleClusteringFullOpt(numClusts=10, transform='None', scale='None', file=None,
                                  colMap='viridis'):
        '''
        colMap - matplotlib colormap name for ensemble clustermaps.
        '''

        #log that the user called ensemble clustering function
        logging.info(': User called Ensemble Clustering function.')

        ENSEMBLE_PARAMS = [
            {'Distance': 'CNS', 'Linkage': 'ward', 'Correlation': 'spearman', 'Optimizer': 'CH'},
            {'Distance': 'CNS', 'Linkage': 'average', 'Correlation': 'spearman', 'Optimizer': 'SIL'},
            {'Distance': 'CS', 'Linkage': 'ward', 'Correlation': 'spearman', 'Optimizer': 'CH'},
            {'Distance': 'CS', 'Linkage': 'average', 'Correlation': 'spearman', 'Optimizer': 'SIL'},
            {'Distance': 'CS', 'Linkage': 'complete', 'Correlation': 'spearman', 'Optimizer': 'DBI'},
            {'Distance': 'CNS', 'Linkage': 'complete', 'Correlation': 'spearman', 'Optimizer': 'DBI'}
                                ]
            
        parameters = pd.DataFrame(ENSEMBLE_PARAMS)
                    
        #optimum number of clusters from validation index.
        sys.setrecursionlimit(10**8)
        file = file or filedialog.askopenfilename()
        data, col_groups = GB.readAndPreProcess(file=file,transform=transform,scale=scale,func='CC')
        metab_data = GB.readAndPreProcess(file=file,transform='None',scale='None',func='Raw')
        del(col_groups)

        # remember the last file used for ensemble so GUI can reuse it
        try:
            config.last_ensemble_file = file
        except Exception:
            pass

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
        # This call generates the ensemble clustergram AND writes cluster files.
        # recClustersPostVal now returns (opt_clusters, output_dir) for downstream use.
        opt_clusters, out_dir = GB.createEnsemDendrogramNew(
            coOcc, metab_data, data,
            norm=0, minMetabs=0, numClusts=len(parameters),
            link='average', dist='euclidean', func="ensemble", colMap=colMap,
        )
        return opt_clusters, out_dir

    def cooccHeatmap(num_clusters, coocc_file, data_file, cmap='coolwarm',
                     col_cluster=False, transform='None', scale='None'):
        '''
        Build clustermaps from EnsembleCoOcc.xlsx, user-chosen K, and a
        group-medians matrix (.xlsx): row 1 = headers (group names in sample
        columns), then one row per feature (same order as ensemble). If the
        medians file has one extra row (Group Medians spacer under the header),
        it is dropped automatically when N = CoOcc + 1.

        transform / scale - same pipeline as ensemble clustering (applied to
        the numeric group-median matrix for the heatmap only; Cluster*.xlsx
        exports stay raw medians).

        cmap - matplotlib colormap name (default coolwarm).
        col_cluster - if True, seaborn also hierarchically clusters columns
        (groups), adding a column dendrogram.

        Each heatmap is shown in a matplotlib window after saving the PDF;
        close the window to continue to the next plot, then the success dialog.
        '''
        logging.info(': User called CoOcc Heatmap.')

        try:
            num_clusters = int(num_clusters)
        except (TypeError, ValueError):
            messagebox.showerror(title='Error', message='Invalid number of clusters.')
            return

        if num_clusters < 1 or num_clusters > 10:
            messagebox.showerror(
                title='Error',
                message='Number of clusters must be between 1 and 10.',
            )
            return

        if not coocc_file or not os.path.isfile(coocc_file):
            messagebox.showerror(
                title='Error',
                message='Select a valid EnsembleCoOcc.xlsx file.',
            )
            return

        if not data_file or not os.path.isfile(data_file):
            messagebox.showerror(
                title='Error',
                message='Select the group medians matrix (.xlsx) for the heatmap.',
            )
            return

        cmap = (cmap or 'coolwarm').strip() or 'coolwarm'
        try:
            plt.get_cmap(cmap)
        except Exception:
            logging.warning(': Unknown colormap %r, using coolwarm.', cmap)
            cmap = 'coolwarm'

        try:
            ext = os.path.splitext(str(coocc_file))[1].lower()
            if ext in ('.xlsx', '.xlsm', '.xltx', '.xltm'):
                co_occ_df = pd.read_excel(
                    coocc_file, sheet_name=0, engine='openpyxl'
                )
            elif ext == '.xls':
                co_occ_df = pd.read_excel(
                    coocc_file, sheet_name=0, engine='xlrd'
                )
            else:
                co_occ_df = pd.read_excel(coocc_file, sheet_name=0)
            coOcc = np.asarray(co_occ_df.values, dtype=float)
        except Exception as e:
            logging.exception(': Failed to read co-occurrence file.')
            messagebox.showerror(
                title='Error',
                message=f'Could not read the co-occurrence matrix:\n{e}',
            )
            return

        if coOcc.ndim != 2 or coOcc.shape[0] != coOcc.shape[1]:
            messagebox.showerror(
                title='Error',
                message='Co-occurrence matrix must be square (N × N).',
            )
            return

        raw_data = GB.fileCheck(file=data_file)
        if raw_data is None:
            return

        n_co = int(coOcc.shape[0])
        base = raw_data.copy().reset_index(drop=True)

        # Medians .xlsx: header row = group names; all following rows should be
        # features. Group Medians output also inserts a spacer row (out.loc[0])
        # so the sheet often has N+1 rows vs N×N CoOcc from ensemble.
        variants = [base]
        if len(base) == n_co + 1:
            variants.append(base.iloc[1:].reset_index(drop=True))
            variants.append(base.iloc[:-1].reset_index(drop=True))

        metab_data = None
        trim_note = None
        for vi, cand in enumerate(variants):
            if len(cand) != n_co:
                continue
            try:
                mz_c, rt_c = GB.detectColumns(cand)
                sc_try = [c for c in cand.columns if c not in (mz_c, rt_c)]
                if not sc_try:
                    continue
                _ = cand[sc_try].astype(float)
                metab_data = cand
                if vi == 1:
                    trim_note = 'top'
                elif vi == 2:
                    trim_note = 'bottom'
                break
            except Exception:
                continue

        if metab_data is None:
            n_rows = len(base)
            messagebox.showerror(
                title='Error',
                message=(
                    f'Co-occurrence matrix is {n_co}×{n_co} (features) but the '
                    f'medians sheet has {n_rows} data row(s).\n\n'
                    'Those counts must match. Common cases:\n'
                    '• Group Medians saves an extra spacer row under the header '
                    '(first row is not a metabolite). If you see N+1 rows, '
                    'delete that row in Excel or regenerate medians.\n'
                    '• Use the same feature list as ensemble: CoOcc was built '
                    'from the raw matrix (sample names in row 1, group labels '
                    'in row 2, then features). Medians should be one row per '
                    'feature after the header row only.\n'
                    '• Wrong or mixed file (different dataset).'
                ),
            )
            return

        if trim_note:
            logging.info(
                ': CoOcc heatmap: dropped %s data row of medians file so row '
                'count matches CoOcc (N=%s).',
                trim_note,
                n_co,
            )

        mz_col, rt_col = GB.detectColumns(metab_data)
        sample_cols = [c for c in metab_data.columns if c not in (mz_col, rt_col)]
        feature_names = metab_data[mz_col].astype(str).tolist()

        sample_labels = GB.sample_labels_for_coocc_group_medians(
            metab_data, mz_col, rt_col,
        )
        if len(sample_labels) != len(sample_cols):
            sample_labels = [str(c) for c in sample_cols]

        k_data = metab_data[sample_cols].astype(float)
        k_data.index = pd.Index(feature_names, name=mz_col)
        k_data.columns = sample_labels

        transform = transform if transform else 'None'
        scale = scale if scale else 'None'
        vals = np.asarray(k_data.values, dtype=np.float64, order='C', copy=True)
        try:
            GB.transformations(vals, transform=transform, scale=scale, first='1')
        except Exception as e:
            logging.exception(': CoOcc heatmap transform/scale failed.')
            messagebox.showerror(
                title='Error',
                message=(
                    'Transform/scale failed on the median matrix.\n\n'
                    f'{e}\n\n'
                    'Try Transform/Scale = None, or match choices to ensemble.'
                ),
            )
            return
        k_vis = pd.DataFrame(
            vals, index=k_data.index, columns=k_data.columns,
        )

        link_mat = linkage(coOcc, method='average')
        labels = fcluster(link_mat, num_clusters, criterion='maxclust')

        out_dir = os.path.dirname(os.path.abspath(coocc_file))
        os.makedirs(out_dir, exist_ok=True)

        uniq = np.unique(labels)
        for i, lab in enumerate(uniq):
            ind = np.where(labels == lab)[0]
            cur = metab_data.iloc[ind].copy()
            sc = [c for c in cur.columns if c not in (mz_col, rt_col)]
            out = cur[[mz_col] + sc + [rt_col]].copy()
            out = out.rename(columns={mz_col: 'Identities', rt_col: 'rt_med'})
            cluster_path = os.path.join(out_dir, f'Cluster{i + 1}.xlsx')
            GB.safeToExcel(out, cluster_path, index=False)

        # Show x labels when few columns, or when column clustering is on
        # (group order changes and labels are needed to read the dendrogram).
        x_lbls = (len(sample_labels) <= 35) or col_cluster

        # Full feature IDs (same as Identities in Cluster*.xlsx); explicit list
        # so labels are not truncated like yticklabels='auto' can do.
        feature_y_labels = metab_data[mz_col].astype(str).tolist()

        try:
            g_feat = sns.clustermap(
                k_vis,
                figsize=(8, 8),
                row_linkage=link_mat,
                col_cluster=col_cluster,
                cmap=cmap,
                linecolor='black',
                yticklabels=feature_y_labels,
                xticklabels=x_lbls,
                cbar_pos=(0.01, 0.8, 0.025, 0.175),
            )
            g_feat.ax_heatmap.set_xlabel('')
            feat_pdf = os.path.join(out_dir, 'feature_labels_heatmap.pdf')
            g_feat.savefig(feat_pdf, dpi=600, transparent=True)
            plt.show()
            plt.close(g_feat.fig)
        except Exception as e:
            logging.warning(f': CoOcc feature-label heatmap failed ({e})')
            messagebox.showerror(
                title='Error',
                message=f'Feature-label heatmap failed:\n{e}',
            )
            return

        try:
            labels_y = [str(x) for x in labels]
            g_cl = sns.clustermap(
                k_vis,
                figsize=(8, 8),
                row_linkage=link_mat,
                col_cluster=col_cluster,
                cmap=cmap,
                linecolor='black',
                yticklabels=labels_y,
                xticklabels=x_lbls,
                cbar_pos=(0.01, 0.8, 0.025, 0.175),
            )
            g_cl.ax_heatmap.set_xlabel('')
            cl_pdf = os.path.join(out_dir, 'cluster_labels_heatmap.pdf')
            g_cl.savefig(cl_pdf, dpi=600, transparent=True)
            plt.show()
            plt.close(g_cl.fig)
        except Exception as e:
            logging.warning(f': CoOcc cluster-label heatmap failed ({e})')
            messagebox.showerror(
                title='Error',
                message=f'Cluster-label heatmap failed:\n{e}',
            )
            return

        extra = ''
        if trim_note == 'top':
            extra = '\n\n(Top data row of the medians file was skipped so it matches CoOcc.)'
        elif trim_note == 'bottom':
            extra = '\n\n(Bottom data row of the medians file was skipped so it matches CoOcc.)'
        messagebox.showinfo(
            title='Success',
            message=(
                f'Saved 2 heatmaps + Cluster1–Cluster{len(uniq)}.xlsx next to:\n'
                f'{coocc_file}\n\n'
                f'(Plots were also shown in figure windows.)'
                f'{extra}'
            ),
        )

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
    
    def validateMono(self, upperLimClust=10, transform='None', scale='None',
                     linkage='average', dissimilarity='euclidean', func='All'):
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
        #create a dictionary of validation metrics the user can call for analysis. 
        valIndex = {
            'Calinski-Harabasz': metrics.calinski_harabasz_score,
            'Silhouette':        metrics.silhouette_score,
            'Davies-Bouldin':    metrics.davies_bouldin_score,
        }

        #log that user called MST
        logging.info(': User called Cluster validation function.')

        filename = filedialog.askopenfilename()
        try:
            data, col_groups = GB.readAndPreProcess(file=filename, transform = transform, scale =scale, func='CC')
        except BaseException:
            logging.error(': Unable to proceed, due to file error!')
            messagebox.showerror(title='Error',message='Unable to proceed, try again or return to homepage!')
            return
        
        # evaluate k from 2..upperLimClust for all three metrics
        ks = list(range(2, upperLimClust + 1))

        # calculate pairwise distances once
        dist = pdist(data, dissimilarity)
        if linkage == 'ward':
            link_mat = ward(dist)
            cluster_func = lambda k: fcluster(link_mat, k, criterion='maxclust')
        else:
            def cluster_func(k):
                agg = AC(
                    n_clusters=k,
                    linkage=linkage,
                    metric='precomputed'
                ).fit(squareform(dist))
                return agg.labels_ + 1

        metrics_to_run = ['Silhouette', 'Davies-Bouldin', 'Calinski-Harabasz']
        scores = {m: [] for m in metrics_to_run}

        for k in ks:
            labels_ = cluster_func(k)
            for m in metrics_to_run:
                s = valIndex[m](data, labels_)
                scores[m].append(s)

        # plot 3-panel summary (no red guideline lines)
        fig, axes = plt.subplots(1, 3, figsize=(15, 4))
        for ax, m in zip(axes, metrics_to_run):
            y = np.asarray(scores[m])
            ax.plot(ks, y, marker='o')
            ax.set_title(m)
            ax.set_xlabel('Number of clusters (k)')
            #ax.grid(True, alpha=0.3)
            
        plt.figtext(1, 1, plt.figtext(0.5, 0.01, 'Optimal clusters: Silhouette & Calinski-Harabasz = max, Davies-Bouldin = min', ha='center', fontsize=9))
        plt.tight_layout()
        plt.savefig('MonoClust_Optimization.png', dpi=300, transparent=True)
        plt.show()
        logging.info(': Successfully completed mono-clustering optimization for all metrics.')
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

        # Keep true MST optimization behavior for the "MST Optimization" button.
        # That button calls this function with func='ensemble'.
        if func == 'ensemble':
            num_groups = data.shape[1]
            pairWise = squareform(pdist(data))
            mstInput = csr_matrix(pairWise)
            mstOut = minimum_spanning_tree(mstInput)
            mstOutInd = mstOut.nonzero()
            dataMST = np.zeros([data.shape[0] - 1, 3])
            for i in range(data.shape[0] - 1):
                dataMST[i, 0] = mstOutInd[0][i]
                dataMST[i, 1] = mstOutInd[1][i]
                dataMST[i, 2] = mstOut[int(dataMST[i, 0]), int(dataMST[i, 1])]

            mstOutDf = pd.DataFrame(dataMST, columns=['index1', 'index2', 'dist']).sort_values(by='dist')
            validationClusters = GB.clustConnect(dataMST, mstOutDf.to_numpy())
            argsMulti = [({0: validationClusters[i]}, data, num_groups)
                         for i in range(len(validationClusters))
                         if i >= max(0, len(validationClusters) - 101)]

            valIndex = [GB.Validate(*a) for a in argsMulti]
            valIndex = np.asarray(
                [np.reshape(np.asarray(v), (2, -1))[:, 0] for v in valIndex],
                dtype=float
            )
            K = valIndex[:, 1]
            y = valIndex[:, 0].astype(float)
            if y.shape[0] > 1:
                y[y.shape[0] - 1] = y[y.shape[0] - 2] * 2
            y = 1.0 / np.where(y == 0, np.nan, y)
            for i in range(y.shape[0]):
                if K[i] > 1 and np.isfinite(y[i]):
                    y[i] = (y[i] * K[i]) / (K[i] - 1)
            if np.all(np.isnan(y)):
                messagebox.showerror(title='Error', message='MST optimization could not determine an optimal K.')
                return
            optK = int(K[int(np.nanargmax(y))])
            messagebox.showinfo(message="MST suggests optimal number of clusters: " + str(optK))
            GB.valPlotting(valIndex, mstOutDf)
            return optK

        max_k = min(25, data.shape[0] - 1)
        
        if max_k < 2:
            messagebox.showerror(title='Error', message='Need at least 3 samples to run cluster validation.')
            return

        ks = list(range(2, max_k + 1))

        def _dunn_index(x, labels):
            uniq = np.unique(labels)
            if uniq.shape[0] < 2:
                return np.nan
            dm = squareform(pdist(x))
            max_intra = 0.0
            for c in uniq:
                idx = np.where(labels == c)[0]
                if idx.shape[0] < 2:
                    continue
                intra = dm[np.ix_(idx, idx)]
                max_intra = max(max_intra, float(np.max(intra)))
            if max_intra == 0:
                return np.nan
            min_inter = float("inf")
            for a in range(len(uniq)):
                ia = np.where(labels == uniq[a])[0]
                for b in range(a + 1, len(uniq)):
                    ib = np.where(labels == uniq[b])[0]
                    inter = dm[np.ix_(ia, ib)]
                    if inter.size > 0:
                        min_inter = min(min_inter, float(np.min(inter)))
            return np.nan if not np.isfinite(min_inter) else float(min_inter / max_intra)

        def _pbm_index(x, labels):
            uniq = np.unique(labels)
            k = uniq.shape[0]
            if k < 2:
                return np.nan
            center_all = np.mean(x, axis=0)
            e1 = float(np.sum(np.linalg.norm(x - center_all, axis=1)))
            ek = 0.0
            centers = []
            for c in uniq:
                cl = x[labels == c]
                ctr = np.mean(cl, axis=0)
                centers.append(ctr)
                ek += float(np.sum(np.linalg.norm(cl - ctr, axis=1)))
            if ek == 0:
                return np.nan
            centers = np.asarray(centers)
            dk = float(np.max(pdist(centers))) if centers.shape[0] > 1 else 0.0
            return float(((e1 * dk) / (ek * k)) ** 2)

        scores = {
            'k-means based': [],
            'Silhouette': [],
            'DBI': [],
            'Dunn': [],
            'PBM': [],
        }

        for k in ks:
            km = KMeans(n_clusters=k, random_state=42, n_init=10)
            labels = km.fit_predict(data)
            scores['k-means based'].append(float(km.inertia_))
            scores['Silhouette'].append(float(metrics.silhouette_score(data, labels)))
            scores['DBI'].append(float(metrics.davies_bouldin_score(data, labels)))
            scores['Dunn'].append(_dunn_index(data, labels))
            scores['PBM'].append(_pbm_index(data, labels))

        metric_meta = {
            'k-means based': {'title': 'KMeans WCSS (Elbow)', 'better': 'min'},
            'Silhouette': {'title': 'Silhouette', 'better': 'max'},
            'DBI': {'title': 'Davies-Bouldin', 'better': 'min'},
            'Dunn': {'title': 'Dunn', 'better': 'max'},
            'PBM': {'title': 'PBM', 'better': 'max'},
        }

        def _best_idx(arr, better):
            vals = np.asarray(arr, dtype=float)
            if np.all(np.isnan(vals)):
                return None
            return int(np.nanargmin(vals) if better == 'min' else np.nanargmax(vals))

        if func == 'All':
            fig, axes = plt.subplots(2, 3, figsize=(14, 9))
            axes = axes.flatten()
            ordered = ['k-means based', 'Silhouette', 'DBI', 'Dunn', 'PBM']
            for i, name in enumerate(ordered):
                y = np.asarray(scores[name], dtype=float)
                ax = axes[i]
                ax.plot(ks, y, marker='o')
                bi = _best_idx(y, metric_meta[name]['better'])
                #if bi is not None:
                    #ax.plot(ks[bi], y[bi], 'r.', markersize=12)
                ax.set_title(metric_meta[name]['title'])
                ax.set_xlabel('Number of clusters (k)')
                #ax.grid(True, alpha=0.3)
            axes[5].axis('off')
            axes[5].text(
                0.5, 0.5,
                "Optimal clusters guide:\n"
                "WCSS = elbow/slope change\n"
                "Silhouette, Dunn, PBM = max\n"
                "DBI = min",
                ha='center', va='center', fontsize=10, wrap=True
            )
            plt.tight_layout()
            plt.savefig('ValidationClusters_AllMetrics.png', dpi=300, transparent=True)
            plt.show()
            pd.DataFrame({'k': ks, **scores}).to_csv('ValidationClusters_AllMetrics.csv', index=False)
            return

        if func not in scores:
            messagebox.showerror(title='Error', message=f'Unknown validation metric: {func}')
            return

        y = np.asarray(scores[func], dtype=float)
        best_i = _best_idx(y, metric_meta[func]['better'])
        plt.figure(figsize=(6, 5), dpi=300)
        plt.plot(ks, y, marker='o')
        #if best_i is not None:
           # plt.plot(ks[best_i], y[best_i], 'r.', markersize=12)
           # plt.annotate(f'k={ks[best_i]}', (ks[best_i], y[best_i]))
        plt.title(metric_meta[func]['title'])
        plt.xlabel('Number of clusters (k)')
        plt.ylabel('Score')
        plt.figtext(1, 1, plt.figtext(0.5, 0.01, 'Optimal clusters: WCSS = elbow/slope change, Silhouette, Dunn, PBM = max, DBI = min,', ha='center', fontsize=9))
        #plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(f"ValidationClusters_{func.replace(' ', '_')}.png", dpi=300, transparent=True)
        plt.show()
        pd.DataFrame({'k': ks, func: y}).to_csv(f"ValidationClusters_{func.replace(' ', '_')}.csv", index=False)
        return

    def peaksToPathways(raw_file=None, clusters_dir=None):
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
        if raw_file is None:
            messagebox.showinfo(title='File selection', message="Please select the original data file submitted for clustering!!")
            filename = filedialog.askopenfilename()
        else:
            filename = raw_file

        dataRaw = GB.fileCheck(file = filename)
        if dataRaw is None:
            #log error and return function to ensure a soft closing of the class
            logging.error(': Error loading the reference Excel sheet.')
            return

        #ask user to select the directory containing the csv files from ensemble clustering output (currently only method available)
        if clusters_dir is None:
            messagebox.showinfo(title="Directory Selection",message="Please select the directory containing the ensemble clustering output files!")
            direct = filedialog.askdirectory()
        else:
            direct = clusters_dir
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
                ensemFiles.append(os.path.join(direct, curCheck))
            elif curCheck[0:5] == 'Clust':
                #append to ensemFiles
                ensemFiles.append(os.path.join(direct, curCheck))

        columnsHead = list(dataRaw.columns)
        columnFirst = columnsHead[0]
        for i in range(len(ensemFiles)):
            dataClust = np.ones((dataRaw.shape[0],3))
            dataClust[:,0] = dataRaw[columnFirst]
            dataClust[:,1] = dataRaw["rtmed"]
            #start the process of reading in and creating the ensemble output files. 
            try:
                dataCur = pd.read_excel(ensemFiles[i], engine='openpyxl')
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

    def heatmapAnalysis(linkFunc, distMet, cmap, norm, colOrder=[], transform='None', scale='None', file=None):
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
        file = file or filedialog.askopenfilename()
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
        This function works by asking the user to select a file with the a column of Cpd.Hits, from which it queries KEGG to determine the associated enzymes. 

        Input:
        xlsx file with a Cpd.Hits column

        Output: 
        xlsx workbook with multiple sheets that containing information on all of the clusters adn their corresponding compounds and enzymes. 
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

    def anovaHM(transform='Log transformation', scale='Auto Scaling', cMap='viridis', file=None):
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
        file = file or filedialog.askopenfilename()

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
        Input:
        A file of original data.
        A file of the t-test results. 

        Output:
        A single file with the confidence intervals for all of the metabolites considered. 
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

    def normalityCheck(transform='None', scale='None', file=None):
        '''
        Input:
        A file in which you would like to check the normality of before performing any downstream analyses. 

        Output:
        A plot showing the comparison of pre and post normalization results. 
        '''

        #have user input the wanted file
        file = file or filedialog.askopenfilename()

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
        Input:
        An output excel workbook from progenesis see example files for appropriate file type

        Output:

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

    def externalCriteria(comp='rand',trans ='None',scale='None'):
        '''
        This function is responsible for allowing users to compare single clustering solutions along with the appropriate optimization metric. 
        
        Input:
        File of data you want to cluster, and a file with the input metrics. 

        Output:
        Plots of the resulting comparisons between the clustering solutions. 

        '''
        # Defomomg ECCO parameters to try
        ENSEMBLE_PARAMS = [
        {'Distance': 'CNS', 'Linkage': 'ward', 'Correlation': 'spearman', 'Optimizer': 'CH'},
        {'Distance': 'CNS', 'Linkage': 'average', 'Correlation': 'spearman', 'Optimizer': 'SIL'},
        {'Distance': 'CS', 'Linkage': 'ward', 'Correlation': 'spearman', 'Optimizer': 'CH'},
        {'Distance': 'CS', 'Linkage': 'average', 'Correlation': 'spearman', 'Optimizer': 'SIL'},
        {'Distance': 'CS', 'Linkage': 'complete', 'Correlation': 'spearman', 'Optimizer': 'DBI'},
        {'Distance': 'CNS', 'Linkage': 'complete', 'Correlation': 'spearman', 'Optimizer': 'DBI'}
                            ]
        
        clusts = pd.DataFrame(ENSEMBLE_PARAMS)
        
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
    # External Comparison Logic
        comp_map = {
            'Rand-index': lambda: GB.randComp(best_labels, distLink),
            'Adjusted Rand-index': lambda: GB.adjRandComp(best_labels, distLink),
            'Normalized Mutual Info.': lambda: GB.mutualInfo(best_labels, distLink, 'norm'),
            'Adjusted Mutual Info.': lambda: GB.mutualInfo(best_labels, distLink, 'adj')
        }
    
        if comp in comp_map:
            comp_map[comp]()
            messagebox.showinfo(title='Completed', message=f'Successfully saved the output of the {comp} comparison.')
        else:
            messagebox.showinfo(message='Currently not supporting this comparison type.')

    def geneToPathway():
        '''Pathway enrichment via gseapy Enrichr + prerank.

        Input CSV: gene ID, log2 fold change, and raw p-value (column names
        matched case-insensitively, e.g. ECCO_Trial_Deglist.csv: Gene, log2FC,
        pvalue). Optional FDR/padj column is kept for export only; significance
        for up/down lists and volcano coloring uses raw p-value and a user-chosen
        cutoff (default 0.05).

        Enrichr gene_sets: GO BP/MF/CC 2021, KEGG_2021_Human (organism-agnostic
        block in this function), and Reactome_2022.
        '''
        try:
            import gseapy as gp
        except ImportError as e:
            messagebox.showerror(
                title="GSEApy import failed",
                message=(
                    "Pathway enrichment needs gseapy with the compiled gseapy.gse module "
                    "(from a PyPI wheel).\n\n"
                    "Remove any conda copy, then reinstall from pip:\n"
                    "  conda remove -y gseapy\n"
                    "  pip uninstall -y gseapy\n"
                    "  pip install --no-cache-dir \"gseapy>=1.1\"\n\n"
                    "If you are on Python 3.14+, use 3.12 or 3.13 until wheels exist.\n\n"
                    f"Error: {e}"
                ),
            )
            return

        sns.set(style="whitegrid", context="talk")

        # Prompt user
        messagebox.showinfo(
            message="Select CSV with Gene, log2FC, and raw p-value "
            "(e.g. pvalue, P.Value). FDR/padj optional. ECCO-style DEG lists supported."
        )
        file = filedialog.askopenfilename()
        if not file:
            return

        raw_df = pd.read_csv(file)
        lm = {str(c).strip().lower(): c for c in raw_df.columns}

        def _pick_col(*names):
            for n in names:
                k = n.lower()
                if k in lm:
                    return lm[k]
            return None

        gene_c = _pick_col("gene", "symbol")
        lfc_c = _pick_col("log2fc", "logfc", "log2foldchange", "log2_fold_change")
        p_c = _pick_col(
            "pvalue",
            "p.value",
            "p_value",
            "p val",
            "pval",
            "raw p",
            "raw_p",
            "pvalue.raw",
        )
        fdr_c = _pick_col("fdr", "padj", "qvalue", "adj.p.val", "adj_p_val", "qval")

        if gene_c is None or lfc_c is None or p_c is None:
            messagebox.showerror(
                message="CSV needs gene names, log2 fold change, and a raw p-value column "
                "(e.g. Gene, log2FC, pvalue).",
            )
            return

        df = raw_df.copy()
        df.rename(columns={gene_c: "Gene", lfc_c: "log2FC", p_c: "pvalue"}, inplace=True)
        df["Gene"] = df["Gene"].astype(str)
        df["log2FC"] = pd.to_numeric(df["log2FC"], errors="coerce")
        df["pvalue"] = pd.to_numeric(df["pvalue"], errors="coerce")
        if fdr_c is not None:
            df["FDR"] = pd.to_numeric(raw_df[fdr_c], errors="coerce")

        df = df.dropna(subset=["Gene", "log2FC", "pvalue"])

        p_cutoff = simpledialog.askfloat(
            "P-value threshold",
            "Significance cutoff (raw p-value):",
            initialvalue=0.05,
            minvalue=0.0,
            maxvalue=1.0,
        )
        if p_cutoff is None:
            return
        if p_cutoff <= 0:
            messagebox.showerror(message="P-value cutoff must be greater than 0.")
            return

        sig = df[df["pvalue"] < p_cutoff]

        threshold = simpledialog.askinteger(
            "Threshold",
            "What is your log2FC threshold?",
            initialvalue=0,
        )
        if threshold is None:
            return
        up = sig[sig["log2FC"] > threshold]
        down = sig[sig["log2FC"] < -threshold]

        up_genes = up["Gene"].astype(str).tolist()
        down_genes = down["Gene"].astype(str).tolist()

        def volcano_plot(df_plot, log2_thresh, p_thresh):
            eps = max(np.finfo(float).tiny, float(p_thresh) * 1e-6)
            ycol = "-log10(P-value)"
            df_v = df_plot.copy()
            df_v[ycol] = -np.log10(df_v["pvalue"].clip(lower=eps))
            df_v["Regulation"] = "NS"
            up_m = (df_v["pvalue"] < p_thresh) & (df_v["log2FC"] > log2_thresh)
            down_m = (df_v["pvalue"] < p_thresh) & (df_v["log2FC"] < -log2_thresh)
            df_v.loc[up_m, "Regulation"] = "Up"
            df_v.loc[down_m, "Regulation"] = "Down"

            plt.figure(figsize=(8, 6))
            palette = {"Up": "red", "Down": "blue", "NS": "gray"}
            for key, group in df_v.groupby("Regulation"):
                plt.scatter(
                    group["log2FC"],
                    group[ycol],
                    c=palette[key],
                    label=key,
                    alpha=0.7,
                    s=20,
                )

            top = df_v[df_v["pvalue"] < p_thresh].nlargest(10, ycol)
            for _, row in top.iterrows():
                plt.text(row["log2FC"], row[ycol], row["Gene"], fontsize=8)

            plt.axhline(-np.log10(p_thresh), linestyle="--", color="black")
            plt.axvline(log2_thresh, linestyle="--", color="black")
            plt.axvline(-log2_thresh, linestyle="--", color="black")

            plt.xlabel("log2 Fold Change")
            plt.ylabel(ycol)
            plt.title("Volcano Plot (P-value cutoff: {p_cutoff})")
            plt.legend()
            plt.tight_layout()
            plt.savefig("volcano_plot.png", dpi=300)
            plt.close()

        volcano_plot(df, threshold, p_cutoff)

        # GO libraries
        organism = "human"
        go_libraries = {
            "BP": "GO_Biological_Process_2021",
            "MF": "GO_Molecular_Function_2021",
            "CC": "GO_Cellular_Component_2021"
        }
        
        kegg_map = {"human": "KEGG_2021_Human"}

        kegg_lib = kegg_map.get(organism, "KEGG_2021_Human")
        
        reactome_lib = "Reactome_2022"

        geneset_lines = [
            f"Organism passed to Enrichr: {organism}",
            "",
            "Gene set libraries (Enrichr names):",
            f"  • {go_libraries['BP']} (GO Biological Process)",
            f"  • {go_libraries['MF']} (GO Molecular Function)",
            f"  • {go_libraries['CC']} (GO Cellular Component)",
            f"  • {kegg_lib}",
            f"  • {reactome_lib}",
            "",
            "GSEA prerank uses this same list (GO ×3, KEGG, Reactome).",
        ]
        if organism not in kegg_map:
            geneset_lines.append(
                f"(Unrecognized organism for KEGG: using {kegg_lib}; "
                "use 'human' or 'mouse' for matched KEGG libraries.)"
            )
        messagebox.showinfo(
            title="Gene set libraries",
            message="\n".join(geneset_lines),
        )

        # Run enrichment
        def run_enrichment(gene_list, library):
            if len(gene_list) == 0:
                return pd.DataFrame()

            enr = gp.enrichr(
                gene_list=gene_list,
                gene_sets=library,
                organism=organism,
                outdir=None,
                cutoff=0.05
            )
            return enr.results

        results = {}

        for name, lib in go_libraries.items():
            results[f"UP_GO_{name}"] = run_enrichment(up_genes, lib)
            results[f"DOWN_GO_{name}"] = run_enrichment(down_genes, lib)

        results["UP_KEGG"] = run_enrichment(up_genes, kegg_lib)
        results["DOWN_KEGG"] = run_enrichment(down_genes, kegg_lib)
        results["UP_Reactome"] = run_enrichment(up_genes, reactome_lib)
        results["DOWN_Reactome"] = run_enrichment(down_genes, reactome_lib)

        #GSEA (Ranked List)
        ranked = df[['Gene','log2FC']].sort_values(by='log2FC', ascending=False)
        ranked = ranked.set_index('Gene')

        gsea_res = gp.prerank(
            rnk=ranked,
            gene_sets=list(go_libraries.values()) + [kegg_lib, reactome_lib],
            outdir="GSEA_results",
            permutation_num=100
        )

        def plot_gsea_curves(pre_res, outdir="GSEA_plots"):
            os.makedirs(outdir, exist_ok=True)

            res2d = pre_res.res2d
            if res2d is None or getattr(res2d, "empty", True):
                return

            def _df_col(df, name):
                n = name.lower()
                for c in df.columns:
                    if str(c).lower() == n:
                        return c
                return None

            nes_c = _df_col(res2d, "nes")
            term_c = _df_col(res2d, "term")
            if nes_c is None or term_c is None:
                return

            res = res2d.sort_values(nes_c, ascending=False, na_position="last")
            high = list(res[term_c].head(3))
            low = list(res[term_c].tail(3))
            top_terms = []
            seen = set()
            for t in high + low:
                if t in seen:
                    continue
                seen.add(t)
                top_terms.append(t)

            def _safe_stem(label):
                s = str(label).replace(" ", "_")
                for ch in '<>:"/\\|?*\n\r\t':
                    s = s.replace(ch, "_")
                return (s[:180] if s else "pathway")

            res_by_term = pre_res.results

            def _gsd_get(gsd, *aliases):
                keymap = {str(k).lower(): k for k in gsd}
                for a in aliases:
                    k = keymap.get(a.lower())
                    if k is not None:
                        return gsd[k]
                raise KeyError(aliases)

            for term in top_terms:
                if term not in res_by_term:
                    continue
                ofn = os.path.join(outdir, f"{_safe_stem(term)}.png")
                try:
                    if hasattr(pre_res, "plot") and callable(pre_res.plot):
                        pre_res.plot(terms=term, ofname=ofn)
                    else:
                        gsd = res_by_term[term]
                        gp.gseaplot(
                            term=term,
                            hits=_gsd_get(gsd, "hits"),
                            nes=float(_gsd_get(gsd, "nes", "NES")),
                            pval=float(_gsd_get(gsd, "pval", "NOM p-val")),
                            fdr=float(_gsd_get(gsd, "fdr", "FDR q-val")),
                            RES=_gsd_get(gsd, "RES"),
                            rank_metric=pre_res.ranking,
                            ofname=ofn,
                        )
                except Exception:
                    continue

        plot_gsea_curves(gsea_res)

        #Reduce redundacy
        def reduce_redundancy(res, overlap_cutoff=0.5):

            if res.empty:
                return res

            res = res.sort_values('Adjusted P-value').copy()
            keep = []
            seen = []

            for _, row in res.iterrows():
                genes = set(row['Genes'].split(';'))

                redundant = False
                for s in seen:
                    overlap = len(genes & s) / len(genes | s)
                    if overlap > overlap_cutoff:
                        redundant = True
                        break

                if not redundant:
                    keep.append(row)
                    seen.append(genes)

            return pd.DataFrame(keep)

        # Plotting function
        def plot_results(res, name):
            if res.empty:
                return

            top = res.sort_values('Adjusted P-value').head(10).copy()

            # Extract gene counts
            top['Gene_Count'] = top['Overlap'].str.split('/').str[0].astype(int)

            # Barplot
            plt.figure(figsize=(8,6))
            sns.barplot(
                data=top,
                y='Term',
                x=-np.log10(top['Adjusted P-value']),
                color='steelblue'
            )
            plt.xlabel("-log10(FDR)")
            plt.ylabel("")
            plt.title(name)

            plt.tight_layout()
            plt.savefig(f"{name}_barplot.png", dpi=300)
            plt.close()

            # Dotplot
            plt.figure(figsize=(8,6))

            scatter = plt.scatter(
                -np.log10(top['Adjusted P-value']),
                top['Term'],
                s=top['Gene_Count'] * 25,
                c=top['Adjusted P-value'],
                cmap='viridis_r',
                edgecolor='black'
            )

            plt.xlabel("-log10(FDR)")
            plt.ylabel("")
            plt.title(name)

            # Colorbar (FDR)
            cbar = plt.colorbar(scatter)
            cbar.set_label("FDR")

            # Size legend
            for size in [5, 10, 20]:
                plt.scatter([], [], s=size*25, c='gray', label=str(size))
            plt.legend(title="Gene Count", loc='lower right')

            plt.tight_layout()
            plt.savefig(f"{name}_dotplot.png", dpi=300)
            plt.close()

        # Enrichment Map
        def enrichment_map(res, name):

            if res.empty:
                return

            top = res.sort_values('Adjusted P-value').head(20).copy()

            # Convert gene lists
            top['Gene_List'] = top['Genes'].apply(lambda x: set(x.split(';')))

            G = nx.Graph()

            for i in range(len(top)):
                G.add_node(top.iloc[i]['Term'])

            # Add edges based on overlap
            for i in range(len(top)):
                for j in range(i+1, len(top)):
                    overlap = len(top.iloc[i]['Gene_List'] & top.iloc[j]['Gene_List'])
                    if overlap >= 3:
                        G.add_edge(top.iloc[i]['Term'], top.iloc[j]['Term'], weight=overlap)
    
            plt.figure(figsize=(8,8))
            pos = nx.spring_layout(G, seed=42)

            nx.draw(
                G, pos,
                with_labels=True,
                node_size=500,
                font_size=6
            )

            plt.title(f"{name} Enrichment Map")
            plt.savefig(f"{name}_network.png", dpi=300)
            plt.close()
            
        # Generate plots
        for name, res in results.items():
            res = reduce_redundancy(res)

            plot_results(res, name)
            enrichment_map(res, name)
            
        # Identify unused genes (not recognized by enrichment tools)
        used_genes = set()
        for res in results.values():
            if not res.empty:
                used_genes.update(
                    sum([g.split(';') for g in res['Genes']], [])
                )

        unused_genes = set(df['Gene'].astype(str)) - used_genes

        # Write to Excel
        with pd.ExcelWriter("Gene_Enrichment_Results.xlsx") as writer:
            df.to_excel(writer, sheet_name="All Genes", index=False)
            up.to_excel(writer, sheet_name="Upregulated", index=False)
            down.to_excel(writer, sheet_name="Downregulated", index=False)

            for name, res in results.items():
                res.to_excel(writer, sheet_name=name[:30], index=False)

            pd.DataFrame(list(unused_genes), columns=["Unused Genes"])\
                .to_excel(writer, sheet_name="Unused Genes", index=False)

        messagebox.showinfo(
            message=f"Analysis complete!\n\nUnused genes: {len(unused_genes)}"
        )

    def RNASeqDeg():
        """
        DESeq2 on a count-matrix CSV; save a DEG list like ECCO_Trial_Deglist.csv.
        Expects a header row (gene id + sample columns), then either:
        - **ECCO-style**: first data row = group label per sample (e.g. condition1...condition1, condition2...condition2),
          following rows = genes × raw counts; or
        - **Simple matrix**: first data row starts with a gene name and numeric counts
          (no separate group row) — not supported; an error explains the Condition1vCondition2 format.

        Output columns: Gene, log2FC, pvalue, FDR (Wald test, BH-adjusted ``padj``).
        """
        if DefaultInference is None:
            messagebox.showerror(
                title="PyDESeq2",
                message="This workflow requires pydeseq2 with DefaultInference "
                "(upgrade pydeseq2 to a recent release, e.g. >= 0.4).",
            )
            return

        messagebox.showinfo(
            message="Select a count matrix CSV (e.g. ECCO_Count_Matrix.csv): "
            "header = gene column + sample IDs; row 2 = one group label per sample "
            "(e.g. Condition1, Condition2); remaining rows = genes and raw integer counts.",
        )
        path = filedialog.askopenfilename(filetypes=[("CSV", "*.csv"), ("All", "*.*")])
        if not path:
            return

        raw = pd.read_csv(path, header=0, dtype=str, low_memory=False)
        if raw.shape[1] < 3:
            messagebox.showerror(message="Need at least one gene column and two sample columns.")
            return

        sample_cols = [str(c) for c in raw.columns[1:]]

        def _cell_float(x):
            try:
                return float(str(x).replace(",", "").strip())
            except ValueError:
                return np.nan

        group_row = None
        data_start = 0
        if raw.shape[0] >= 2:
            top0 = _cell_float(raw.iloc[0, 1])
            top1 = _cell_float(raw.iloc[1, 1])
            if np.isnan(top0) and not np.isnan(top1):
                group_row = [str(x).strip() for x in raw.iloc[0, 1:].tolist()]
                data_start = 1

        if group_row is None:
            messagebox.showerror(
                message="Could not find a group-label row after the header.\n\n"
                "Use a format like ECCO_Count_Matrix.csv: row 1 = gene_id, samples…; "
                "row 2 = condition per sample (e.g. Condition1, Condition1, …, Condition2, Condition2, …); "
                "row 3+ = genes and counts.",
            )
            return

        if len(group_row) != len(sample_cols):
            messagebox.showerror(
                message="Number of group labels does not match the number of sample columns.",
            )
            return

        body = raw.iloc[data_start:].copy()
        body.index = body.iloc[:, 0].astype(str).str.strip()
        body = body.iloc[:, 1:]
        body.columns = sample_cols
        body = body.apply(pd.to_numeric, errors="coerce").fillna(0.0)
        body = body.astype(np.int64)
        body = body[body.index.astype(str).str.len() > 0]
        if body.index.duplicated().any():
            body = body[~body.index.duplicated(keep="first")]

        # PyDESeq2: samples × genes
        counts = body.T
        counts.index = counts.index.astype(str)

        metadata = pd.DataFrame({"group": group_row}, index=counts.index)
        metadata["group"] = metadata["group"].astype(str).str.strip()
        _g = metadata["group"]
        if _g.isna().any() or (_g == "").any():
            messagebox.showerror(message="Missing or empty group label for some samples.")
            return

        levels = sorted(set(metadata["group"].tolist()))
        if len(levels) != 2:
            messagebox.showerror(
                message=f"DESeq2 comparison needs exactly two groups; found: {levels}",
            )
            return

        if set(levels) == {"Condition1", "Condition2"}:
            test_g, ref_g = "Condition2", "Condition1"
        else:
            test_g, ref_g = levels[1], levels[0]

        test_g = simpledialog.askstring(
            "Contrast",
            "Log2 fold-change = log2(test / reference).\n"
            "Enter test group (numerator, e.g. Condition2):",
            initialvalue=test_g,
        )
        if test_g is None:
            return
        test_g = test_g.strip()
        ref_g = simpledialog.askstring(
            "Contrast",
            "Enter reference group (denominator, e.g. Condition1):",
            initialvalue=ref_g,
        )
        if ref_g is None:
            return
        ref_g = ref_g.strip()
        if test_g == ref_g:
            messagebox.showerror(message="Test and reference groups must differ.")
            return

        inference = DefaultInference(n_cpus=None)

        try:
            dds = DeseqDataSet(
                counts=counts,
                metadata=metadata,
                design="~group",
                refit_cooks=True,
                inference=inference,
                quiet=True,
            )
            dds.deseq2()
            stats = DeseqStats(
                dds,
                contrast=["group", test_g, ref_g],
                inference=inference,
                quiet=True,
            )
            stats.summary()
        except Exception as exc:
            messagebox.showerror(title="DESeq2 failed", message=str(exc))
            return

        res = stats.results_df
        out = pd.DataFrame(
            {
                "Gene": res.index.astype(str),
                "log2FC": res["log2FoldChange"],
                "pvalue": res["pvalue"],
                "FDR": res["padj"],
            }
        )

        default_out = os.path.splitext(path)[0] + "_deglist.csv"
        out_path = filedialog.asksaveasfilename(
            title="Save DEG list (ECCO-style)",
            defaultextension=".csv",
            filetypes=[("CSV", "*.csv"), ("All", "*.*")],
            initialfile=os.path.basename(default_out),
            initialdir=os.path.dirname(path) or ".",
        )
        if not out_path:
            return

        out.to_csv(out_path, index=False)
        messagebox.showinfo(
            title="Done",
            message=f"Saved {out.shape[0]} genes to:\n{out_path}\n\n"
            f"Contrast: {test_g} vs {ref_g} (reference).",
        )