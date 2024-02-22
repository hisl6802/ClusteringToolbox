import seaborn as sns
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list
import matplotlib.pyplot as plt
import GuiBackground as GB


#put in transform scale variables for testing
transform = 'Log transformation'
scale = 'Auto Scaling'
link = 'ward'
dist = 'euclidean'
file = '/Users/bradyhislop/Documents/GitHub/ClusteringGUI/ExampleFiles/HeatmapAnalysis/HM_analysis.xlsx'

#log that the user called the Create Clustergram function
print(': User called the Create Clustergram Function.')
#check that the file the user selects is appropriate

#send the data off to the readAndPreProcess function for analysis. 
data, col_groups = GB.readAndPreProcess(file=file,transform=transform,scale=scale,func="CC")

del(col_groups)
#create messagebox explaining to users how they need to select clusters.
print('Select clusters of interest, cluster and peak to pathway files will be automatically generated!')

#Create the appropriate plt figure to allow for the comparison of linkage functions
fig, axes = plt.subplots(1,1,figsize=(8,8))

#find the linkages
linkageOne = linkage(data,link,metric=dist)

if len(linkageOne[:,2]) == len(np.unique(linkageOne[:,2])):
    print('No need to jitter data!')

else:
    print(': Matching distance need to jitter distances')
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

sns.heatmap(dataFinal, annot=True,fmt='float', cmap="viridis",yticklabels=False)
plt.show()