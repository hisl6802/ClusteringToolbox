
import pandas as pd
import GuiBackground as GB
func = 'CC'
#ask the user for the input excel workbook needs to contain two sheets
file = "/Users/bradyhislop/Downloads/Example_AHM.xlsx"

#open sheet 0 - containing the original data
#open sheet 1 - containing all ANOVA outcomes or the pre-truncated ANOVA results
try:
    anovaRes = pd.read_excel(file,sheet_name=1)
except:
    print('Unable to open file!')



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
dataOrig = GB.readAndPreProcess(file=file,transform='Log transformation',scale='Auto Scaling',func='ANHM')

dataUpdated = dataOrig[metabNames.isin(anovaRes[anovaResC[0]])]
# del(dataOrig)

#send data to be transformed
dataUpdated = pd.DataFrame(dataUpdated,index = anovaRes[anovaResC[0]],columns=colHeaders)
# dataUpdated.to_excel('RawTopANOVA.xlsx')
# g=sns.clustermap(dataUpdated,cMap=cMap)
# plt.show()
