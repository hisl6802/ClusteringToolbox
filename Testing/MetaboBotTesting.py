from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.chrome.options import Options
import time
import os
import zipfile

#file name after analysis
name = 'Testing1'
file = '/Users/bradyhislop/Documents/GitHub/ClusteringGUI/ExampleFiles/Uni.csv'

#make it run headerless
# chrome_options = Options()
# chrome_options.add_argument("--headless") 
#get the appropriate driver for chrome, and go to statistical analysis homepage
driver = webdriver.Chrome()#options=chrome_options
driver.get('https://www.metaboanalyst.ca/MetaboAnalyst/upload/StatUploadView.xhtml')


#wait the homepage loads, if it doesn't load in 10 seconds quit.  
try:
    #wait for element to show up
    element = WebDriverWait(driver,10).until(
        EC.presence_of_element_located((By.XPATH,'//*[@id="j_idt12:j_idt20"]/div[3]/div'))
    )
    #with element detected click it.
    driver.find_element(By.XPATH,'//*[@id="j_idt12:j_idt20"]/div[3]/div').click()
    driver.find_element(By.XPATH,'//*[@id="j_idt12:j_idt26_input"]').send_keys(file)
    driver.find_element(By.XPATH,'//*[@id="j_idt12:j_idt27"]').click()
except:
    driver.quit()
    print('I quit on homepage')

#wait for the data check to complete  
try:
    #wait for element to show up
    element = WebDriverWait(driver,10).until(
        EC.presence_of_element_located((By.XPATH,'//*[@id="form1:j_idt18"]'))
    )
    #proceed with data analysis
    driver.find_element(By.XPATH,'//*[@id="form1:j_idt18"]').click()
except:
    driver.quit()
    print('I quit on data processing page')

#check for the correct button then move on.
varianceFilter ='SD'
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


    #always need to submit and proceed. 
    driver.find_element(By.XPATH,'//*[@id="j_idt14:j_idt41"]').click()
    driver.find_element(By.XPATH,'//*[@id="j_idt14:j_idt42"]').click()

except:
    driver.quit()
    #add message box here if you end up here.
    print('I quit on the standard deviation page')



#get the data normalized
sampNorm = ('None','Sum','Median','Quantile')
sampleNorm = sampNorm[0]
dataTrans = ('None','Log10','Sqrt','Cube')
trans = dataTrans[0]
scale = 'Mean-Center'
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
    if trans == 'Log10':
        #click on log transformation
        driver.find_element(By.XPATH,'//*[@id="form1:j_idt74"]/div[2]/span').click()

    elif trans == 'Sqrt':
        #click on square-root tranformation
        driver.find_element(By.XPATH,'//*[@id="form1:j_idt78"]/div[2]/span').click()

    elif trans == 'Cube':
        #click on the cube-root transformation
        driver.find_element(By.XPATH,'//*[@id="form1:j_idt82"]/div[2]/span').click()

    #select the dataScale = ('None','Mean-Center','Auto-scale','Pareto-scale','Range-scale')correct data scaling
    if scale == 'Mean-Center':
        #click on mean centering
        driver.find_element(By.XPATH,'//*[@id="form1:j_idt91"]/div[2]/span').click()

    elif scale == 'Auto-scale':
        #click on auto-scaling
        driver.find_element(By.XPATH,'//*[@id="form1:j_idt94"]/div[2]/span').click()

    elif scale == 'Pareto-scale':
        #click on paret-scaling
        driver.find_element(By.XPATH,'//*[@id="form1:j_idt97"]/div[2]/span').click()

    elif scale =='Range-scale':
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
    print('I quit on the data normalization page') 

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
    print(' I quit on the preprocessing page')


####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################
######################################################################## MULTIVARIATE ANALYSIS CODE ################################################################################

# #ANOVA
# try:
#     #the element is present but not clickable. 
#     #wait for element to show up
#     element = WebDriverWait(driver,15).until(
#         EC.presence_of_element_located((By.XPATH,'//*[@id="j_idt12"]/table/tbody/tr[2]/td/table/tbody/tr[2]/td/table/tbody/tr[2]/td/a'))
#     )
#     driver.find_element(By.XPATH,'//*[@id="j_idt12"]/table/tbody/tr[2]/td/table/tbody/tr[2]/td/table/tbody/tr[2]/td/a').click()  

# except:
#     driver.quit()
#     #add appropriate message
#     print(' I quit on the preprocessing page')


# #principal components analysis.
# try:
#     #wait for element
#     element = WebDriverWait(driver,10).until(
#         EC.presence_of_element_located((By.XPATH,'//*[@id="treeForm:j_idt92:3_7"]/div'))
#     )
#     driver.find_element(By.XPATH,'//*[@id="treeForm:j_idt92:3_7"]/div').click()
#     #wait for element
#     element = WebDriverWait(driver,10).until(
#         EC.presence_of_element_located((By.XPATH,'//*[@id="ac"]/ul/li[3]/a'))
#     )
#     driver.find_element(By.XPATH,'//*[@id="ac"]/ul/li[3]/a').click()
# except:
#     driver.quit()
#     print('I quit on the principal components analysis page')


# #PLS-DA
# try:
    
#     #wait for element
#     element = WebDriverWait(driver,10).until(
#         EC.presence_of_element_located((By.XPATH,'//*[@id="treeForm:j_idt167:3_8"]/div'))
#     )
#     driver.find_element(By.XPATH,'//*[@id="treeForm:j_idt167:3_8"]/div').click()

# except:
#     driver.quit()
#     #add needed message.
#     print('I quit on initial PLS-DA')

# #PLS-DA extra's
# try:
#     time.sleep(2)
#     #wait for element
#     element = WebDriverWait(driver,10).until(
#         EC.presence_of_element_located((By.XPATH,'//*[@id="ac"]/ul/li[2]/a'))
#     )
    
#     driver.find_element(By.XPATH,'//*[@id="ac"]/ul/li[2]/a').click()
#     #wait for element
#     element = WebDriverWait(driver,10).until(
#         EC.presence_of_element_located((By.XPATH,'//*[@id="ac"]/ul/li[4]/a'))
#     )
#     driver.find_element(By.XPATH,'//*[@id="ac"]/ul/li[4]/a').click()

# except:
#     driver.quit()
#     #add appropriate message
#     print('I quit on PLS-DA extras')

# #create a dendrogram 
# try:
#     element = WebDriverWait(driver,10).until(
#         EC.presence_of_element_located((By.XPATH,'//*[@id="treeForm:j_idt220:3_13"]/div'))
#     )
#     driver.find_element(By.XPATH,'//*[@id="treeForm:j_idt220:3_13"]/div').click()

# except:
#     driver.quit()
#     #add appropriate message
#     print('I quit at the dendrogram')

# #create a heatmap
# try:
#     element = WebDriverWait(driver,10).until(
#         EC.presence_of_element_located((By.XPATH,'//*[@id="treeForm:j_idt35:3_14"]/div'))
#     )
#     driver.find_element(By.XPATH,'//*[@id="treeForm:j_idt35:3_14"]/div').click()

# except:
#     driver.quit()
#     #add appropriate message
#     print('I quit at the heatmap')


####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################
######################################################################## UNIVARIATE ANALYSIS CODE ##################################################################################



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
    print('I quit on the fold-change selection page')

#peform a t-test
try:
    element = WebDriverWait(driver,10).until(
        EC.presence_of_element_located((By.XPATH,'//*[@id="treeForm:j_idt78:3_1"]/div'))
    )
    driver.find_element(By.XPATH,'//*[@id="treeForm:j_idt78:3_1"]/div').click()
except:
    driver.quit()
    print('I quit on the t-test page')


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
    print('I quit on the volcano plot analyses page')



#principal components analysis.
try:
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
    print('I quit on the principal components analysis page')


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
    print('I quit on initial PLS-DA')

#PLS-DA extra's
try:
    time.sleep(3)
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
    print('I quit on PLS-DA extras')

#create a dendrogram 
try:
    element = WebDriverWait(driver,10).until(
        EC.presence_of_element_located((By.XPATH,'//*[@id="treeForm:j_idt220:3_13"]/div'))
    )
    driver.find_element(By.XPATH,'//*[@id="treeForm:j_idt220:3_13"]/div').click()

except:
    driver.quit()
    #add appropriate message
    print('I quit at the dendrogram')

#create a heatmap
try:
    element = WebDriverWait(driver,10).until(
        EC.presence_of_element_located((By.XPATH,'//*[@id="treeForm:j_idt35:3_14"]/div'))
    )
    driver.find_element(By.XPATH,'//*[@id="treeForm:j_idt35:3_14"]/div').click()

except:
    driver.quit()
    #add appropriate message
    print('I quit at the heatmap')

### go to the download page
try:
    element = WebDriverWait(driver,10).until(
        EC.presence_of_element_located((By.XPATH,'//*[@id="treeForm:j_idt127:4"]/div'))
    )
    driver.find_element(By.XPATH,'//*[@id="treeForm:j_idt127:4"]/div').click()

except:
    #driver.quit()
    #add appropriate message
    print('I didn''t make it to the download page' )


#download results
try:
    print('Here I am, in the downloading phase')
    element = WebDriverWait(driver,10).until(
        EC.presence_of_element_located((By.XPATH,'//*[@id="ac:form1:j_idt20_data"]/tr[1]/td[1]/a'))
    )
    time.sleep(2)
    driver.find_element(By.XPATH,'//*[@id="ac:form1:j_idt20_data"]/tr[1]/td[1]/a').click()
    print('I clicked!')
    
except:
    driver.quit()
    #Add appropriate message
    print('I didn''t download the data. ')

time.sleep(10)
# #rename zip file, but first check that it has downloaded.
# rename = name+'.zip'
# basepath = os.path.expanduser('~')
# basepath +='/Downloads/Download.zip'

# print(basepath)
# start = time.time()
# downloadTime = 0
# downloaded = False
# print('about to start')
# while downloadTime < 60:
#     #checking for download.
#     downloadTime = time.time() - start
#     if os.path.exists(basepath):
#         #file has downloaded
#         downloaded = True
#         break

# if downloaded:
#     #rename the file and unzip into the appropriate folder
#     os.rename(basepath,rename)
#     directory = os.getcwd()
#     zip_dir = directory + '/'+ name
#     with zipfile.ZipFile(rename,'r') as zip_ref:
#         zip_ref.extractall(zip_dir)
    
#     print('complete appropriately')