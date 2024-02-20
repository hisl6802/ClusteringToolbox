from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.chrome.options import Options
import time
import os
import zipfile



#set up the browser runner, update the directory, get all the csv files within the folder of interest
file = '/Users/bradyhislop/Documents/GitHub/ClusteringGUI/Testing/ExampleP2P.csv'

curDir = os.getcwd()
driver = webdriver.Chrome()



## Navigate to the MetaboAnalyst Website for the one-factor stats analysis
driver.get('https://www.metaboanalyst.ca/MetaboAnalyst/upload/PeakUploadView.xhtml')


#determine the mode used for analysis. 
lcMode = 'Negative'
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
        print('Nope')


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
    print('Didn''t proceed')

pval = 0.05
db = 'Mouse (KEGG)'
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
    print('Not working yet')


try:
    #making sure the download tab is available
    element = WebDriverWait(driver,10).until(
        EC.presence_of_element_located((By.XPATH,'//*[@id="treeForm:j_idt77:4"]/div'))
    )
    driver.find_element(By.XPATH,'//*[@id="treeForm:j_idt77:4"]/div').click()

except:
    print('Nope')


try:
    time.sleep(2)
    element = WebDriverWait(driver,10).until(
        EC.presence_of_element_located((By.XPATH,'//*[@id="ac:form1:j_idt20_data"]/tr[1]/td[1]/a'))
    )
    driver.find_element(By.XPATH,'//*[@id="ac:form1:j_idt20_data"]/tr[1]/td[1]/a').click()

except:
    print('Nope')

