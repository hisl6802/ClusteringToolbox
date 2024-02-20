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
file = '/Users/bradyhislop/Documents/GitHub/ClusteringGUI/ExampleFiles/Multi.csv'

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
try:
    #wait for element to show up
    element = WebDriverWait(driver,10).until(
        EC.presence_of_element_located((By.XPATH,'//*[@id="j_idt14:j_idt24"]/div[2]/div/div/div[2]/span'))
    )
    #standard deviation filtering into pre-processing
    driver.find_element(By.XPATH,'//*[@id="j_idt14:j_idt24"]/div[2]/div/div/div[2]/span').click()
    driver.find_element(By.XPATH,'//*[@id="j_idt14:j_idt41"]').click()
    driver.find_element(By.XPATH,'//*[@id="j_idt14:j_idt42"]').click()

except:
    driver.quit()
    #add message box here if you end up here.
    print('I quit on the standard deviation page')

#get the data normalized
try:
    #wait for element to show up
    element = WebDriverWait(driver,10).until(
        EC.presence_of_element_located((By.XPATH,'//*[@id="form1:j_idt74"]/div[2]/span'))
    )
    #click on log-transformation
    driver.find_element(By.XPATH,'//*[@id="form1:j_idt74"]/div[2]/span').click()

    #wait for element to show up
    element = WebDriverWait(driver,10).until(
        EC.presence_of_element_located((By.XPATH,'//*[@id="form1:j_idt94"]/div[2]/span'))
    )
    driver.find_element(By.XPATH,'//*[@id="form1:j_idt94"]/div[2]/span').click()
    #wait for element to show up
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
    print(' I quit on the preprocessing page')


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
