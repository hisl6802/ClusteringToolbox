import pandas as pd
import seaborn as sns
import numpy as np

from matplotlib import pyplot as plt

#read in data
data = pd.read_excel('ControlYInjuredY.xlsx')

#drop labels
data = data.drop(0,axis=0)
data.head()

#drop mz
data = data.drop('Unnamed: 0',axis=1)
data.head()
#drop rt
data = data.drop('rtmed',axis=1)
data.head()

data = np.reshape(data,(data.size,1))

sns.kdeplot(data)

plt.show()