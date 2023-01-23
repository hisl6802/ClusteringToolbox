#default number of threads for 
numThreads = 1
curUser = "Someone"
colorNum = 0
#colorlist used across many of the different functionalities of the GUI
colorList = ('viridis', 'plasma', 'inferno', 'magma', 'cividis','Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
                      'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
                      'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn','Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
                      'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
                      'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn','PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu', 'RdYlBu',
                      'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic','twilight', 'twilight_shifted', 'hsv')
#transform list, scale list, linkage list, distance metric list, normarilization list
transformList = ('None','Log transformation', 'Square root transformation', 'Cube root transformation')
scaleList = ('None', 'Mean centering', 'Auto Scaling', 'Pareto Scaling', 'Range Scaling')
linkageList = ('single','ward','complete','average')
distList = ('euclidean','seuclidean','sqeuclidean','cosine','chebyshev','correlation','canberra','braycurtis','minkowski','cityblock')
normList = ('Normalize','Do not')


#parameters for bootstrapping
numReSamp = 0
numPerSamp = 0

#current parameters for function being called
curTrans = 'None'
curScale = 'None'

