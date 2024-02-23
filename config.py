import inspect
import os
import GuiBackground as GB



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
#transform list, scale list, linkage list, distance metric list, normarilization list, optimization metrics, and metrics for cluster comparison
transformList = ('None','Log transformation', 'Square root transformation', 'Cube root transformation')
scaleList = ('None', 'Mean centering', 'Auto Scaling', 'Pareto Scaling', 'Range Scaling')
linkageList = ('single','ward','complete','average')
distList = ('euclidean','seuclidean','sqeuclidean','cosine','chebyshev','correlation','canberra','braycurtis','minkowski','cityblock')
normList = ('Normalize','Do not')
optList = ('Silhouette','Davies-Bouldin','Calinski-Harabasz')
metrics = ('Rand-index', 'Adjusted Rand-index', 'Mono-clustering comparison','Normalized Mutual Info.','Adjusted Mutual Info.')

###mummichog DB's
mummidbs = ('Human (MFN)','Mouse (KEGG)','Human (BioCyc)','Human (KEGG)','Mouse (BioCyc)','Rat (KEGG)','Cow (KEGG)')
modes = ('Positive','Negative')
pval = (0.05)


###Look-up list location (Glycan's)
module_file_path = inspect.getfile(GB)
directory = os.path.dirname(module_file_path)
keggGlycanLoc = (directory+'/KEGG_Compound_Glycans.xlsx')

#list of the options for MetaboAnalyst
typeAnalysis = ('Uni', 'Multi')
varianceFilters = ('IQR', 'SD', 'MAD','RSD','MAD_m')
sampNorm = ('None','Sum','Median','Quantile')

#keyword arguments for grid, simply add more as need. 
grid_kwargs = {
    0: {'column': 0, 'row': 0, 'columnspan': 4},
    1: {'column': 1, 'row': 1, 'sticky': 'nsew'},
    2: {'column': 2, 'row': 1, 'sticky': 'nsew'},
    3: {'column': 3, 'row': 1, 'sticky': 'nsew'},
    4: {'column': 1, 'row': 2, 'sticky': 'nsew'},
    5: {'column': 2, 'row': 2, 'sticky': 'nsew'},
    6: {'column': 3, 'row': 2, 'sticky': 'nsew'},
    7: {'column': 1, 'row': 3, 'sticky': 'nsew'},
    8: {'column': 2, 'row': 3, 'sticky': 'nsew'},
    9: {'column': 3, 'row': 3, 'sticky': 'nsew'},
    10: {'column': 1, 'row': 4, 'sticky': 'nsew'},
    11: {'column': 2, 'row': 4, 'sticky': 'nsew'},
    12: {'column': 3, 'row': 4, 'sticky': 'nsew'},
    13: {'column': 1, 'row': 5, 'sticky': 'nsew'},
    14: {'column': 2, 'row': 5, 'sticky': 'nsew'},
    15: {'column': 3, 'row': 5, 'sticky': 'nsew'},
    16: {'column': 1, 'row': 6, 'sticky': 'nsew'},
    17: {'column': 2, 'row': 6, 'sticky': 'nsew'},
    18: {'column': 3, 'row': 6, 'sticky': 'nsew'},
    19: {'column': 1, 'row': 7, 'sticky': 'nsew'},
    20: {'column': 2, 'row': 7, 'sticky': 'nsew'},
    21: {'column': 3, 'row': 7, 'sticky': 'nsew'},
    22: {'column': 1, 'row': 8, 'sticky': 'nsew'},
    23: {'column': 2, 'row': 8, 'sticky': 'nsew'},
}

#parameters for bootstrapping
numReSamp = 0
numPerSamp = 0

#current parameters for function being called
curTrans = 'None'
curScale = 'None'

