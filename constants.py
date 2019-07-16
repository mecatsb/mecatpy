from os.path import join
from os.path import normcase

PREFIX = '..'
MECAT_BASE = join(PREFIX, "Models")                                 # MECAT base folder (for models)
DATAPATH_KEGG = join(PREFIX, "Data", "2019")                        # folder with KEGG data (put raw data here)
DATAPATH_KEGG_BASE = join(PREFIX, "MECAT", "Data")                  # base folder for KEGG data
DATAPATH_THERMODYNAMICS = normcase(join(PREFIX, 'Thermodynamics'))  # folder with thermodynamics data
DATAPATH_STATISTICS_KEGG = join(DATAPATH_KEGG, "Statistik")         # folder for KEGG statistics
FOR_MATLAB = True                                                   # indices start with 1 instead of 0
