import zipfile 
from urllib import urlopen 
from StringIO import StringIO 

import pandas as pd 
import numpy as np 

def get_data_pwt(path=None, version=71, date='11302012', extract=True, panel=True):
    """Loads the Penn World Tables (PWT) data. If you have a local copy of the
    PWT data that you wish to work with, then you can directly specify the path
    arugment to load your local copy directly.  Alternatively, you can download 
    a copy of the PWT data by specify a version number and date corresponding
    to the specific (default settings will grab the most recent update of the 
    most recent version of PWT). The required date format is 'MMDDYYYY'. 

    By default, any downloaded data is extracted into your current working 
    directory. If you don't want a local copy of the PWT data, then set extract 
    to false.
     
    Returns a pandas Panel object where the major axis is country (technically 
    country iso codes) and the minor axis is year.  Setting panel to false will 
    return a DataFrame object with the pandas default index.

    """
    # first check for a local copy of PWT
    if path != None and panel == False:
        pwt = pd.read_csv(path)

    elif path != None and panel == True:
        pwt = pd.read_csv(path, index_col=['year', 'isocode'])
        pwt = pwt.to_panel()

    # otherwise, download the appropriate zip file 
    elif path == None:
        url = 'http://pwt.econ.upenn.edu/Downloads/pwt' + str(version) + \
              '/pwt' + str(version) +'_' + date + 'version.zip' 
        archive = zipfile.ZipFile(StringIO(urlopen(url).read()), 'r') 

        # to extract or not to extract...
        tmp_file = 'pwt' + str(version) + '_wo_country_names_wo_g_vars.csv'
        if extract == True:
            archive.extractall()
            pwt = pd.read_csv(tmp_file)
        else:
            pwt = archive.read(tmp_file)
            pwt = pd.read_csv(StringIO(pwt)) 

        # Do you want a Panel or Dataframe? 
        if panel == True:
            # inelegant! Should use MultiIndex instead...
            pwt = pd.read_csv(tmp_file, index_col=['year', 'isocode'])
            pwt = pwt.to_panel()
        else:
            pass         
            
    return pwt
     

# TODO: Work out how to generate the panel groupings

#def intervals(x, n=5): 
#    """Maps years to n year groups, eg., 1960-n -> 1960""" 
#    return 1960 + 5 * int(np.floor((x - 1960.)/5)) 
#
#def panel_average(x): 
#    panel_grp = x.groupby(lambda x : x[1])
#    return panel_grp.aggregate('mean') # ignores NaNs 
#
#yr_group = pwt.groupby(intervals, level=0) 
#pwt_new = yr_group['rgdpl'].apply(panel_average)
#pwt_new = pwt_new.T.stack(dropna=False) 
#pwt_new.name = 'rgdpl'
