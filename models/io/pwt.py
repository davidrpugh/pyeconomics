from __future__ import division
import zipfile 
from urllib import urlopen 
from StringIO import StringIO 

import pandas as pd

def get_pwt_data_old(version=71, date='11302012', extract=False):
    """
    Load older (i.e., < 8.0) versions of the Penn World Tables data. If no local
    copy of the data is found, then a copy will be downloaded from the web.

    Arguments:
 
        version: (int) Version number for PWT data. Default is 71 (which is the 
                 most recent version).
        date:    (str) Date of the version update of the PWT data. Format is 
                 'mmddyyyy'. Default is '11302012' (which is the most recent 
                 update of version 71 of PWT).
        extract: (boolean) Whether or not you wish to save local copy of the 
                 PWT data in your working directory. Default is False.

    Returns:

        pwt:    A Pandas Panel object containing the Penn World Tables 
                Data along with the Solow residuals.

    TODO: Convert time index to datetime object.

    """        
    # first check for a local copy of PWT
    try: 
        path = 'pwt' + str(version) + '_wo_country_names_wo_g_vars.csv'
        pwt = pd.read_csv(path, index_col=['year', 'isocode'])

    # otherwise, download the appropriate zip file 
    except IOError:  
        url = ('http://pwt.econ.upenn.edu/Downloads/pwt' + str(version) + 
               '/pwt' + str(version) +'_' + date + 'version.zip') 
        archive = zipfile.ZipFile(StringIO(urlopen(url).read()), 'r') 

        # to extract or not to extract...
        tmp_file = 'pwt' + str(version) + '_wo_country_names_wo_g_vars.csv'
        
        if extract == True:
            archive.extractall()
            pwt = pd.read_csv(tmp_file, index_col=['year', 'isocode'])
        else:
            pwt = archive.read(tmp_file)
            pwt = pd.read_csv(StringIO(pwt), index_col=['year', 'isocode'])         
    
    # convert to Pandas Panel object
    pwt = pwt.to_panel()
    
    return pwt

def get_pwt_data(version=80, extract=True):
    """
    Load the Penn World Tables data. If no local copy of the data is found, then
    a copy will be downloaded from the web.

    Arguments:
 
        version: (int) Version number for PWT data. Default is 80 (which is the 
                 most recent version).
    
        extract: (boolean) Whether or not you wish to save local copy of the 
                 PWT data in your working directory. Default is True.

    Returns:

        pwt:    A Pandas Panel object containing the Penn World Tables data.

    TODO: Convert time index to datetime object.

    """        
    # first check for a local copy of PWT
    try: 
        path = 'pwt' + str(version) + '.dta'
        pwt = pd.read_stata(path)

    # otherwise, download the appropriate zip file 
    except IOError:  
        url = ('http://www.rug.nl/research/GGDC/data/pwt/V' + str(version) +
               '/pwt' + str(version) + '.zip') 
        archive = zipfile.ZipFile(StringIO(urlopen(url).read()), 'r') 

        # to extract or not to extract...
        tmp_file = 'pwt' + str(version) + '.dta'
        
        if extract == True:
            archive.extractall()
            pwt = pd.read_stata(tmp_file)
        else:
            pwt = archive.read(tmp_file)
            pwt = pd.read_stata(StringIO(pwt))         
    
    # convert to Pandas Panel object
    pwt.index = [pwt['year'], pwt['countrycode']]
    pwt = pwt.to_panel()
    
    return pwt
