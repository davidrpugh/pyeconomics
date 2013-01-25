import zipfile 
from urllib import urlopen 
from StringIO import StringIO 

import numpy as np
import pandas as pd

def get_SolowResiduals(rgdppc, rgdppw, g0, delta, alpha, h=10, **kwargs):
    """Computes Solow residuals!

    Required arguments:

        1. rgdppc: Must specify a valid measure of real gdp per capita. 
           Valid options are rgdpl, rgdpl2, rgdpch. See Penn World 
           Tables documentation for definitions of these variables.
        2. rgdppw: Must specify a valid measure of real gdp per worker. 
           Valid options are rgdpwok, rgdpl2wok, rgdpl2pe, rgdpl2te. 
           See Penn World Tables documentation for definitions of 
           these variables.
        3. g0: Initial guess for the growth rate of technology. 
           Required in order to pin down an initial estimate of the 
           capital stock.
        4. delta: Estimated rate of capital decpreciation rate. 
           Required in order to compute the capital stock. 
        5. alpha: Estimated share of income/output going to capital. 
           Required in order to decompose the growth rates in order to 
           impute A and g.
        6. h (default=10): Must specify the amount of smoothing to be 
           applied to computed growth rates (i.e., those of labor 
           force, investment share, and technology). 

    """
        
    # optional keywords args
    path    = kwargs.get('path', None)
    version = kwargs.get('version', 71)
    date    = kwargs.get('date', '11302012')
    extract = kwargs.get('extract', True)
        
    # first check for a local copy of PWT
    if path != None:
        pwt = pd.read_csv(path, index_col=['year', 'isocode'])

    # otherwise, download the appropriate zip file 
    elif path == None:
        url = 'http://pwt.econ.upenn.edu/Downloads/pwt' + \
              str(version) + '/pwt' + str(version) +'_' + date + \
              'version.zip' 
        archive = zipfile.ZipFile(StringIO(urlopen(url).read()), 'r') 

        # to extract or not to extract...
        tmp_file = 'pwt' + str(version) + \
                   '_wo_country_names_wo_g_vars.csv'
        if extract == True:
            archive.extractall()
            pwt = pd.read_csv(tmp_file, index_col=['year', 'isocode'])
        else:
            pwt = archive.read(tmp_file)
            pwt = pd.read_csv(StringIO(pwt), index_col=['year', 'isocode'])         
    
    # convert to Pandas Panel object
    pwt = pwt.to_panel()
    
    # to compute Solow residuals need to a measure of real GDP...
    pwt['realGDP'] = pwt[rgdppc] * pwt.POP
    
    # ...and a measure of labor force growth rates 
    pwt['laborForce'] = (pwt[rgdppc] / pwt[rgdppw]) * pwt.POP
    pwt['laborForceGrowth'] = pwt.laborForce.pct_change() 

    # annual growth rates are noisy, smooth them! Note backshift!
    tmp_data = pd.rolling_mean(pwt.laborForceGrowth, 
                               window=10, min_periods=h).shift(-h)
    pwt['smoothedLaborForceGrowth'] = tmp_data
    
    # convert ki to proportion
    pwt['investmentShare'] = pwt.ki / 100.
    
    # investment shares are also noisy, smooth them! Note backshift!
    tmp_data = pd.rolling_mean(pwt.investmentShare, window=10, 
                               min_periods=h).shift(-h)
    pwt['smoothedInvestmentShare'] = tmp_data
    
    # initial data on imputed K is just made up of K0 (and junk!)
    break_even_s = (1 + pwt.smoothedLaborForceGrowth) * (1 + g0) - (1 - delta)
    actual_s = pwt.smoothedInvestmentShare 
    pwt['imputedK'] = pwt.realGDP * (actual_s / break_even_s)
    
    # impute capital stock for each country
    for ctry, value in pwt.imputedK.iteritems():
        for year in pwt.imputedK.index[:-1]: # no data beyond 2010
            if np.isnan(value.ix[year]) == False:
                k = pwt.imputedK[ctry].ix[year]
                s = pwt.investmentShare[ctry].ix[year]
                y = pwt.realGDP[ctry].ix[year]
                k_plus = (1 - delta) * k + s * y
                pwt.imputedK[ctry].ix[year + 1] = k_plus
            else:
                pass
            
    # create capital-output ratio
    pwt['capitalOutputRatio'] = pwt.imputedK / pwt.realGDP
    
    # compute the implied level of technology
    pwt['technology'] = pwt[rgdppw] / (pwt.capitalOutputRatio**(alpha / (1 - alpha)))
    
    # finally, compute the growth in technology
    pwt['technologyGrowth'] = pwt['technology'].pct_change()

    return pwt
