# energy profile fixer v0.1 beta
# by Haohao Fu (fhh2626_at_gmail.com)
#
# this script fixes the result of a 2d energy profile (usually a CV-space pair interaction)
# it normalizes the data by letting the minimal value to zero
# and fixes the missing values
#
# it can also do smoothing
#
# usage:
#    python energyProfileFix xxxxxx.energy, xxxxx.energy 0
#    (the last number 0/1 determines whether one wants to do smoothing
#

import sys

import numpy as np
import pandas as pd
import scipy.signal as signal
from sklearn.impute import KNNImputer

def fixFile(files, filt=False):

    for oneFile in files:

        outputFile = oneFile + '.fixed'
        
        # read data
        df = pd.read_table(oneFile, sep = ' ', comment='#', names=['x', 'y', 'z'])
        df2 = df.pivot(index='y', columns='x', values='z')
        df2 = df2.replace('-nan(ind)', np.nan)
        df2 = df2.astype(float)
        
        # subtract the minimal value
        minvalue = float(df2.min().min())
        df2 = df2 - minvalue
        
        # guess the missing value
        values = df2.values
        imputer = KNNImputer(n_neighbors=2, weights="uniform")
        newValues = imputer.fit_transform(values)
        
        # filt
        if filt:
            newValues = signal.medfilt(newValues,(3,3))
            
            # when doing this, the values at margin will be set to zero
            # hence we need to make them as "missing values"
            # and do the steps above again
            
            df2.values[:] = newValues
            
            df2 = df2.replace(0, np.nan)
            minvalue = float(df2.min().min())
            df2 = df2 - minvalue
            values = df2.values
            
            imputer = KNNImputer(n_neighbors=2, weights="uniform")
            newValues = imputer.fit_transform(values)
        
        df2.values[:] = newValues
        
        # reformat
        result = df2.unstack().reset_index()
        np.savetxt(outputFile, result.values, fmt='%g')

if __name__ == "__main__":
    filt = bool(int(sys.argv[-1]))
    fixFile(sys.argv[1:-1], filt)