#!/usr/bin/env python

import sys
import numpy as np
import pandas as pd


if __name__ == '__main__':
    data = np.loadtxt("Scan.dat")
    data = data[:,:2]
    data = data[data[:,0].argsort()]
    df = pd.DataFrame({"angle": data[:,0], "energy" : data[:,1]})
    lowest = df.groupby('angle', group_keys=False).apply(lambda x: x.loc[x.energy.idxmin()])
    np.savetxt("Scan_lowest.dat", lowest.values, fmt="%6.2f %12.6e")
