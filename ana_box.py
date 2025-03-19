#!/usr/bin/env python

import os
import csv
import numpy as np
import pandas as pd
import argparse as arg
import MDAnalysis as mda


def options():
    '''Defines the options of the script.'''

    parser = arg.ArgumentParser(
                formatter_class=arg.ArgumentDefaultsHelpFormatter)

    #
    # Input Options
    #
    inp = parser.add_argument_group("Input Data")

    inp.add_argument('-c', '--crd', type=str, dest='CrdFile',
                     required=True, help='''Coordinates File.''')

    inp.add_argument('-t', '--trj', type=str, dest='TrjFile',
                     default=None, help='''Trajectory File.''')

    #
    # Output Options
    #
    out = parser.add_argument_group("Output Options")

    out.add_argument('-o', '--output', default=None, type=str, dest='OutFile',
                     help='''Output File Prefix.''')

    args = parser.parse_args()
    Opts = vars(args)

    return Opts


if __name__ == '__main__':

    Opts = options()

    u = mda.Universe(Opts["CrdFile"], Opts["TrjFile"])
    
    if Opts["OutFile"]:
        outfile = Opts["OutFile"] + "_box.csv"
    else:
        outfile = "box.csv"

    boxes = []
    times = []
    for ts in u.trajectory:
        times.append(ts.time)
        box = u.dimensions.copy()
        boxes.append(box)

    times = np.array(times)
    boxes = np.array(boxes)
    data = np.c_[ times, boxes ]
    df = pd.DataFrame(columns=["t", "a", "b", "c", "alpha", "beta", "gamma"], data=data)
    df.to_csv(f'{outfile}', index=False, quoting=csv.QUOTE_NONNUMERIC)    
