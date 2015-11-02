#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import argparse
import pandas
import sys
import math

TEXT_RED="\033[31m%s"
TEXT_GREEN="\033[32m%s"

def Square(x):
    return x*x

def Compare(path1, path2):
    df1 = pandas.read_table(path1, names=("t", "CD", "_", "CL", "Pin", "Pout"))
    df2 = pandas.read_table(path2, names=("t", "CD", "_", "CL", "Pin", "Pout"))

    error = math.sqrt(sum(map(Square, df1["CD"] - df2["CD"])))
    return Judge(error)

def Judge(error):
    THRESHOLD = 1.0e-6
    print("ERROR: %e" % error)
    print("THRESHOLD: %e" % THRESHOLD)
    if error < THRESHOLD:
        print(TEXT_GREEN % "TEST PASSED!")
        return 0
    else:
        print(TEXT_RED % "TEST FAILED!")
        return 1

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="This script checks the error of two aerodynamic loads")
    parser.add_argument("data1",
                        type=str,
                        help="path to data")
    parser.add_argument("data2",
                        type=str,
                        help="path to data")
    args = parser.parse_args()
    sys.exit(Compare(args.data1, args.data2))
