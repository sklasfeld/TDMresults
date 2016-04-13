import argparse
import re
import numpy as np


parser = argparse.ArgumentParser(description="Build Accuracy Table")
parser.add_argument('input_list', help='tab delimited text file containing name of \
                    accuracy files and prefix for the type of standardization \
                    that makes that file')
#parser.add_argument('out', help='The name of the output file')
args = parser.parse_args()

#class MagicString(str):
#	magicSplit = str.split
prcnt = 0.0
inp = open(args.input_list, 'r')
for line2 in inp:
    line2 = line2.strip()
    params = line2.split("\t")
    fa = open(params[0], 'r')
    #print (params)
    for line1 in fa:
        line1 = line1.strip()
        if(re.search('^Accuracy\s+', line1)):
            #print(line1)
            n_line = next(fa)
            data = n_line.split()
            #print (data)
            print ("%s\t%s\t%s\t%s\t%s" % (params[1], prcnt, data[0], data[2], data[3]))
            prcnt = np.round(prcnt + 1e-1,decimals=1)
            if prcnt > 1:
                prcnt = 0.0
