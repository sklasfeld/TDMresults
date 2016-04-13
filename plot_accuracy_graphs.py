import argparse
import numpy as np
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser(description="Plot Accuracy Graph")
parser.add_argument('table', help='Accuracy table must be tab-delimited format \
                                   with the following fields in this \
                                   particular order: Quanitification method, \
                                   Percent of RNA seq in training Data, \
                                   Accuracy, Accuracy Lower, Accuracy Upper')
parser.add_argument('title', help='title of graph')
parser.add_argument('out', help='The name of the output image')
args = parser.parse_args()

#data_arrays
logx = []
logy = []
logmin = []
logmax = []
npnx = []
npny = []
npnmin = []
npnmax = []
qnx = []
qny = []
qnmin = []
qnmax = []
tdmx = []
tdmy = []
tdmmin = []
tdmmax = []
unx = []
uny = []
unmin = []
unmax = []

#fill out data arrays
tbl = open(args.table, 'r')
for lines in tbl:
    line = lines.split("\t")
    if line[0] == "log":
        logx.append(float(line[1]))
        logy.append(float(line[2]))
        logmin.append(float(line[3]))
        logmax.append(float(line[4]))
    if line[0] == "npn":
        npnx.append(float(line[1]))
        npny.append(float(line[2]))
        npnmin.append(float(line[3]))
        npnmax.append(float(line[4]))
    if line[0] == "qn":
        qnx.append(float(line[1]))
        qny.append(float(line[2]))
        qnmin.append(float(line[3]))
        qnmax.append(float(line[4]))
    if line[0] == "tdm":
        tdmx.append(float(line[1]))
        tdmy.append(float(line[2]))
        tdmmin.append(float(line[3]))
        tdmmax.append(float(line[4]))
    if line[0] == "un":
        unx.append(float(line[1]))
        uny.append(float(line[2]))
        unmin.append(float(line[3]))
        unmax.append(float(line[4]))

#convert arrays to numpy arrays
logx = np.array(logx)
logy = np.array(logy)
logmin = np.array(logy) - np.array(logmin)
logmax = np.array(logmax) - np.array(logy)
npnx = np.array(npnx)
npny = np.array(npny)
npnmin = np.array(npny) - np.array(npnmin)
npnmax = np.array(npnmax) - np.array(npny)
qnx = np.array(qnx)
qny = np.array(qny)
qnmin = np.array(qny) - np.array(qnmin)
qnmax = np.array(qnmax) - np.array(qny)
tdmx = np.array(tdmx)
tdmy = np.array(tdmy)
tdmmin = np.array(tdmy) - np.array(tdmmin)
tdmmax = np.array(tdmmax) - np.array(tdmy)
unx = np.array(unx)
uny = np.array(uny)
unmin = np.array(uny) - np.array(unmin)
unmax = np.array(unmax) - np.array(uny)

# plot Graph
plt.figure()
line1 = plt.errorbar(logx, logy, yerr=[logmin, logmax], fmt='--o', color='r')
line2 = plt.errorbar(npnx, npny, yerr=[npnmin, npnmax], fmt='--o', color='b')
line3 = plt.errorbar(qnx, qny, yerr=[qnmin, qnmax], fmt='--o', color='g')
line4 = plt.errorbar(tdmx, tdmy, yerr=[tdmmin, tdmmax], fmt='--o', color='m')
line5 = plt.errorbar(unx, uny, yerr=[unmin, unmax], fmt='--o', color='c')
first_legend = plt.legend([line1,line2,line3,line4,line5], ['log','npn','qn','tdm','un'], loc=4)


plt.suptitle(args.title)
plt.savefig(args.out, bbox_inches='tight')
