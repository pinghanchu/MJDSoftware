#!/usr/bin/python
import os
import re
import matplotlib.pyplot as plt

numbers = re.compile(r'(\d+)')
def numericalSort(value):
    parts = numbers.split(value)
    parts[1::2] = map(int, parts[1::2])
    return parts

prefixed = [filename for filename in os.listdir("/Users/daq/Desktop/"
        "PANDA-HPGe-Diode/DeadLayerTests/"
        "Script/08-11to14-2017_wFill2/output_files") if filename.startswith("output")]

sortedpf = []
for i in sorted(prefixed, key=numericalSort):
    sortedpf.append(i)

n = 0
ratios = []
disp = []
errors = []
for i in sortedpf:
    fin = open(i)
    for line in fin:
        line = line.rstrip()
        array = line.split(" ")
        if (float(array[4]) < 0):
            array[4] = 0.000
        if float(array[5]) < 0:
            array[5] = float(array[5])*(-1)
        if (float(array[4]) < 1) and (float(array[5]) < 1):
            ratios.append(float(array[4]))
            errors.append(float(array[5]))
            disp.append(n)
    fin.close()
    n = n + 1

print len(ratios),' ',len(disp),' ',len(errors)

### Plot ###

#plt.plot(ratios)
plt.errorbar(disp, ratios, xerr=0.1, yerr=errors)
plt.ylabel('Ratio of 81keV/356keV Counts')
plt.xlabel('Displacement from Center [mm]')
plt.title('Ratio Plot')
plt.savefig('RatioPlot.pdf')
plt.show()
