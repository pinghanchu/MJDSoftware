#!/usr/bin/python
#
# This code reads the information from the NNDC database and
# extract beta decay information. The most useful information
# is the end-point energy of the beta decay and the intensity.
#
import re
import json
import urllib
import string
from BeautifulSoup import *

def betalist(nuc):
    serviceurl = 'http://www.nndc.bnl.gov/chart/decaysearchdirect.jsp?'
    url = serviceurl+urllib.urlencode({'nuc':nuc})
    #print 'Retrieving', url
    uh = urllib.urlopen(url)
    count1 = list()
    count2 = list()
    count = 0
    for line in uh:
        if re.search('^<p><u>Beta-</u>',line):
            count1.append(count)
        if re.search('<p>Mean beta- energy',line):
            count2.append(count)
        count = count+1
#    print count1, count2
    j = 0
    blist=list()
    while j< len(count1):
#        print j
        data = list()
        i = 0
        uh1 = urllib.urlopen(url)
        for line in uh1:
            if i in range(count1[j],count2[j]):
                data.append(line)
            i = i+1 
        count = 0
        y = list()
        for line1 in data:
            x = re.findall('.&nbsp;([0-9*]....[0-9*])|.&nbsp;([0-9*]...[0-9*])|.&nbsp;([0-9*]..[0-9*])|.&nbsp;([0-9*].[0-9*])|.&nbsp;([0-9*][0-9*])',line1)
            if len(x) > 0 :
                for k in range(0,len(x[0])):
                    if x[0][k] != '' :
                        z=float(x[0][k])
                        y.append(float(x[0][k]))
                        #print j, count, z  
            count = count +1
        j= j+1
        blist.append(y)
    return blist

def isotope(name):
    # Output "Energy", "End-Point Energy", "Intensity", "Dose"
    ll = betalist(name)
    with open('./data/'+name+'_betalist.txt','w') as target:
        for i in range(0,len(ll)):
            for j in range(0,len(ll[i])/4):
            #print i, j, ll[i][j*4+0],ll[i][j*4+1],ll[i][j*4+2],ll[i][j*4+3]
                target.write(str(i) +' ' + str(ll[i][j*4+0]) +' '+ str(ll[i][j*4+1]) + ' '+str(ll[i][j*4+2])+' '+str(ll[i][j*4+3])+'\n')
    target.close()


def isotopelist(list):
    name = []
    with open('./data/n'+list+'_list.txt') as f1:
        for line in f1:
            data = line.split()
            name.append(data[0])
        for i in range(0,len(name)):
            xx = name[i]
            disotop = isotope(xx)

isotopelistname = []
isotopelistname.append('Ge70')
isotopelistname.append('Ge72')
isotopelistname.append('Ge73')
isotopelistname.append('Ge74')
isotopelistname.append('Ge76')
for i in range(0,len(isotopelistname)):
    isotopelist(isotopelistname[i])
