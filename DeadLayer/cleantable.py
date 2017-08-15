#!/usr/bin/python
#
user_input = raw_input("Please input data file:")
fin1 = open(user_input)
rawdata = []
data = []
count = 0
startline = 0
endline = 0
flag = 0
for line in fin1:
    line = line.rstrip()
    array = line.split(" ")
    rawdata.append(array[0])
    if array[0] == "$DATA:" :
        startline = count
    if array[0] == "$ROI:" :
        endline = count
    count = count + 1
fout = open("inputdata.txt","w")
for i in range(startline+1,endline-1):
    fout.write(rawdata[i]+"\n")
fout.close()
