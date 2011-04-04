#!/usr/bin/env python
# convertTabularDataToEdgeList.py
# Converts tab-delimited text files to edge lists
#
# Call as follows:
#   convertTabularDataToEdgeList.py input.csv output.adj key.txt
# Where input.csv is composed of tab-delimited data, with headings
# specifying which column corresponds to the start (CUE) and end
# (TARGET) nodes of a directed edge with weight (VALUE).
# 
# example input.csv:
#   HeadA   HeadB   HeadC   HeadD
#   USA     UK      15      Miscelanneous extra information
#   UK      USA     5       Random stuff
#   TW      USA     7
# Then CUE="HeadA", TARGET="HeadB", and VALUE="HeadC"
#
# output.adj:
#   1   2   15
#   2   1   5
#   3   1   7
#
# key.txt
#   1   USA
#   2   UK
#   3   TW
#
# ***NO overwrite checks are in place!!! Make certain that output.adj
# and key.txt are not files you need!!!***
#
# Whitespace in the key values is *not* encouraged and will probably break the script.

# Parameters:

DELIMITER='\t' # Tabs are the delimiter used in the source document

#(need to match column headings):
CUE='CUE'
TARGET='TARGET'
VALUE='FSG'



#######################################################################

import csv
import sys

if len(sys.argv) < 4:
    print 'convertTabularDataToEdgeList.py needs at least 3 arguments to function:'
    print '\t convertTabularDataToEdgeList.py input.csv output.adj key.txt'
    sys.exit(2)

keyDict={}
def key2num(key):
    if key in keyDict:
        return keyDict[key]
    else:
        keyDict[key]=len(keyDict)+1
        return keyDict[key]


wordReader = csv.DictReader(open(sys.argv[1],'rb'), delimiter=DELIMITER, quotechar='\\')
adjWriter = csv.writer(open(sys.argv[2],'wb'), delimiter=' ', quotechar='\\')
keyWriter = csv.writer(open(sys.argv[3],'wb'), delimiter='\t', quotechar='\\')
for row in wordReader:
    #print row
    # You can select only particular rows to materialize as edges by using the other columns as controls if you want
    #if row['NORMED?']=='YES':
        adjWriter.writerow([key2num(row[CUE]),key2num(row[TARGET]),row[VALUE]])

# Sort the keyDict by value
from operator import itemgetter
sortedKeyDict = sorted(keyDict.items(), key=itemgetter(1))
for (k,v) in sortedKeyDict:
    keyWriter.writerow([v,k])

