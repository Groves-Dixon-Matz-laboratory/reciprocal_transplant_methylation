#!/usr/bin/env python
##assemble_htseq_counts.py
##written 2/25/16 by Groves Dixon
ProgramName = 'assemble_htseq_counts.py'
LastUpdated = '2/25/16'
By = 'Groves Dixon'
VersionNumber = '1.0'
print "\nRunning Program {}...".format(ProgramName)

VersionString = '{} version {} Last Updated {} by {}'.format(ProgramName, VersionNumber, LastUpdated, By)

Description = '''
Description:
Takes a set of gff counts files output from htseq and outputs
them as a combined table.

'''


##Import Modules 

import time
import argparse
from sys import argv
from sys import exit
import numpy as np
Start_time = time.time() ##keeps track of how long the script takes to run


##Set Up Argument Parsing
parser = argparse.ArgumentParser(description=Description) ##create argument parser that will automatically return help texts from global variables above
parser.add_argument('-i', required = True, dest = 'input', nargs = "+", help = 'Glob to the input files (eg *counts.gff)')
parser.add_argument('-o', required = True, dest = 'out', help = 'The desired name for the output file')
parser.add_argument('-pos', required = False, default = None, dest = 'position', help = 'Optional argument for position of informative part of file names. For exampel in this file sample1a__S_counts.gff, if all you need is ample1a, then put 1')
parser.add_argument('-delim', required = False, default = '_', dest = 'delim', help = 'The delmiter to use when picking a position for the informative part of the file name. Default = "_"')
args = parser.parse_args()

#Assign Arguments
infileList = args.input
infileList.sort()
outfileName = args.out
position = args.position
delimit = args.delim

def read_file(infileList):
    '''read in the files and build a dictionary
    '''
    print("\nReading in the following {} files:".format(len(infileList)))
    for file in infileList:
        print(file)
    headerList = ['geneID']
    dataDict = {}
    for file in infileList:
        if position:
            index = int(position) - 1
            headerList.append(file.split(delimit)[index])
        else:
            headerList.append(file)
        with open(file, 'r') as infile:
            for line in infile:
                line = line.strip("\n").split("\t")
                gene = line[0]
                count = line[1]
                try:
                    dataDict[gene].append(count)
                except KeyError:
                    dataDict[gene] = [count]
    return headerList, dataDict

def output(headerList, dataDict):
    '''Outputs all counts data as a combined table'''
    print "\nOutputting counts to file {}...".format(outfileName)
    with open(outfileName, 'w') as out:
        geneList = dataDict.keys()
        geneList.sort()
        out.write("\t".join(headerList))
        for gene in geneList:
            out.write("\n" + gene + "\t" + "\t".join(dataDict[gene]))
                
        


headerList, dataDict = read_file(infileList)
output(headerList, dataDict)

#return time to run
Time = time.time() - Start_time
print('\nTime took to run: {}'.format(Time))


