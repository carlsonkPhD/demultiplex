#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  4 14:14:00 2022

@author: carlsokr
"""
import argparse
import os

parser = argparse.ArgumentParser(description='demultiplex forward reads from fastq')
parser.add_argument('--IDfile', '-i', dest='IDfile', type=str, nargs=1, help='name of file containing sample and barcode names in form "SAMPLE:BARCODE".', required = True)
parser.add_argument('--SearchString', '-s', dest='SearchString', type=str, nargs=1, default= 'AGATCGGAAGAGCACACGTCTGAA', help='sequence at 3\' end of barcode to be used to identify barcode in each read')
parser.add_argument('--Forwardfile', '-f', dest='Forwardfile', type=str, nargs=1, help='name of fastq file with forward reads', required = True)
parser.add_argument('--Reversefile', '-r', dest='Revfile',type=str, nargs=1, help='name of fastq file storing reverse reads', required = True)
args = parser.parse_args()

def demultiplex(Forwardfile, Revfile, IDfile, SearchString):
    #parse IDfile to obtain sample names and barcodes
    print("serching of sequence '"+SearchString+"' in forward reads")
    goodCodes = []
    R1_handles_string = {}
    R1_handles = {}
    R2_handles_string = {}
    R2_handles = {}

    with open(IDfile, "r") as ids:
        for line in ids:
            sample = line.rstrip().split(":", 1)[0]
            barcode = line.rstrip().split(":", 1)[1]
            if barcode in goodCodes:
                sys.exit("Duplicate barcodes detected in IDfile\nProper demultiplexing requires a unique barcode for each sample\n")
            else:
                goodCodes.append(barcode)
                R1_handles_string[barcode] = sample
                R2_name = sample + "R"
                R2_handles_string[barcode] = R2_name


    #open connection for each output file
    cwd = os.getcwd()
    numFiles = len(R1_handles_string)
    R1_count = 0
    R2_count = 0
    
    for key in R1_handles_string:
        my_handle = R1_handles_string[key]
        R1_handles[key] = open(my_handle+".fastq", "w")
    for key in R2_handles_string:
        my_handle = R2_handles_string[key]
        R2_handles[key] = open(my_handle+".fastq", "w")
        
    #print(R1_handles)
    #print(R2_handles)
    
    #initialize lists and counter
    count = 0
    codes = []
    myLines = []
    myR = []
    
    #open input fastq files and parse reads
    with open(Forwardfile, "r") as fp, open(Revfile, "r") as pf:
        for x, y in zip(fp,pf):
        
            #add line to myLines list
            myLines.append(x)
            myR.append(y)
            count += 1
            if count % 4 == 0:
                #search for sequence 3' of barcode in adapter "AGATCGGAAGAGCACACGTCTGAA" was previously used
                dex = myLines[1].find(SearchString)
                #if substring found, find() will return index, otherwise returns -1
                if dex != -1:
                    barcode = myLines[1][(dex-5):dex] #slice barcode
                    if barcode in goodCodes: #check that barcode is valid
                        #add barcode to codes list then write to file based on identity of barcode
                        codes.append(barcode)
                        my_handle = R1_handles[barcode]
                        my_handle.write(myLines[0])
                        my_handle.write(myLines[1])
                        my_handle.write(myLines[2])
                        my_handle.write(myLines[3])
                        #write reverse reads to file
                        my_Rhandle = R2_handles[barcode]
                        my_Rhandle.write(myR[0])
                        my_Rhandle.write(myR[1])
                        my_Rhandle.write(myR[2])
                        my_Rhandle.write(myR[3])
                #clear myLines object to free up memory        
                myLines = [] 
                myR = []
    
    #Close file connections
    for item in R1_handles:
        R1_handles[item].close()
    for item in R2_handles:
        R2_handles[item].close()

    #print counts of each barcode
    for key in goodCodes:
        print(key+" : "+str(sum(1 for i in codes if i == key)))
    print("total barcodes identified: " + str(len(codes)))

demultiplex(args.Forwardfile[0], args.Revfile[0], args.IDfile[0], args.SearchString)

