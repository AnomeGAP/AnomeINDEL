#!/bin/env python
#
# @note Copyright (C) 2015, Anome Incorporated. All Rights Reserved.
#       This program is an unpublished copyrighted work which is proprietary to
#       Anome Incorporated and contains confidential information that is not to
#       be reproduced or disclosed to any other person or entity without prior
#       written consent from Anome, Inc. in each and every instance.
#
# @warning Unauthorized reproduction of this program as well as unauthorized
#          preparation of derivative works based upon the program or distribution of
#          copies by sale, rental, lease or lending are violations of federal copyright
#          laws and state trade secret laws, punishable by civil and criminal penalties.
#
# @file    ExtractCandidateRefGene.py
#
# @brief   Extracting candidate from gene list
#
# @author  Chung-Tsai Su(chungtsai_su@anome.com)
#
# @date    2015/04/15
#
# @version 1.0
#
# @remark
#

import sys
import getopt
import os.path
import math
import re
from collections import defaultdict

def Usage():
    print("ExtractCandidateRefGene.py -i <Input RefGene file> -f <Candidate gene list> -o <Output RefGene file>") 
    print("Argument:")
    print("\t-h: Usage")
    print("\t-i: Input RefGene file")
    print("\t-f: Candidate gene list")
    print("\t-o: Output RefGene file")
    print("Usage:")
    print("\t./ExtractCandidateRefGene.py -i ../example/refGene.txt -f ../example/hearingloss_gene.list -o ../example/refGene.txt.hearingloss")

    return

def ExtractCandidateRefGene(ifn, ffn, ofn):
    rgfile = open(ifn, "r")
    ffile = open(ffn, "r")
    output = open(ofn, "w")
    num_reads = 0
    num_unmapped = 0
    hGeneList = {}

    for line in ffile:
        if re.match("^#", line):
           print "YES"
           continue

        line = line.strip()
        hGeneList[line] = 0
        #print "[%s]" % (line)

    for line in rgfile:
        line = line.strip()
        items = line.split("\t")
        name = items[12]
        if name in hGeneList:
            #print line
            output.write("%s\n" % (line))
            hGeneList[name] += 1

    for k in hGeneList.keys():
        if (hGeneList[k] == 0):
            print "%s\t%d" % (k, hGeneList[k])

    output.close()
    ffile.close()
    rgfile.close()
    return

def main(argv):
    ifile = ""
    ffile = ""
    ofile = ""

    try:
        opts, args = getopt.getopt(argv,"hi:f:o:")
    except getopt.GetoptError:
        Usage()
        sys.exit(1)

    for opt, arg in opts:
        if opt == '-h':
            Usage()
            sys.exit()
        elif opt in ("-i"):
            ifile = arg
        elif opt in ("-f"):
            ffile = arg
        elif opt in ("-o"):
            ofile = arg

    #error handling for input parameters
    if (ifile == ""):
        print("Error: '-i' is required")
        Usage()
        sys.exit(2)
    elif (os.path.isfile(ifile) == False): 
        print("Error: input file(%s) is not existed" % (ifile))
        Usage()
        sys.exit(3)

    if (ffile == ""):
        print("Error: '-f' is required")
        Usage()
        sys.exit(2)
    elif (os.path.isfile(ffile) == False): 
        print("Error: input file(%s) is not existed" % (ffile))
        Usage()
        sys.exit(3)

    if (ofile == ""):
        ofile = "%s.filtered" % (ifile)

    #Main Function
    ExtractCandidateRefGene(ifile, ffile, ofile)

    return

if __name__ == '__main__':
    main(sys.argv[1:])

