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
# @file    FindLongDeletion.py
#
# @brief   Finding Long Deletion from BAM file
#
# @author  Chung-Tsai Su(chungtsai_su@anome.com)
#
# @date    2015/02/12
#
# @version 1.0
#
# @remark
#

import sys
import getopt
import os.path
import math
from collections import defaultdict
import pysam

#Default Parameter
MIN_SEGMENTLENGTH = 500
MAX_SEGMENTLENGTH = 50000
MIN_QUALITY = 10
MIN_COVERAGE = 0.9
MIN_SUPPORT = 5


def Usage():
    print("FindLongDeletion.py -i <Input BAM> -o <Output SAM> -m <Minimal Segment Length> -M <Maximal Segment Length> -q <Minimal Quality> -c <Minimal Coverage> -s <Minimal Number of Supports>") 
    print("Argument:")
    print("\t-h: Usage")
    print("\t-i: Input BAM file")
    print("\t-o: Output SAM file which contains all read pairs within specific segment length [Default=<inputfile>.ld.sam]")
    print("\t-m: Minimal segment length [Default=%d]" % (MIN_SEGMENTLENGTH) )
    print("\t-M: Maximal segment length [Default=%d]" % (MAX_SEGMENTLENGTH) )
    print("\t-q: Minimal read quality [Default=%d]" % (MIN_QUALITY) )
    print("\t-c: Minimal ratio of overlapping reads [Default=%f]" % (MIN_COVERAGE) )
    print("\t-s: Minimal number of supports [Default=%d]" % (MIN_SUPPORT) )
    print("Usage:")
    print("\t./FindLongDeletion.py -i ../../github/pydata-book/TSC320_S34_L001.srtd.reh.ddup.bam -o TSC320.sam -m 500 -M 50000 -q 10 -c 0.9 -s 5")

    return

def CalCoverage(r1, r2):
    ''' Description: Calcuate the percentage of overlapped region between two read pairs.
        Argument: rl - read-pair 1
                  r2 - read-pair 2
        Return: percentage of overlapped regions
    '''    
    percentage = 0.0
    p1 = min(r1.pos, r2.pos)
    p2 = max(r1.pos, r2.pos)
    p3 = min(r1.mpos, r2.mpos)
    p4 = max(r1.mpos, r2.mpos)
    if (p1 > p4):
        percentage = 0
    else:
        percentage = float(p3-p2) / float(p4 - p1)
    #print "(%d,%d,%d,%d)=>(%d,%d,%d,%d) %f" % (r1.pos,r1.mpos,r2.pos,r2.mpos,p1,p2,p3,p4,percentage)
    
    return percentage

def FindLongDeletion(cl, min_support, min_coverage):
    ''' Description: find the long deletion sorted by their support.
        Argument: cl - list of candidated reads
                  min_support - number of minimal supports
                  min_coverage - percentage of minimal coverage
        Return: list of long deletion with supported reads
    '''
    output = []
    num = len(cl)
    if num < 2:
        return output
    max_support = 0
    max_idx = 0
    items = []
    for i in range(0,num-1):
        item = {}
        item['start'] = cl[i].pos
        item['end'] = cl[i].mpos
        item['length'] = cl[i].mpos - cl[i].pos
        item['support'] = 1
        item['neighbor'] = []
        for j in range(i+1, num):
            if (cl[i].rname != cl[j].rname):
                continue
            if (cl[i].mpos < cl[j].pos):
                continue
            if (CalCoverage(cl[i], cl[j]) >= min_coverage):
                item['support'] += 1
                item['neighbor'].append(j)
                #print "(%s,%d,%d,%d,%d,%d,%s)" % (cl[i].qname, cl[i].rname, cl[i].pos, cl[i].mpos, cl[j].pos, cl[j].mpos, cl[j].qname)
        items.append(item)

    flag = [0] * num
    for i in range(0,num-1):
        if ((items[i]['support'] < min_support) or (flag[i] == 1)):
            continue
        
        item = {}
        item['support'] = items[i]['support']
        item['reads'] = [cl[i]]
        flag[i] = 1
        for j in items[i]['neighbor']:
            item['reads'].append(cl[j])
            flag[j] = 1
        output.append(item)
        
    return output

def LongDeletionIdentification(ifn, ofn):
    samfile = pysam.AlignmentFile(ifn, "rb")
    output = pysam.AlignmentFile(ofn, "wh", template=samfile)
    logfn = "%s_len_%d_%d.cov_%.2f.sup_%d.qual_%d.log" % (ofn, MIN_SEGMENTLENGTH, MAX_SEGMENTLENGTH, MIN_COVERAGE, MIN_SUPPORT, MIN_QUALITY)
    outputlog = open (logfn, "w")

    num_reads = 0
    num_longsegments = 0
    num_candidates = 0
    dicts = defaultdict(int) #values will initialize to 0
    candidate_list = []
    read_list = []

    for read in samfile.fetch():
        if ((read.rname != read.mrnm) or (read.is_unmapped == True)):
            continue
        
        gap = abs(read.pos-read.mpos)
        if (gap >= MIN_SEGMENTLENGTH):
            if (gap > MAX_SEGMENTLENGTH):
                gap = MAX_SEGMENTLENGTH 
            else:
                if ((read.pos <= read.mpos) and (read.mapq >= MIN_QUALITY)):
                    read_list.append(read)
                    num_candidates += 1
                num_longsegments += 1
                output.write(read)
    
        num_reads += 1
        dicts[gap] += 1

    candidate_list = FindLongDeletion(read_list, MIN_SUPPORT, MIN_COVERAGE)

    print "Number of reads = %d, Number of long segment reads=%d, Number of candidated reads=%d, Number of long deletions=%d" % ( num_reads, num_longsegments, num_candidates, len(candidate_list) )

    #Print the result
    print("SegmentLength=[%d,%d], Coverage=%.2f, Support=%d, Quality=%d\n" % (MIN_SEGMENTLENGTH, MAX_SEGMENTLENGTH, MIN_COVERAGE, MIN_SUPPORT, MIN_QUALITY))
    outputlog.write("SegmentLength=[%d,%d], Coverage=%.2f, Support=%d, Quality=%d\n" % (MIN_SEGMENTLENGTH, MAX_SEGMENTLENGTH, MIN_COVERAGE, MIN_SUPPORT, MIN_QUALITY))
    for i in range(len(candidate_list)):
        print( "#Support=%d\t\t\tRead Name\t\t\tChr\t\tStart\t\tEnd\t\tLength\t\tQuality\n" % (candidate_list[i]['support']))
        outputlog.write( "#Support=%d\t\t\tRead Name\t\t\tChr\t\tStart\t\tEnd\t\tLength\t\tQuality\n" % (candidate_list[i]['support']))
        for j in candidate_list[i]['reads']:
            print( "\t\t%s\t%d\t%16d\t%16d\t%16d\t%16d\n" % (j.qname,j.rname+1, j.pos, j.mpos, j.mpos-j.pos, j.mapq))
            outputlog.write( "\t\t%s\t%d\t%16d\t%16d\t%16d\t%16d\n" % (j.qname,j.rname+1, j.pos, j.mpos, j.mpos-j.pos, j.mapq))

    outputlog.close()
    output.close()
    samfile.close()
    return

def main(argv):
    global MIN_SEGMENTLENGTH
    global MAX_SEGMENTLENGTH
    global MIN_QUALITY
    global MIN_COVERAGE
    global MIN_SUPPORT
    ifile = ""
    ofile = ""

    try:
        opts, args = getopt.getopt(argv,"hi:o:m:M:q:c:s:")
    except getopt.GetoptError:
        Usage()
        sys.exit(1)

    for opt, arg in opts:
        if opt == '-h':
            Usage()
            sys.exit()
        elif opt in ("-i"):
            ifile = arg
        elif opt in ("-o"):
            ofile = arg
        elif opt in ("-m"):
            MIN_SEGMENTLENGTH = int(arg)
        elif opt in ("-M"):
            MAX_SEGMENTLENGTH = int(arg)
        elif opt in ("-q"):
            MIN_QUALITY = int(arg)
        elif opt in ("-c"):
            MIN_COVERAGE = float(arg)
        elif opt in ("-s"):
            MIN_SUPPORT = int(arg)

    #error handling for input parameters
    if (ifile == ""):
        print("Error: '-i' is required")
        Usage()
        sys.exit(2)
    elif (os.path.isfile(ifile) == False): 
        print("Error: input file(%s) is not existed" % (ifile))
        Usage()
        sys.exit(3)
    if (ofile == ""):
        ofile = "%s.ld.sam" % (ifile)

    #Main Function
    LongDeletionIdentification(ifile, ofile)

    return

if __name__ == '__main__':
    main(sys.argv[1:])

