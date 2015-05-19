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
# @file    SegmentedSuffix.py
#
# @brief   Given a DNA sequence (FASTA file), segmenting them by a user defined length and 
#          then generate all of suffix for each segment.
#
# @author  Chung-Tsai Su(chungtsai_su@anome.com.tw)
#
# @date    2015/05/14
#
# @version 0.1
#
# @remark
#           2015/5/14 (0.1) generate rowkey and its column


import sys
import getopt
import os.path
import math
import re
from collections import defaultdict
from datetime import datetime
from sets import Set

SEGMENT_LENGTH=200
MINIMAL_LENGTH=10
gStartTtime = datetime.now()

def TimeLog_init():
    global gStartTtime 
    gStartTtime = datetime.now()
    return

def TimeLog_write(msg):
    global gStartTtime
    lEndTime = datetime.now()
    print('{}: %s'.format(lEndTime - gStartTtime) % (msg))
    gStartTtime = lEndTime
    return

def Usage():
    print("SegmentedSuffix.py -i <Input FASTA file> -l <segment length> -d <minimal length of suffix> -o <Output CSV file>") 
    print("Argument:")
    print("\t-h: Usage")
    print("\t-i: Input FASTA file")
    print("\t-l: Segmentation length [Default:%d]" % (SEGMENT_LENGTH))
    print("\t-d: Minimal length of suffix [Default:%d]" % (MINIMAL_LENGTH))
    print("\t-o: Output CSV file")
    print("Usage:")
    print("\t./SegmentedSuffix.py -i ../example/chr1.fa -l 200 -d 10 -o ../example/chr1.suffix.csv")

    return

def SuffixExtend(items, segment, min_len, idx):
    length = len(segment)-min_len+1
    for i in range(0,length):
        # remove the prefix 'N'
        if segment[i] == 'N': continue

        # remove the suffix 'N'
        e = segment.find('N', i+1)
        if e == -1: e = len(segment)
        #print "e=%d [%s]" % (e, segment[i:])

        if i+min_len > e: continue

        item = "%s %d" % (segment[i:e], (idx+i))
        items.append(item)
        #print item

    return

def AddCommonPrefix(Dprefix, DSamePrefix, R2P, P2S, value):
    item = "%s %s" % (R2P, P2S)
    #print "[%s] %s" % (item, value)
    Dprefix[item] = Dprefix[item].union(eval(value))
    
    ## padding position information if same prefix
    key = ".%s" % (R2P)
    DSamePrefix[key].add(item)
    #print "%s add [%s]" % (key, item)
    
    return

def UpdatePosition(Dprefix, DSamePrefix, seq, value):
    if (seq == ""): return

    ##debug
    #return

    length = len(seq) - 1
    while length >= 0:
        item = ".%s" % (seq[:length])
        #print "check [%s]" % (item)
        if item in DSamePrefix:
            for key in DSamePrefix[item]:
                #print "Add %s into %s [%s][%s]" % (value, key, seq, item)
                Dprefix[key] = Dprefix[key].union(value)

        length = length - 1
    
    return

def CommonPrefixExtract(items):
    Dprefix = defaultdict(Set)
    DSamePrefix = defaultdict(Set)

    #items = ['AAAAAC 0', 'AAAAC 1', 'AAAC 2', 'AAAC 10', 'AAC 3', 'AC 4', 'AC 20', 'C 5', 'C 60', 'GAAA 100', 'GACC 200', 'TAAA 300', 'TTTT 400']
    #items = ['AAAAA 1', 'AAAAA 10', 'AACAA 15', 'AACAA 33', 'AAGGG 25', 'CAAAA 9', 'CAACA 27', 'CAACA 52', 'CACCC 37', 'GGGGG 70']

    ##Round 1: aggregate the same suffix sequence together
    previous = ""
    value = []
    aggregated = []
    for item in items:
        (seq, idx) = item.split(' ', 1)
        if (previous == ''):
            previous = seq
            value.append(idx)
        elif (seq != previous):
            temp = "%s %s" % (previous, value)
            #print temp
            aggregated.append(temp)

            #reset
            previous = seq
            value = [idx]
        else:
            value.append(idx)
    if (len(value) > 0):
        temp = "%s %s" % (previous, value)
        #print temp
        aggregated.append(temp)
    TimeLog_write("     Round 1: aggregate the same suffix")

    ##Round 2: Extract common prefix 
    previous = ""
    pre_idx = []
    for item in aggregated:
        (seq, idx) = item.split(' ', 1)
        #print "%s %s" % (seq, idx)
        if (previous != ''):
            length = min(len(seq), len(previous))
            i = 0
            while (i < length):
                if (seq[i] != previous[i]): break
                i+=1
            AddCommonPrefix(Dprefix, DSamePrefix, previous[:i], previous[i:], pre_idx)
            AddCommonPrefix(Dprefix, DSamePrefix, seq[:i], seq[i:], idx)
        else:
            AddCommonPrefix(Dprefix, DSamePrefix, "", seq, idx)
        AddCommonPrefix(Dprefix, DSamePrefix, seq, "", idx) 
 
        #reset
        previous = seq
        pre_idx = idx
    TimeLog_write("     Round 2: Extract common prefix")

    ##Round 3: update position information
    items = Dprefix.keys()
    items.sort()
    i = 0
    for item in items:
        #if (i % 100 = 0):
        #    print "%d" % (i)
        i += 1
        (part1, part2) = item.split(' ', 1)
        #print "[%s] [%s] [%s]" % (part1, part2, Dprefix[item])
        ## Update position information if same prefix for each previous record
        if (part2 != ""):
            UpdatePosition(Dprefix, DSamePrefix, part1, Dprefix[item])
    TimeLog_write("     Round 3: update position information")


    ## Round 4: Generate prefix list
    ## Covert dict to list
    keys = Dprefix.keys()
    keys.sort()
    prefixlist = []
    for k in keys:
        prefixlist.append("%s %s" % (k, Dprefix[k]))

    print "Total of items for common prefix list : %d" % (len(Dprefix))    
    TimeLog_write("     Round 4: Generate prefix list")

    return prefixlist

def AddRowkey(rowkeylist, rowkey, R2P, P2S, Schildern, Spos):
    item = "[%s] [%s] [%s] [%s] %d [%s]" % (rowkey, R2P, P2S, Schildern, len(Spos), Spos)
    #print item
    rowkeylist.append(item)
    return

def RowkeyExtract(items):
    Lrowkey = []

    ##Round 3: Extract rowkey 
    pre_R2P = ""
    pre_part1 = ""
    Schildren = Set()
    R2P = ""
    P2S = ""
    Spos = Set()
    #print "Rowkey:"
    for item in items:
        ## common suffix, uncommon suffix, position list
        (part1, part2, Lpos) = item.split(' ', 2)
        Lpos = eval(Lpos)
        #print "%s %s$ %s" % (part1, part2, Lpos)

        if (part2 == ""): part2 = "$"

        if (pre_part1 == part1):
            Schildren.add( part2[:1] )
            Spos = Spos.union(Lpos)
        else:
            ##output
            rowkey = "%s.%s" % (R2P, P2S[:1])
            AddRowkey(Lrowkey, rowkey, R2P, P2S, Schildren, Spos)
            ##reset
            idx = part1.find(pre_part1)
            if (idx != 0):
                #restart from root
                length = min(len(part1), len(pre_part1))
                i = 0
                while (i < length):
                    if (part1[i] != pre_part1[i]): break
                    i+=1
                idx = i
            else:
                idx += len(pre_part1)

            #pre_R2P = R2P
            pre_part1 = part1
            R2P = part1[:idx]
            P2S = part1[idx:]
            Schildren = Set([part2[:1]])
            Spos = Set(Lpos)
            #print "idx=%d pre_R2P=[%s] pre_part1=[%s] R2P=[%s], P2S=[%s]" % (idx, pre_R2P, pre_part1, R2P, P2S)
    if (pre_part1 != ""):
        rowkey = "%s.%s" % (R2P, P2S[:1])
        AddRowkey(Lrowkey, rowkey, R2P, P2S, Schildren, Spos)

    return Lrowkey

def OutputFile(ofn, items):
    ofd = open(ofn, "w")
    for item in items:
        ofd.write("%s\n" % item)

    ofd.close()
    return    

def SegmentedSuffix(ifn, slength, min_len, ofn):
    total_length = 0
    idx = 0
    seq = ""
    items = []
    moving_window = int(slength / 2)

    ifd = open(ifn, "r")
    ##Read sequence and segment them with length of slength
    for line in ifd:
        if re.match("^>", line):
            #print "Skip comment: %s" % (line)
            continue

        line = line.strip()
        #print "line : %s" % (line)
        seq += line.upper()
        length  = len(seq)
        while length > slength:
            segment = seq[:slength]
            SuffixExtend(items, segment, min_len, idx) 

            #move to the next segment
            seq = seq[moving_window:]
            idx += moving_window
            total_length += moving_window
            length = len(seq)

    ##Extend the last segment of a given sequence
    if (length > 0):
        SuffixExtend(items, seq, min_len, idx)
        total_length += len(seq)
    print "Total length of input sequence = %d" % (total_length)
    print "Number of suffix = %d" % (len(items))

    TimeLog_write("SuffixExtend()")
    ifd.close()

    #Sort all of suffix
    items.sort()
    TimeLog_write("items.sort()")

    ##Extract all of common prefix by pairwise comparison
    prefixlist = CommonPrefixExtract(items)
    TimeLog_write("CommonPrefixExtract()")

    #Sort all of rowkey
    ##prefixlist.sort()
    TimeLog_write("prefixlist.sort()")

    ###Extract all of rowkey by sorted prefix list
    rowkeylist = RowkeyExtract(prefixlist)
    TimeLog_write("RowkeyExtract()")

    #Output rowkeys and their columns
    OutputFile(ofn, rowkeylist)
    TimeLog_write("OutputFile()")
    return

def main(argv):
    ifile = ""
    slength = SEGMENT_LENGTH
    min_len = MINIMAL_LENGTH
    ofile = ""

    TimeLog_init()

    try:
        opts, args = getopt.getopt(argv,"hi:l:d:o:")
    except getopt.GetoptError:
        Usage()
        sys.exit(1)

    for opt, arg in opts:
        if opt == '-h':
            Usage()
            sys.exit()
        elif opt in ("-i"):
            ifile = arg
        elif opt in ("-l"):
            slength = int(arg)
        elif opt in ("-d"):
            min_len = int(arg)
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

    if (ofile == ""):
        ofile = "%s.suffix.csv" % (ifile)


    #Main Function
    SegmentedSuffix(ifile, slength, min_len, ofile)

    return

if __name__ == '__main__':
    main(sys.argv[1:])

