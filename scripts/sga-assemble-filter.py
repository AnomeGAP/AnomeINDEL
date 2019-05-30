#!/bin/env python3
#
# @note Copyright (C) 2016, Anome Incorporated. All Rights Reserved.
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
# @file    VCFQualityFilter.py
#
# @brief   An example to filter Quality and Coverage on VCF file (Delly)
#
# @author  A-Tsai Su(chungtsai_su@anome.com.tw)
#
# @date    2016/01/23
#
# @version 1.0
#
# @remark
#

import sys
import getopt
import os.path
import collections

def Usage():
    ''' Description: Program Usage
        Argument:    NONE
        Return:	     NONE
    '''
    print("sga-assemble-filter.py -i <SGA contig file> -l <minimal length> -o <output file>")
    print("Argument:")
    print("\t-h: Usage")
    print("\t-i: <SGA contig file>")
    print("\t-l: <minimal contig length>")
    print("\t-o: <output fasta file>")
    print("Usage:")
    print("\t./sga-assemble-filter.py -i sga-contig.fa -l 1000 -o sga-contig-1000.fa")
    print("\tpython ./sga-assemble-filter.py -i all-bp2/all-bp-assemble-contigs.fa -l 142 -o "
          "all-bp2/all-bp-assemble-contigs.142.fa")

    return


def sga_assemble_filter(ifn, ofn, length):

    ifd = open(ifn, "r")
    ofd = open(ofn, "w")
    b_write = False
    h_cnt = collections.defaultdict(int)
    max_len = 0
    num_filtered = 0
    num_pass = 0
    total_bases = 0
    l_length = []

    for line in ifd:
        line = line.strip()
        if line.startswith(">"):
            items = line.split(' ')
            contig_len = int(items[1])
            if max_len < contig_len:
                max_len = contig_len
            h_cnt[int(contig_len/100)] += 1
            if contig_len >= length:
                b_write = True
                ofd.write("%s\n" % line)
                num_pass += 1
                total_bases += contig_len
                l_length.append(contig_len)
            else:
                b_write = False
                num_filtered += 1

        elif b_write:
            ofd.write("%s\n" % line)

    ifd.close()
    ofd.close()

    for k in h_cnt.keys():
        print("%d\t%d\t%d" % (k*100, (k+1)*100-1, h_cnt[k]))

    print("Maximal contig length = %d" % max_len)
    print("Total bases = %d" % total_bases)
    print("Number of passed contigs = %d" % num_pass)
    print("Mean of contig = %d" % (int(total_bases/num_pass)))
    print("There are %d contig removed (due to < %d)" % (num_filtered, length))

    l_length.sort(reverse=True)
    threshold = int(total_bases/2)
    cnt = 0
    i = 0
    while i < len(l_length):
        #print("%d" % l_length[i])
        cnt += l_length[i]
        if cnt > threshold:
            print("N50 of contig = %d" % l_length[i])
            i = len(l_length)
        i += 1
    return


def main(argv):
    infile = ""
    outfile = ""
    length = 0

    try:
        opts, args = getopt.getopt(argv,"hi:o:l:")
    except getopt.GetoptError:
        Usage()
        sys.exit(1)

    for opt, arg in opts:
        if opt == '-h':
            Usage()
            sys.exit()
        elif opt in "-i":
            infile = arg
        elif opt in "-l":
            length = int(arg)
        elif opt in "-o":
            outfile = arg

    if not infile:
        Usage()
        sys.exit(2)
    elif not os.path.isfile(infile):
        print("Error: input file(%s) is not existed" % infile)
        Usage()
        sys.exit(3)

    #Main Function
    sga_assemble_filter(infile, outfile, length)

    return

if __name__ == '__main__':
    main(sys.argv[1:])
