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
import re

MIN_COVERAGE=10
MIN_QUALITY=10

def Usage():
    ''' Description: Program Usage
        Argument:    NONE
        Return:	     NONE
    '''
    print("VCFQualityFilter.py -i <VCF file> -c <integer> -q <integer>")
    print("Argument:")
    print("\t-h: Usage")
    print("\t-i: <VCF file>")
    print("\t-c: <Minimal number of coverage(read depth)")
    print("\t-q: <Minimal number of quality")
    print("Usage:")
    print("\t./VCFQualityFilter.py -i TSC181K_11132015_bwamem.del.vcf -c 10 -q 40")

    return


def VCFQualityFilter(ifn, min_cov, min_qual):

    ifd = open(ifn, "r")
    ##Read sequence and segment them with length of slength
    for line in ifd:
        if re.match("^#", line):
            #print "Skip comment: %s" % (line)
            continue

        line = line.strip()
        items = line.split('\t')
        #print(items[7]) #INFO
        vlist = items[7].split(';')
        cov = qual = 0
        for vitem in vlist:
            values = vitem.split('=', 1)
            if ( values[0] == 'MAPQ'):
                qual = int(values[1])
            elif (values[0] == 'PE'):
                cov = int(values[1])
        if (cov >= min_cov) and (qual >= min_qual):
            print("%s\t%s\t%s" %(cov,qual,line))
    ifd.close()

    return

def main(argv):
    infile = ""
    cov = MIN_COVERAGE
    qual = MIN_QUALITY

    try:
        opts, args = getopt.getopt(argv,"hi:c:q:")
    except getopt.GetoptError:
        Usage()
        sys.exit(1)

    for opt, arg in opts:
        if opt == '-h':
            Usage()
            sys.exit()
        elif opt in ("-i"):
            infile = arg
        elif opt in ("-c"):
            cov = int(arg)
        elif opt in ("-q"):
            qual = int(arg)
    if not infile:
        Usage()
        sys.exit(2)
    elif (os.path.isfile(infile) == False): 
        print("Error: input file(%s) is not existed" % (infile))
        Usage()
        sys.exit(3)

    #Main Function
    VCFQualityFilter(infile, cov, qual)

    return

if __name__ == '__main__':
    main(sys.argv[1:])
