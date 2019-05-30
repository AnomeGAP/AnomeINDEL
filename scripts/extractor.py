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
    print("extractor.py -i <Line-based file>")
    print("Argument:")
    print("\t-h: Usage")
    print("\t-i: <Line-based file>")
    print("Usage:")
    print("\t./extractor.py -i test.data")

    return
def Extractor(ifn):

    ifd = open(ifn, "r")
    ##Read sequence and segment them with length of slength
    for line in ifd:
        if re.match("^#", line):
            #print "Skip comment: %s" % (line)
            continue
        if len(line) < 2:
            continue

        line = line.strip()
        items = line.split(',')
        #print(items[7]) #INFO
        vrecord = items[1].split(':')
        vitem = items[2].split(':')
        vtime = items[3].split(':')
        print("%s\t%s\t%s" %(vrecord[1], vitem[1], vtime[1].rstrip("s")))
    ifd.close()

    return

def main(argv):
    infile = ""


    try:
        opts, args = getopt.getopt(argv,"hi:")
    except getopt.GetoptError:
        Usage()
        sys.exit(1)

    for opt, arg in opts:
        if opt == '-h':
            Usage()
            sys.exit()
        elif opt in ("-i"):
            infile = arg

    if not infile:
        Usage()
        sys.exit(2)
    elif (os.path.isfile(infile) == False):
        print("Error: input file(%s) is not existed" % (infile))
        Usage()
        sys.exit(3)

    #Main Function
    Extractor(infile)

    return

if __name__ == '__main__':
    main(sys.argv[1:])
