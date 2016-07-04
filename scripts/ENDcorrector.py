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
# @file    ENDcorrector.py
#
# @brief   An example to correct END column for VCFGenerator
#
# @author  A-Tsai Su(chungtsai_su@anome.com.tw)
#
# @date    2016/02/29
#
# @version 1.0
#
# @remark
#

import sys
import getopt
import os.path
import re


def Usage():
    ''' Description: Program Usage
        Argument:    NONE
        Return:	     NONE
    '''
    print("VCFQualityFilter.py -i <VCF file> ")
    print("Argument:")
    print("\t-h: Usage")
    print("\t-i: <VCF file>")
    print("Usage:")
    print("\t./ENDcorrector.py -i /anomegap/2015summer/rawkey/vcfkey/2/vcf2.part0.vcf")

    return
def ENDcorrector(ifn):

    ifd = open(ifn, "r")
    ##Read sequence and segment them with length of slength
    for line in ifd:
        if re.match("^#", line):
            #print(line)
            continue

        #22      46720879        .       A       <CN0    60      PASS    END=0
        line = line.strip()
        items = line.split('\t')
        #print(items[7]) #INFO
        start = items[1]
        ref = items[3]
        vlist = items[7].split(';')
        for vitem in vlist:
            values = vitem.split('=', 1)
            if values[0] == 'END':
                end = int(values[1])
        if (end == 0) or (end < start):
            if ref == "":
                end = int(start)
            else:
                end = int(start) + len(ref) - 1
            print("%s\t%s\t%s\t%s\t%s\t%s\t%s\tEND=%d" % (items[0], items[1], items[2], items[3], items[4], items[5], items[6], end))

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
    ENDcorrector(infile)

    return

if __name__ == '__main__':
    main(sys.argv[1:])
