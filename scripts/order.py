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
# @file    CliqueStatistic.py
#
# @brief   A program to sort the order of RS and RP
#
# @author  A-Tsai Su(chungtsai_su@anome.com.tw)
#
# @date    2017/11/05
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
    """ Description: Program Usage
        Argument:    NONE
        Return:	     NONE
    """
    print("order.py -i <input file> ")
    print("Argument:")
    print("\t-h: Usage")
    print("\t-i: <input file>")
    print("Usage:")
    print("\tpython ~/src/github/AnomeINDEL/scripts/order.py -i ~/data/edges.log")

    return


def order(ifn1):

    total = 0
    #print(ifn1)
    ifd = open(ifn1, "r")
    # Read sequence and segment them with length of slength
    for line in ifd:
        if re.match("^#", line):
            # print "Skip comment: %s" % (line)
            continue

        line = line.strip()
        items = line.split("\t")
        elms = items[1].split(" ")
        if int(elms[0]) > int(elms[1]):
            elms[0],elms[1] = elms[1], elms[0]
            elms[2],elms[5] = elms[5], elms[2]
            elms[3],elms[6] = elms[6], elms[3]
            elms[4],elms[7] = elms[7], elms[4]
        print("%s\t%s" %(items[0], " ".join(elms)))
        total += 1

    ifd.close()

    # output
    #print("Number of edges = %ld" % total)

    return


def main(argv):
    infile1 = ""

    try:
        opts, args = getopt.getopt(argv, "hi:")
    except getopt.GetoptError:
        Usage()
        sys.exit(1)

    for opt, arg in opts:
        if opt == '-h':
            Usage()
            sys.exit()
        elif opt in "-i":
            infile1 = arg
    if not infile1 :
        Usage()
        sys.exit(2)
    elif not os.path.isfile(infile1):
        print("Error: input file(%s) is not existed" % infile1)
        Usage()
        sys.exit(3)

    # Main Function
    order(infile1)

    return


if __name__ == '__main__':
    main(sys.argv[1:])
