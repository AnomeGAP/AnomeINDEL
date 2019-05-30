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
# @file    fix.py
#
# @brief   A program to fix data
#
# @author  A-Tsai Su(chungtsai_su@anome.com.tw)
#
# @date    2018/11/23
#
# @version 1.0
#
# @remark
#

import sys
import getopt
import os.path
import re


def usage():
    """ Description: Program Usage
        Argument:    NONE
        Return:	     NONE
    """
    print("fix.py -i <input file> ")
    print("Argument:")
    print("\t-h: Usage")
    print("\t-i: <input file>")
    print("Usage:")
    print("\tpython ./fix.py -i ./ED1.txt > ED.txt")

    return


def fix(ifn1):

    total = 0
    print(ifn1)
    ifd = open(ifn1, "r")
    # Read sequence and segment them with length of slength
    for line in ifd:
        if len(line) <= 1:
            # print "Skip comment: %s" % (line)
            continue

        print("%s" % line)

    ifd.close()

    return


def main(argv):
    infile1 = ""

    try:
        opts, args = getopt.getopt(argv, "hi:j:")
    except getopt.GetoptError:
        usage()
        sys.exit(1)

    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit()
        elif opt in "-i":
            infile1 = arg
    if not infile1 :
        usage()
        sys.exit(2)
    elif not os.path.isfile(infile1):
        print("Error: input file(%s) is not existed" % infile1)
        usage()
        sys.exit(3)

    # Main Function
    fix(infile1)

    return


if __name__ == '__main__':
    main(sys.argv[1:])
