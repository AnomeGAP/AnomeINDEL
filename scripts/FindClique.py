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
# @file    FindClique.py
#
# @brief   A program to verify and filter clique data
#
# @author  A-Tsai Su(chungtsai_su@anome.com.tw)
#
# @date    2016/04/18
#
# @version 1.0
#
# @remark
#

import sys
import getopt
import os.path
import re

MAX_MISSING = 1
MIN_LEFT = 1
MIN_RIGHT = 1

def Usage():
    """ Description: Program Usage
        Argument:    NONE
        Return:	     NONE
    """
    print("FindClique.py -i <Clique file> -c <integer> -l <integer> -r <integer>")
    print("Argument:")
    print("\t-h: Usage")
    print("\t-i: <Clique file>")
    print("\t-c: <Maximal number of missing edge>")
    print("\t-r: <Minimal number of left nodes>")
    print("\t-c: <Minimal number of right nodes>")
    print("Usage:")
    print("\tpython ~/src/github/AnomeINDEL/scripts/FindClique.py -i TF2SITE.complete.out -c 1 -l 3 -r 2")

    return


def cliqueFilter(ifn, max_missing, min_left, min_right):
    ifd = open(ifn, "r")
    # Read sequence and segment them with length of slength
    for line in ifd:
        if re.match("^#", line):
            # print "Skip comment: %s" % (line)
            continue

        line = line.strip()
        items = line.split('\t')
        lLeft = items[0].split(";")
        num_left = len(lLeft) - 1
        if (num_left < min_left):
            continue

        lRight = items[1].split(';')
        if len(lRight) < min_right + 1:
            continue
        found = 1
        for item in lRight:
            if item == "":
                break
            val = item.split(',', 1)
            if int(val[1]) + max_missing < num_left:
                found = 0
                break
        if found == 1:
            print("%s" % line)
    ifd.close()

    return


def main(argv):
    infile = ""
    missing = MAX_MISSING
    min_left = MIN_LEFT
    min_right = MIN_RIGHT

    try:
        opts, args = getopt.getopt(argv, "hi:l:r:c:")
    except getopt.GetoptError:
        Usage()
        sys.exit(1)

    for opt, arg in opts:
        if opt == '-h':
            Usage()
            sys.exit()
        elif opt in "-i":
            infile = arg
        elif opt in "-c":
            missing = int(arg)
        elif opt in "-l":
            min_left = int(arg)
        elif opt in "-r":
            min_right = int(arg)
    if not infile:
        Usage()
        sys.exit(2)
    elif not os.path.isfile(infile):
        print("Error: input file(%s) is not existed" % (infile))
        Usage()
        sys.exit(3)

    # Main Function
    cliqueFilter(infile, missing, min_left, min_right)

    return


if __name__ == '__main__':
    main(sys.argv[1:])
