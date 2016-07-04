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
# @brief   A program to statistic clique data
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


def Usage():
    """ Description: Program Usage
        Argument:    NONE
        Return:	     NONE
    """
    print("CliqueStatistic.py -i <Clique file> -c <integer>")
    print("Argument:")
    print("\t-h: Usage")
    print("\t-i: <Clique file>")
    print("Usage:")
    print("\tpython ~/src/github/AnomeINDEL/scripts/CliqueStatistic.py -i ~/src/clique-tf/output/TF2SITE.complete.out")

    return

MAX_COUNT = 20


def CliqueStatistic(ifn):
    _num_cliques = 0
    _max_lnode = 0
    _lnode_name = ""
    _max_neighbors = 0
    _max_neighbor_name = ""
    _max_rnode = 0
    _rnode_name = ""
    _hcount = [[0 for x in range(MAX_COUNT+1)] for x in range(MAX_COUNT+1)]

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
        for item in lLeft:
            if item == "":
                break
            val = item.split(',', 1)
            if int(val[1]) > _max_neighbors:
                _max_neighbors = int(val[1])
                _max_neighbor_name = val[0]

        lRight = items[1].split(';')
        num_right = len(lRight) - 1

        if num_left > _max_lnode:
            _max_lnode = num_left
            _lnode_name = line

        if num_right > _max_rnode:
            _max_rnode = num_right
            _rnode_name = line

        if num_left > MAX_COUNT:
            num_left = MAX_COUNT
        if num_right > MAX_COUNT:
            num_right = MAX_COUNT

        _hcount[num_left][num_right] += 1
        _num_cliques += 1

    ifd.close()

    # output
    print("Total cliques = %d" % _num_cliques)
    print("Maximal neighbors=%d %s" % (_max_neighbors,_max_neighbor_name))
    print("Maximal left node = %d %s" % (_max_lnode, _lnode_name))
    print("Maximal right node = %d %s" % (_max_rnode, _rnode_name))
    for i in range(MAX_COUNT):
            print("%s" % "\t".join([str(x) for x in _hcount[i]]))

    return


def main(argv):
    infile = ""

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
            infile = arg
    if not infile:
        Usage()
        sys.exit(2)
    elif not os.path.isfile(infile):
        print("Error: input file(%s) is not existed" % (infile))
        Usage()
        sys.exit(3)

    # Main Function
    CliqueStatistic(infile)

    return


if __name__ == '__main__':
    main(sys.argv[1:])
