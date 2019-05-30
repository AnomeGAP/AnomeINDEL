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
# @file    sga-ed-extractByLayer.py
#
# @brief   An example to transform SGA ED record to graph record for Gephi
#
# @author  A-Tsai Su(chungtsai_su@anome.com.tw)
#
# @date    2018/05/16
#
# @version 1.0
#
# @remark
#

import sys
import getopt
import os.path


def usage():
    """ Description: Program Usage
        Argument:    NONE
        Return:	     NONE
    """
    print("sga-ed-extractByLayer.py -i <SGA ED file> -r <rid> -l <number of layer> -o <SGA ED visual file>")
    print("Argument:")
    print("\t-h: Usage")
    print("\t-i: <SGA ED file>")
    print("\t-r: read ID (XXX for RS of the first edge")
    print("\t-l: number of layers")
    print("\t-p: <previous SGA ED visual file>")
    print("\t-0: <SGA ED visual file>")
    print("Usage:")
    print("\tpython /sga-ed-extractByLayer.py -i m85.graphseq.ed.txt -r XXX -l 10 -p XXX.19.tsv -o XXX.tsv")

    return


def sga_ed_extractbylayer(ifn, rid, nlayer, prefile, ofn):
    hnext = {}
    if prefile != "":
        pfd = open(prefile, "r")
        for line in pfd:
            line = line.strip()
            items = line.split('\t')
            if items[0] != "Source":
                hnext[items[0].split('R')[0]] = 1
                hnext[items[1].split('R')[0]] = 1

    for i in range(0, nlayer):
        hbag = hnext
        hnext = {}
        cnt = 0
        ofd = open(ofn + "." + str(i), "w")
        ofd.write("Source\tTarget\tValue\n")
        ifd = open(ifn, "r")
        # ED    123 456 100 150 151
        for line in ifd:
            line = line.strip()
            items = line.split('\t')
            res = items[1].split(' ')
            if rid == 'XXX' and prefile == "":
                rid = res[0]
            hbag[rid] = 1
            # print("ADD %s\n" % res[0])
            # print("check %s %s\n" % (res[0], res[1]))
            if res[0] in hbag or res[1] in hbag:
                hnext[res[0]] = 1
                hnext[res[1]] = 1
                # print("add %s\n" % res[0])
                # print("add %s\n" % res[1])
                if res[2] == '0':
                    if res[5] == '0':
                        ofd.write("%sRC\t%s\t%d\n" % (res[0], res[1], abs(int(res[2]) - int(res[3])) + 1))
                        ofd.write("%sRC\t%s\t%d\n" % (res[1], res[0], abs(int(res[2]) - int(res[3])) + 1))
                    else:
                        ofd.write("%sRC\t%sRC\t%d\n" % (res[0], res[1], abs(int(res[2]) - int(res[3])) + 1))
                        ofd.write("%s\t%s\t%d\n" % (res[1], res[0], abs(int(res[2]) - int(res[3])) + 1))
                else:
                    if res[5] == '0':
                        ofd.write("%s\t%s\t%d\n" % (res[0], res[1], abs(int(res[2]) - int(res[3])) + 1))
                        ofd.write("%sC\t%sRC\t%d\n" % (res[1], res[0], abs(int(res[2]) - int(res[3])) + 1))
                    else:
                        ofd.write("%s\t%sRC\t%d\n" % (res[0], res[1], abs(int(res[2]) - int(res[3])) + 1))
                        ofd.write("%s\t%sRC\t%d\n" % (res[1], res[0], abs(int(res[2]) - int(res[3])) + 1))
                cnt += 2
        print("the %d-th round has %d edges\b" % (i, cnt))
        ifd.close()
        ofd.close()

    return


def main(argv):
    infile = ""
    outfile = ""
    prefile = ""
    rid = "XXX"
    nlayer = 1

    try:
        opts, args = getopt.getopt(argv, "hi:l:r:p:o:")
    except getopt.GetoptError:
        usage()
        sys.exit(1)

    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit()
        elif opt in "-i":
            infile = arg
        elif opt in "-l":
            nlayer = int(arg)
        elif opt in "-r":
            rid = arg
        elif opt in "-p":
            prefile = arg
        elif opt in "-o":
            outfile = arg

    if not infile:
        usage()
        sys.exit(2)
    elif not os.path.isfile(infile):
        print("Error: input file(%s) is not existed" % infile)
        usage()
        sys.exit(3)

    # Main Function
    sga_ed_extractbylayer(infile, rid, nlayer, prefile, outfile)

    return


if __name__ == '__main__':
    main(sys.argv[1:])
