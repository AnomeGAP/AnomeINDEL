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
# @file    SeqGraph-Debug.py
#
# @brief   A program to figure out where the difference is between two String Graph
#
# @author  A-Tsai Su(chungtsai_su@atgenomix.com)
#
# @date    2017/12/07
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
    print("SeqGraph-Debug.py -i <VT file1> -j <ED file1> -k <VT file2> -l <ED file2>")
    print("Argument:")
    print("\t-h: Usage")
    print("\t-i: <VT file1>")
    print("\t-j: <ED file1>")
    print("\t-k: <VT file2>")
    print("\t-l: <ED file2>")
    print("Usage:")
    print("\tpython ~/src/github/AnomeINDEL/scripts/SeqGraph-Debug.py -i ~/data/VT1.log -j ~/data/ED1.log "
          "-k ~/data/VT2.log -l ~/data/ED2.log")

    return


def SeqGraph_Debug(VT1, ED1, VT2, ED2):

    cnt = 0
    hash_seq = {}
    ifd = open(VT1, "r")
    for line in ifd:
        if re.match("^#", line):
            # print "Skip comment: %s" % (line)
            continue
        line = line.strip()
        if len(line) <= 0:
            continue

        items = line.split("\t")
        hash_seq[items[2]] = items[1]
        #print("%s\t%s" % (items[2], items[1]))
        cnt += 1
        if (cnt % 1000000) == 0:
            print("VT1\t%d" % cnt)
    ifd.close()

    cnt = 0
    hash_vt = {}
    ifd = open(VT2, "r")
    for line in ifd:
        if re.match("^#", line):
            # print "Skip comment: %s" % (line)
            continue
        line = line.strip()
        if len(line) <= 0:
            continue

        items = line.split("\t")
        #print("%s\t%s" % (items[2], items[1]))
        hash_vt[hash_seq[items[2]]] = items[1]
        cnt += 1
        if (cnt % 1000000) == 0:
            print("VT2\t%d" % cnt)
    ifd.close()

    cnt = 0
    hash_ed1 = {}
    ifd = open(ED1, "r")
    for line in ifd:
        if re.match("^#", line):
            # print "Skip comment: %s" % (line)
            continue
        line = line.strip()
        if len(line) <= 0:
            continue

        val = line.split("\t")
        items = val[1].split(" ")
        #print("%s\t%s" % (items[0], items[1]))
        if items[0] in hash_ed1:
            hash_ed1[items[0]] += 1
        else:
            hash_ed1[items[0]] = 1
        if items[1] in hash_ed1:
            hash_ed1[items[1]] += 1
        else:
            hash_ed1[items[1]] = 1
        cnt += 1
        if (cnt % 1000000) == 0:
            print("ED1\t%d" % cnt)
    ifd.close()

    cnt = 0
    hash_ed2 = {}
    ifd = open(ED2, "r")
    for line in ifd:
        if re.match("^#", line):
            # print "Skip comment: %s" % (line)
            continue
        line = line.strip()
        if len(line) <= 0:
            continue

        val = line.split("\t")
        items = val[1].split(" ")
        #print("%s\t%s" % (items[0], items[1]))
        if items[0] in hash_ed2:
            hash_ed2[items[0]] += 1
        else:
            hash_ed2[items[0]] = 1
        if items[1] in hash_ed2:
            hash_ed2[items[1]] += 1
        else:
            hash_ed2[items[1]] = 1
        cnt += 1
        if (cnt % 1000000) == 0:
            print("ED2\t%d" % cnt)
    ifd.close()

    # output
    #for key in hash_vt.keys():
    #    print("%s\t%s" % (key, hash_vt[key]))
    ofd = open("debug.log", "w")

    for key in hash_ed1.keys():
        if hash_ed1[key] != hash_ed2[hash_vt[key]]:
            ofd.write("ERROR:%s\t%d\t%s\t%d\n" % (key, hash_ed1[key], hash_vt[key], hash_ed2[hash_vt[key]]))

    ofd.close()
    return


def main(argv):
    inVT1 = ""
    inVT2 = ""
    inED1 = ""
    inED2 = ""

    try:
        opts, args = getopt.getopt(argv, "hi:j:k:l:")
    except getopt.GetoptError:
        Usage()
        sys.exit(1)

    for opt, arg in opts:
        if opt == '-h':
            Usage()
            sys.exit()
        elif opt in "-i":
            inVT1 = arg
        elif opt in "-j":
            inED1 = arg
        elif opt in "-k":
            inVT2 = arg
        elif opt in "-l":
            inED2 = arg
    if not inVT1 or not inVT2 or not inED1 or not inED2:
        Usage()
        sys.exit(2)
    elif not os.path.isfile(inVT1):
        print("Error: input file(%s) is not existed" % inVT1)
        Usage()
        sys.exit(3)
    elif not os.path.isfile(inVT2):
        print("Error: input file(%s) is not existed" % inVT2)
        Usage()
        sys.exit(4)
    elif not os.path.isfile(inED1):
        print("Error: input file(%s) is not existed" % inED1)
        Usage()
        sys.exit(5)
    elif not os.path.isfile(inED2):
        print("Error: input file(%s) is not existed" % inED2)
        Usage()
        sys.exit(6)

    # Main Function
    SeqGraph_Debug(inVT1, inED1, inVT2, inED2)

    return


if __name__ == '__main__':
    main(sys.argv[1:])
