#!/bin/env python3
#
# @note Copyright (C) 2018, Atgenomix Incorporated. All Rights Reserved.
#       This program is an unpublished copyrighted work which is proprietary to
#       Atgenomix Incorporated and contains confidential information that is not to
#       be reproduced or disclosed to any other person or entity without prior
#       written consent from Atgenomix, Inc. in each and every instance.
#
# @warning Unauthorized reproduction of this program as well as unauthorized
#          preparation of derivative works based upon the program or distribution of
#          copies by sale, rental, lease or lending are violations of federal copyright
#          laws and state trade secret laws, punishable by civil and criminal penalties.
#
# @file    nregions_from_reference.py
#
# @brief   identify continuous N regions from reference genome
#
# @author  A-Tsai Su(chungtsai_su@atgenomix.com)
#
# @date    2019/02/26
#
# @version 1.0
#
# @remark
#

import sys
import getopt
import os.path


# Parameter setting


def _usage():
    """ Description: Program Usage
        Argument:    NONE
        Return:	     NONE
    """
    print(
        "nregions_from_reference.py -i <Reference Genome> -o <Output Interval File>")
    print("Argument:")
    print("\t-h: Usage")
    print("\t-i: <Reference Genome>")
    print("\t-o: <Output Interval file>")

    print("Usage:")
    print("\tpython ./nregions_from_reference.py -i /seqslab/mnt/reference/38/HG/ref.fa -o ~/hg38.fa.n.interval.tsv")

    return


def nregions_from_reference(ifn, ofn):
    ofd = open(ofn, "w")
    start = 0
    end = 0
    pos = 0
    name = ""
    num_n = 0
    num_segments = 0

    ifd = open(ifn, "r")
    for line in ifd:
        if line.startswith(">"):
            if start > 0:
                ofd.write("%s\t%d\t%d\t%d\n" % (name, start, pos, pos - start + 1))
                num_n += pos - start + 1
                num_segments += 1
            if name != "":
                print("%s\t%d\t%d\t%d" % (name, num_n, num_segments, pos))
            name = line.strip().split(">")[1].split(" ")[0]
            pos = 0
            num_n = num_segments = 0
            start = 0
        else:
            read = []
            read += line.strip()
            for i in range(len(read)):
                if start == 0 and read[i] == 'N':
                    start = pos + i + 1
                elif start > 0 and read[i] != 'N':
                    end = pos + i
                    ofd.write("%s\t%d\t%d\t%d\n" % (name, start, end, end-start+1))
                    num_n += end-start+1
                    num_segments += 1
                    start = 0
            pos += len(read)
    if start > 0:
        ofd.write("%s\t%d\t%d\t%d\n" % (name, start, pos, pos - start + 1))
        num_n += pos - start + 1
        num_segments += 1
    print("%s\t%d\t%d\t%d" % (name, num_n, num_segments, pos))

    ifd.close()
    ofd.close()

    return


def main(argv):
    infn = ""
    outfn = ""

    try:
        opts, args = getopt.getopt(argv, "hi:o:")
    except getopt.GetoptError:
        _usage()
        sys.exit(1)

    for opt, arg in opts:
        if opt == '-h':
            _usage()
            sys.exit()
        elif opt in "-i":
            infn = arg
            outfn = infn + ".n.interval.tsv"
        elif opt in "-o":
            outfn = arg
    if not infn or not os.path.isfile(infn):
        print("Error: input file(%s) is not existed" % infn)
        _usage()
        sys.exit(2)

    # Main Function
    nregions_from_reference(infn, outfn)

    return


if __name__ == '__main__':
    main(sys.argv[1:])
