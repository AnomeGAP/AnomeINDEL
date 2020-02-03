#!/bin/env python
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
# @file    FASTQFormat.py
#
# @brief   Reformat FASTQ Header
#
# @author  Chung-Tsai Su(chungtsai_su@atgenomix.com)
#
# @date    2019/07/20
#
# @version 1.0
#
# @remark
#

import sys
import getopt
import os.path
import re
from collections import defaultdict
import fileinput

# Default Parameter


def usage():
    print("zcat <Input FASTQ> | FASTQFormat.py -o <Output FASTQ>")
    print("Argument:")
    print("\t-h: Usage")
    print("Usage:")
    print("\tzcat ../data/NA19240/fastq/SRR7782669_pass_1.fastq.gz | python ./FASTQFormat.py -i 1 -o ../data/NA19240/fastq/SRR7782669_new_r1.fq")
    print("\tzcat ../data/NA19240/fastq/SRR7782669_pass_2.fastq.gz | python ./FASTQFormat.py -i 2 -o ../data/NA19240/fastq/SRR7782669_new_r2.fq")

    return


def analyzer(id, ofn):
    i = 0
    idx = 0
    ofd = open(ofn, "w")
    for line in sys.stdin:
        if i % 4 == 0:
            idx += 1
            ofd.write("@CC%d:%d\n" % (idx, id))
        elif i % 4 == 2:
            ofd.write("+\n")
        elif i % 4 == 3 and str(line).startswith("@"):
            # temp = "F" * len(line.strip())
            # ofd.write("%s\n" % temp)
            temp = "F" * len(line.strip())
            ofd.write("A%s\n" % line.strip()[1:])
        else:
            ofd.write(line)
        i += 1

    ofd.close()

    return


def main(argv):
    ofile = "output.fq"
    id = 1

    try:
        opts, args = getopt.getopt(argv, "hi:o:")
    except getopt.GetoptError:
        usage()
        sys.exit(1)

    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit()
        elif opt == '-i':
            id = int(arg)
        elif opt == '-o':
            ofile = arg

    analyzer(id, ofile)

    return


if __name__ == '__main__':
    main(sys.argv[1:])
