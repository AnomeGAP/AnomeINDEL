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
# @file    EC_Benchmarking.py
#
# @brief   An evaluation tool for Error Correction TSV data
#
# @author  A-Tsai Su(chungtsai_su@anome.com.tw)
#
# @date    2018/09/13
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
from enum import IntEnum


# Parameter setting
class CLASS(IntEnum):
    SINGLE = 0
    MULTIPLE = 1
    END = 2


class TYPE(IntEnum):
    SNP = 0
    INS = 1
    DEL = 2
    END = 3


class ANSWER(IntEnum):
    TP = 0
    FP = 1
    FN = 2
    END = 3


def _usage():
    """ Description: Program Usage
        Argument:    NONE
        Return:	     NONE
    """
    print(
        "EC_Benchmarking.py -i <Answer TSV file> -j <Predicted TSV file>")
    print("Argument:")
    print("\t-h: Usage")
    print("\t-i: <FASTQ file>")
    print("\t-j: <FASTQ file>")

    print("Usage:")
    print("\tpython3 ./EC_Benchmarking.py -i new_chr22_error.tsv -j new_chr22_sga_error.tsv")

    return


def cal_ec(ifn, i2fn):
    h_cnt = defaultdict(int)
    h_type = defaultdict(int)
    h_ref = defaultdict(str)
    h_error = defaultdict(str)
    h_total = [[0 for y in range(TYPE.END)] for x in range(CLASS.END)]
    h_answer = [[0 for y in range(ANSWER.END)] for x in range(CLASS.END*TYPE.END)]
    total = 0

    ifd = open(ifn, "r")
    for line in ifd:
        if re.match("#", line):
            continue
        items = list(line.strip().split("\t"))
        idx = int(items[0])
        pos = int(items[1])
        ref = items[2]
        err = items[3]
        h_cnt[idx] += 1
        key = "%d-%d" % (idx, pos)
        h_ref[key] = ref
        h_error[key] = err
        if ref == '-':
            h_type[key] = TYPE.INS
        elif err == '-':
            h_type[key] = TYPE.DEL
        else:
            h_type[key] = TYPE.SNP
        total += 1
    ifd.close()

    for k in h_type.keys():
        items = str(k).split("-")
        idx = int(items[0])
        if h_cnt[idx] == 1:
            h_total[CLASS.SINGLE][h_type[k]] += 1
        else:
            h_total[CLASS.MULTIPLE][h_type[k]] += 1

    print("CLASS\tTYPE\tAmount")
    for i in range(CLASS.END):
        for j in range(TYPE.END):
            print("%d\t%d\t%d" % (i, j, h_total[i][j]))
    print("Total=%d" % total)

    total = 0
    num_corrected = 0
    i2fd = open(i2fn, "r")
    for line in i2fd:
        if re.match("#", line):
            continue
        items = list(line.strip().split("\t"))
        idx = int(items[0])
        pos = int(items[1])
        err = items[2]
        ref = items[3]
        key = "%d-%d" % (idx, pos)
        if key in h_error and h_error[key] != err:
            print("ERROR!!!! ID:%d POS:%d ERROR:%s != %s" % (idx, pos, err, h_error[key]))
        i = (CLASS.MULTIPLE if h_cnt[idx] != 1 else CLASS.SINGLE) * TYPE.END + h_type[key]
        if key in h_error and h_ref[key] == ref:
            h_answer[i][ANSWER.TP] += 1
            num_corrected += 1
        else:
            h_answer[i][ANSWER.FP] += 1
        total += 1
    i2fd.close()

    print("CLASS\tTYPE\tANSWER\tAmount")
    TP = 0
    FP = 0
    FN = 0
    for i in range(CLASS.END):
        for j in range(TYPE.END):
            h_answer[i*TYPE.END+j][ANSWER.FN] = h_total[i][j] - h_answer[i*TYPE.END+j][ANSWER.TP]
            for k in range(ANSWER.END):
                print("%d\t%d\t%d\t%d" % (i, j, k,  h_answer[i*TYPE.END+j][k]))
            TP += h_answer[i*TYPE.END+j][ANSWER.TP]
            FP += h_answer[i*TYPE.END+j][ANSWER.FP]
            FN += h_answer[i*TYPE.END+j][ANSWER.FN]
    print("Total Prediction=%d" % total)
    print("Precision = %5.2f" % float(TP*100 / (TP + FP)))
    print("Sensitivity = %5.2f" % float(TP*100 / (TP + FN)))
    print("TP=%d" % TP)
    print("FP=%d" % FP)
    print("FN=%d" % FN)

    return


def main(argv):
    infile = ""
    in2file = ""

    try:
        opts, args = getopt.getopt(argv, "hi:j:")
    except getopt.GetoptError:
        _usage()
        sys.exit(1)

    for opt, arg in opts:
        if opt == '-h':
            _usage()
            sys.exit()
        elif opt in "-i":
            infile = arg
        elif opt in "-j":
            in2file = arg
    if not infile or not os.path.isfile(infile):
        print("Error: input file(%s) is not existed" % infile)
        _usage()
        sys.exit(2)

    # Main Function
    cal_ec(infile, in2file)

    return


if __name__ == '__main__':
    main(sys.argv[1:])
