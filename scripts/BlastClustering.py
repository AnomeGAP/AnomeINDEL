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
# @file    BlastClustering.py
#
# @brief   An Genome Identification by using graph-based clustering
#
# @author  A-Tsai Su(chungtsai_su@atgenomix.com)
#
# @date    2018/09/26
#
# @version 1.0
#
# @remark
#

import sys
import getopt
import os.path
import random
from math import *
from collections import defaultdict

# Constant
UNKNOWN = "Unclassified"
SEPERATOR = ":"
REFERENCE_GENOME_PREFIX = "CR_"

# #Parameter setting
THRESHOLD = 5  # <= (topscore - score)
MIN_BITSCORE = 150
MIN_SIMILARITY = 96  # percentage
# THRESHOLD = 10  # <= (topscore - score)
# MIN_BITSCORE = 0
# MIN_SIMILARITY = 0  # percentage
MAX_INSERTION_SIZE = 500
PERCENTAGE_CONCURRENCE = 0.9
CONTAMINANT = ["CR_195", "CR_487", "CR_475"]


def _usage():
    """ Description: Program Usage
        Argument:    NONE
        Return:	     NONE
    """
    print(
        "BlastClustering.py -i <Alignment result for read1> -j <Alignment result for read2> "
        "-o <Genome Identification Output file> -m <Method>")
    print("Argument:")
    print("\t-h: Usage")
    print("\t-i: <Blast alignment file>")
    print("\t-j: <Blast alignment file>")
    print("\t-o: <Output file>")
    print("\t-m: <method> (Default: 1)")
    print("\t\tMethod 1 (Aggressive Diversity) - first record (priority : r1 > r2)")
    print("\t\tMethod 2 (Conservative Diversity) - first properly paired record (priority: r1)")
    print("\t\tMethod 3 (Aggressive Focus) - all high score record and then clustering")
    print("\t\tMethod 4 (Conservative Focus) - all high score properly paired records and the clustering")

    print("Usage:")
    print("\tpython3 ./BlastClustering.py -i ../result/Hello_World_R1.tsv "
          "-j ../result/Hello_World_R2.tsv -m 1 -o ../result/Hello_World")
    print("\tpython3 ./BlastClustering.py -i ../result/Hello_World_R1.tsv.G1 "
          "-j ../result/Hello_World_R2.tsv.G1 -m 1 -o ../result/Hello_World")
    return


def scoring(alength, mismatches, gaps, bitscore):
    # score = alength - mismatches - gaps * 2 + bitscore - bitscore
    score = bitscore
    return score


def read_blast(fn):
    name = ""
    hits = 0
    nohits = 0
    cnt = 0
    (topscore, score) = (0.0, 0.0)
    candidates = set()
    h_ref = defaultdict(int)
    h_relation = defaultdict(int)
    (both_hit, human_hit, cr_hit, human_only, cr_only, total) = (0, 0, 0, 0, 0, 0)

    # # BLASTN 2.7.1+
    # # Query: Example.1/1
    # # Database: ../data/ALL_mask.db
    # # Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end,
    # #         s. start, s. end, evalue, bit score
    # # 7 hits found
    # Example.1/1     CR_471_Contig_0 100.000 151     0       0       1       151     1076722 1076872 1.34e-73       279
    # Example.1/1     CR_466_Contig_0 100.000 151     0       0       1       151     1088294 1088444 1.34e-73       279
    # Example.1/1     CR_411_Contig_0 100.000 151     0       0       1       151     1026010 1025860 1.34e-73       279
    # Example.1/1     CR_369_Contig_0 100.000 151     0       0       1       151     441918  442068  1.34e-73       279
    # Example.1/1     CR_303_Contig_1 100.000 151     0       0       1       151     783717  783567  1.34e-73       279
    # Example.1/1     CR_276_Contig_0 100.000 151     0       0       1       151     1760748 1760898 1.34e-73       279
    # Example.1/1     CR_203_Contig_1 100.000 151     0       0       1       151     413778  413628  1.34e-73       279

    ifd = open(fn, "r")
    for line in ifd:
        if line.startswith("#"):
            if line.startswith("# Query:"):
                if hits != cnt:
                    print("ERROR: %s has %d hits (expect %d hits)" % (name, cnt, hits))
                lc = list(candidates)
                for i in range(len(lc)):
                    candidates.add(i)
                    for j in range(len(lc)):
                        if i != j:
                            h_relation[lc[i] + ":" + lc[j]] += 1
                            h_relation[lc[j] + ":" + lc[i]] += 1
                if human_hit == 1 and cr_hit == 1:
                    both_hit += 1
                elif human_hit == 1 and cr_hit == 0:
                    human_only += 1
                elif human_hit == 0 and cr_hit == 1:
                    cr_only += 1
                items = list(line.strip().split(" "))
                name = items[2]
                # print("name = %s" % name)
            elif line.endswith("found\n"):
                items = list(line.strip().split(" "))
                hits = int(items[1])
                # print("hits= %d" % hits)
                if hits == 0:
                    nohits += 1
                    # print("%s has %d hits" % (name, hits))
                cnt = 0
                candidates.clear()
                human_hit = 0
                cr_hit = 0
                total += 1
        else:
            items = list(line.strip().split("\t"))
            if items[0] != name:
                print("ERROR: %s : %s" % (name, line))
            elif float(items[11]) >= MIN_BITSCORE and float(items[2]) >= MIN_SIMILARITY:
                score = scoring(int(items[3]), int(items[4]), int(items[5]), float(items[11]))
                if cnt == 0:
                    topscore = score
                if topscore - score <= THRESHOLD:
                    if items[1].startswith(REFERENCE_GENOME_PREFIX):
                        lc = items[1].split("_")
                        contig_name = lc[0] + "_" + lc[1]
                        h_ref[contig_name] += 1
                        candidates.add(contig_name)
                        cr_hit = 1
                    else:
                        human_hit = 1
            cnt += 1
    lc = list(candidates)
    for i in range(len(lc)):
        candidates.add(i)
        for j in range(len(lc)):
            if i != j:
                h_relation[lc[i] + ":" + lc[j]] += 1
                h_relation[lc[j] + ":" + lc[i]] += 1

    ifd.close()
    print("%s has %d reads" % (fn, total))
    print("%s has %d reads with no hits" % (fn, nohits))
    print("%s has %d reads both hit HUMAN and CR" % (fn, both_hit))
    print("%s has %d reads hit HUMAN only but not CR" % (fn, human_only))
    print("%s has %d reads hit CR only but not HUMAN" % (fn, cr_only))
    return h_ref, h_relation


def read_blast_by_aggressive(fn, fn2, h_mapping):
    name = ""
    cnt = 0
    h_read = defaultdict(str)

    ifd = open(fn2, "r")
    for line in ifd:
        if line.startswith("#"):
            if line.startswith("# Query:"):
                items = list(line.strip().split(" "))
                name = items[2].split("/")[0]
                h_read[name] = UNKNOWN
                # print("name = %s" % name)
                cnt = 0
        else:
            items = list(line.strip().split("\t"))
            if not items[0].startswith(name):
                print("ERROR: %s : %s" % (name, line))
            elif cnt == 0 and items[1].startswith(REFERENCE_GENOME_PREFIX) and \
                    (float(items[11]) >= MIN_BITSCORE and float(items[2]) >= MIN_SIMILARITY):
                lc = items[1].split("_")
                contig_name = lc[0] + "_" + lc[1]
                if bool(h_mapping):
                    h_read[name] = h_mapping[contig_name]
                else:
                    h_read[name] = contig_name
            if random.randrange(2):
                cnt += 1

    ifd.close()

    ifd = open(fn, "r")
    for line in ifd:
        if line.startswith("#"):
            if line.startswith("# Query:"):
                items = list(line.strip().split(" "))
                name = items[2].split("/")[0]
                cnt = 0
        else:
            items = list(line.strip().split("\t"))
            if not items[0].startswith(name):
                print("ERROR: %s : %s" % (name, line))
            elif cnt == 0 and items[1].startswith(REFERENCE_GENOME_PREFIX) and \
                (float(items[11]) >= MIN_BITSCORE and float(items[2]) >= MIN_SIMILARITY):
                lc = items[1].split("_")
                contig_name = lc[0] + "_" + lc[1]
                if bool(h_mapping):
                    h_read[name] = h_mapping[contig_name]
                else:
                    h_read[name] = contig_name
            if random.randrange(2):
                cnt += 1

    ifd.close()

    return h_read


def read_blast_by_conservative(fn, fn2, h_mapping):
    name = ""
    cnt = 0
    (topscore, score) = (0.0, 0.0)
    h_read = defaultdict(str)
    h_paired = defaultdict(str)

    ifd = open(fn2, "r")
    for line in ifd:
        if line.startswith("#"):
            if line.startswith("# Query:"):
                items = list(line.strip().split(" "))
                name = items[2].split("/")[0]
                cnt = 0
        else:
            items = list(line.strip().split("\t"))
            if not items[0].startswith(name):
                print("ERROR: %s : %s" % (name, line))
            elif items[1].startswith(REFERENCE_GENOME_PREFIX) and \
                    (float(items[11]) >= MIN_BITSCORE and float(items[2]) >= MIN_SIMILARITY):
                score = scoring(int(items[3]), int(items[4]), int(items[5]), float(items[11]))
                if cnt == 0:
                    topscore = score
                if topscore - score <= THRESHOLD:
                    lc = items[1].split("_")
                    contig_name = lc[0] + "_" + lc[1]
                    h_paired[name + SEPERATOR + contig_name] = int(items[9])
            if random.randrange(2):
                cnt += 1

    ifd.close()

    ifd = open(fn, "r")
    for line in ifd:
        if line.startswith("#"):
            if line.startswith("# Query:"):
                items = list(line.strip().split(" "))
                name = items[2].split("/")[0]
                h_read[name] = UNKNOWN
                # print("name = %s" % name)
                cnt = 0
                topscore = 0
        else:
            items = list(line.strip().split("\t"))
            if not items[0].startswith(name):
                print("ERROR: %s : %s" % (name, line))
            elif cnt == 0 and (float(items[11]) >= MIN_BITSCORE and float(items[2]) >= MIN_SIMILARITY):
                score = scoring(int(items[3]), int(items[4]), int(items[5]), float(items[11]))
                if topscore == 0:
                    topscore = score
                if topscore - score <= THRESHOLD and items[1].startswith(REFERENCE_GENOME_PREFIX):
                    lc = items[1].split("_")
                    contig_name = lc[0] + "_" + lc[1]
                    if name + SEPERATOR + contig_name in h_paired \
                            and abs(h_paired[name + SEPERATOR + contig_name] - int(items[9])) <= MAX_INSERTION_SIZE:
                        if bool(h_mapping):
                            h_read[name] = h_mapping[contig_name]
                        else:
                            h_read[name] = contig_name
                        # print("name = %s\t%s\t%d\t%d" % (name, contig_name, int(items[9]),
                        #                                  h_mapping[name + SEPERATOR + contig_name]))
                        if random.randrange(2):
                            cnt += 1

    ifd.close()

    return h_read


def output(h_read, ifn, ofn):
    h_cnt = defaultdict(int)
    total = 0
    unknown = 0
    ofd = open(ofn + "_Reads_Identification.tsv", "w")
    # write header
    ofd.write("readID\trefID\tscore\n")
    ifd = open(ifn, "r")
    for line in ifd:
        if line.startswith("# Query:"):
            items = list(line.strip().split(" "))
            name = items[2].split("/")[0]
            # print("name = %s\t%d\t%d" % (name, total, unknown))
            if h_read[name] in str(CONTAMINANT):
                ofd.write("%s\t%s\t%f\n" % (items[2], UNKNOWN, 0.5))
                h_cnt[UNKNOWN] += 1
                unknown += 1
            else:
                ofd.write("%s\t%s\t%f\n" % (items[2], h_read[name], 0.5))
                h_cnt[h_read[name]] += 1
                if h_read[name] == UNKNOWN:
                    unknown += 1
            total += 1
    ifd.close()
    ofd.close()

    max_confidence = 0.01
    for i in range(518):
        ref = "CR_%d" % i
        if h_cnt[ref] / (total - unknown) > max_confidence:
            max_confidence = h_cnt[ref] / (total - unknown)

    gifd = open(ofn + "_Genome_Identification.tsv", "w")
    gqfd = open(ofn + "_Genome_Quantification.tsv", "w")
    for i in range(518):
        ref = "CR_%d" % i
        gifd.write("%s\t%f\n" % (ref, sqrt(h_cnt[ref] / (total - unknown) / max_confidence)))
        gqfd.write("%s\t%f\n" % (ref, h_cnt[ref] / total))
    gqfd.write("%s\t%f\n" % (UNKNOWN, unknown / total))
    gqfd.close()
    gifd.close()
    return


def blastclustering(ifn, i2fn, method, ofn):

    # Strategy:
    # Aggressive/Conservative - more UNKNOWN for Conservative (Properly paired)
    # Diversity/Focus - more species for diversity
    # Method 1 (Aggressive Diversity) - first record (priority : r1 > r2)
    # Method 2 (Conservative Diversity) - first properly paired record (priority: r1)
    # Method 3 (Aggressive Focus) - all high score record and then clustering
    # Method 4 (Conservative Focus) - all high score properly paired records and the clustering
    if method == 1:
        h_read = read_blast_by_aggressive(ifn, i2fn, {})
        output(h_read, ifn, ofn)
    elif method == 2:
        h_read = read_blast_by_conservative(ifn, i2fn, {})
        output(h_read, ifn, ofn)
    elif method == 3 or method == 4:
        h_mapping = defaultdict(str)
        (h_ref, h_relation) = read_blast(ifn)
        (h_ref2, h_relation2) = read_blast(i2fn)
        for k in h_ref2.keys():
            h_ref[k] += h_ref2[k]
        for k in h_relation2.keys():
            h_relation[k] += h_relation2[k]

        for key, value in sorted(h_ref.items(), key=lambda item: (item[1], item[0]), reverse=True):
            if value >= 300:
                print("%s: %s" % (key, value))

        for key, value in sorted(h_relation.items(), key=lambda item: (item[1], item[0]), reverse=True):
            if value >= 800:
                print("%s: %s" % (key, value))
            (id1, id2) = key.strip().split(":")
            if h_ref[id1] >= h_ref[id2]:
                if id1 not in h_mapping:
                    if id2 not in h_mapping:
                        h_mapping[id1] = id1
                        if value >= PERCENTAGE_CONCURRENCE * h_ref[id2]:
                            h_mapping[id2] = h_mapping[id1]
                        else:
                            h_mapping[id2] = id2
                    elif h_ref[id1] <= h_ref[h_mapping[id2]]:
                        h_mapping[id1] = h_mapping[id2]
                else:
                    if id2 not in h_mapping:
                        h_mapping[id2] = h_mapping[id1]
                    elif h_ref[h_mapping[id2]] > h_ref[h_mapping[id1]]:  # and value >= 400:
                        for k in h_mapping.keys():
                            if h_mapping[k] == h_mapping[id1] and k != id1:
                                print("WARN: recursive update %s: %s=>%s" % (k, h_mapping[k], h_mapping[id2]))
                                h_mapping[k] = h_mapping[id2]
                        h_mapping[id1] = h_mapping[id2]
        for key in h_mapping.keys():
            print("%s: %s" % (key, h_mapping[key]))
        if method == 3:
            h_read = read_blast_by_aggressive(ifn, i2fn, h_mapping)
        else:
            h_read = read_blast_by_conservative(ifn, i2fn, h_mapping)

        output(h_read, ifn, ofn)

    return


def main(argv):
    infile = ""
    in2file = ""
    outfile = ""
    method = 1

    try:
        opts, args = getopt.getopt(argv, "hi:j:o:m:")
    except getopt.GetoptError:
        _usage()
        sys.exit(1)

    for opt, arg in opts:
        if opt == '-h':
            _usage()
            sys.exit()
        elif opt in "-i":
            infile = arg
            outfile = infile + "_output.tsv"
        elif opt in "-j":
            in2file = arg
        elif opt in "-m":
            method = int(arg)
        elif opt in "-o":
            outfile = arg
    if not infile or not os.path.isfile(infile):
        print("Error: input file(%s) is not existed" % infile)
        _usage()
        sys.exit(2)
    elif not in2file or not os.path.isfile(in2file):
        print("Error: input file(%s) is not existed" % in2file)
        _usage()
        sys.exit(3)

    # Main Function
    blastclustering(infile, in2file, method, outfile)

    return


if __name__ == '__main__':
    main(sys.argv[1:])
