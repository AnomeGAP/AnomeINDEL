#!/usr/bin/env python3
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
# @file    VCFComparison.py
#
# @brief   An example to compare VCF generated by SNP Array and other sources
#
# @author  A-Tsai Su(chungtsai_su@anome.com.tw)
#
# @date    2017/05/17
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
    ''' Description: Program Usage
        Argument:    NONE
        Return:	     NONE
    '''
    print("VCFComparison.py -i <Answer-based VCF file> -j <first VCF file> -k <second VCF file> [-p]")
    print("Argument:")
    print("\t-h: Usage")
    print("\t-i: <Answer-based VCF file>")
    print("\t-j: <first VCF file>")
    print("\t-k: <second VCF file>")
    print("\t-p: check position only (skip REF/ALT)")
    print("Usage:")
    print("\t./VCFComparison.py -p -i snp_array.vcf -j seqslab.vcf -k multi-threading.vcf")

    return

def MakeKey(pos, prefix,  items):
    if pos == 1:
        key = prefix + items[0] + ":" + items[1]
    else:
        key = prefix + items[0] + ":" + items[1] + ":" + items[3] + ":" + items[4]
    return key

def VCFComparision(ifn, ifn1, ifn2, pos):
    answer = {}
    counter = 0

    ifd = open(ifn, "r")
    for line in ifd:
        if re.match("^#", line):
            #print "Skip comment: %s" % (line)
            continue

        line = line.strip()
        items = line.split('\t')

        if items[9] == "./.":
            continue

        key = MakeKey(pos, "chr",  items)       ## SNP array data has no "chr" prefix string

        #print(key) #INFO
        answer[key] = 1
        counter += 1
        if counter % 100000 == 0:
            print("%d %s" % (counter, key))

    ifd.close()
    print("%s has %d variants" % (ifn, counter))

    counter = 0
    ifd = open(ifn1, "r")
    for line in ifd:
        if re.match("^#", line):
            #print "Skip comment: %s" % (line)
            continue

        line = line.strip()
        items = line.split('\t')
        key = MakeKey(pos, "", items)

        if key in answer:
            answer[key] += 2
        else:
            answer[key] = 2
        #print(key + " = " + str(answer[key]))  # INFO
        counter += 1
        if counter % 1000000 == 0:
            print("%d %s" % (counter, key))
    ifd.close()
    print("%s has %d variants" % (ifn1, counter))

    counter = 0
    ifd = open(ifn2, "r")
    for line in ifd:
        if re.match("^#", line):
            #print "Skip comment: %s" % (line)
            continue

        line = line.strip()
        items = line.split('\t')
        key = MakeKey(pos, "", items)

        if key in answer:
            answer[key] += 4
        else:
            answer[key] = 4
        #print(key + " = " + str(answer[key]))  # INFO
        counter += 1
        if counter % 1000000 == 0:
            print("%d %s" % (counter, key))

    ifd.close()
    print("%s has %d variants" % (ifn2, counter))

    ## statistical result
    amount = [0, 0, 0, 0, 0, 0, 0, 0]
    for key in answer.keys():
        amount[0] += 1
        amount[answer[key]] += 1
    for i in range(0,8):
        print("%d\t%d" % (i, amount[i]))

    counter = 0
    ## data extraction
    ifd = open(ifn, "r")
    ofd = open(ifn+".1.vcf", "w")
    for line in ifd:
        if re.match("^#", line):
            #print "Skip comment: %s" % (line)
            ofd.write("%s" % line)
            continue

        line = line.strip()
        items = line.split('\t')
        items[3], items[4] = items[4], items[3]
        key = MakeKey(pos, "", items)

        if key in answer:
            ofd.write("%s\n" % line)
            counter += 1
    ofd.close()
    ifd.close()
    print("%s has %d variants" % (ifn+".1.vcf", counter))


    counter = 0
    ## data extraction
    ofd = open(ifn1+".3.vcf", "w")
    for key in answer.keys():
        if answer[key] == 3:
            ofd.write("%s\n" % line)
            counter += 1
    ofd.close()
    print("%s has %d variants" % (ifn1+".3.vcf", counter))

    counter = 0
    ## data extraction
    ofd = open(ifn2+".5.vcf", "w")
    for key in answer.keys():
        if answer[key] == 5:
            ofd.write("%s\n" % line)
            counter += 1
    ofd.close()
    print("%s has %d variants" % (ifn2+".5.vcf", counter))

    return


def main(argv):
    infile = ""
    infile1 = ""
    infile2 = ""
    check_position = 0

    try:
        opts, args = getopt.getopt(argv,"hi:j:k:p")
    except getopt.GetoptError:
        Usage()
        sys.exit(1)

    for opt, arg in opts:
        if opt == '-h':
            Usage()
            sys.exit()
        elif opt in ("-i"):
            infile = arg
        elif opt in ("-j"):
            infile1 = arg
        elif opt in ("-k"):
            infile2 = arg
        elif opt in ("-p"):
            check_position = 1
    if not infile or not infile1 or not infile2:
        Usage()
        sys.exit(2)
    elif (os.path.isfile(infile) == False):
        print("Error: input file(%s) is not existed" % (infile))
        Usage()
        sys.exit(3)
    elif (os.path.isfile(infile1) == False):
        print("Error: input file(%s) is not existed" % (infile1))
        Usage()
        sys.exit(4)
    elif (os.path.isfile(infile2) == False):
        print("Error: input file(%s) is not existed" % (infile2))
        Usage()
        sys.exit(5)

    #Main Function
    VCFComparision(infile, infile1, infile2, check_position)

    return

if __name__ == '__main__':
    main(sys.argv[1:])
