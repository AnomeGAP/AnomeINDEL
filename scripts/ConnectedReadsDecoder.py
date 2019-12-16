#!/bin/env python
#
# @note Copyright (C) 2019, Atgenomix Incorporated. All Rights Reserved.
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
# @file    ConnectedReadsDecoder.py
#
# @brief   Decode the contigs of ConnectedReads whose quality string are encoded by QD encoding algorithm
#
# @author  Chung-Tsai Su(chungtsai.su@atgenomix.com)
#
# @date    2019/06/27
#
# @version 1.0
#
# @remark
#

import sys
import getopt
import os.path
import re
import gzip
import io

# CONST
TOTAL_PHRED_SCORES = 42
DECODER_FULLY = 1
DECODER_ONLY_DEPTH = 2
DECODER_RESET_QUALITY = 3
DECODER_QUALITY_CHECK = 4

PHRAD33_LAYERS = 42
PHRAD33_BASE = '!'
PHRAD33_MAX = chr(ord(PHRAD33_BASE) + PHRAD33_LAYERS - 1)
DEPTH_NORMALIZATION_LEN = 51    #must be odd

# Default Parameter
DEFAULT_ENCODER = "QD-3-14"
DEFAULT_READ_LEN = 150
DEFAULT_QUALITY = 'F'
ENCODERS = {'QD-1-42': 1,
            'QD-2-21': 2,
            'QD-3-14': 3,
            'QD-6-7': 6,
            'QD-7-6': 7,
            'QD-14-3': 14,
            'QD-21-2': 21}
DECODERS = {'QD-DECODER': DECODER_FULLY,
            'DEPTH-ONLY': DECODER_ONLY_DEPTH,
            'QUALITY-SETTING': DECODER_RESET_QUALITY,
            'QUALITY-CHECK': DECODER_QUALITY_CHECK}


def usage():
    print("ConnectedReadsDecoder.py -i <Input FASTQ.gz> -e <ENCODER Type> -d <Decoding Function> -l <Read Length> "
          "-o <Output FASTQ.gz>")
    print("Argument:")
    print("\t-h: Usage")
    print("\t-i: Input FASTQ.gz file")
    print("\t-e: the encoder type (Default: %s)" % (ENCODERS[DEFAULT_ENCODER]))
    for k, v in ENCODERS.items():
        print("\t\t%s: %d levels for quality and %d levels for depth" % (k, ENCODERS[k], 42/ENCODERS[k]))
    print("\t-d: <Decoding Function> ")
    print("\t\t%s: both quality and depth are fully decoded" % 'QD-DECODER')
    print("\t\t%s: only depth are decoded" % 'DEPTH-ONLY')
    print("\t\t%s: quality string will be reset by '-q' parameter" % 'QUALITY-SETTING')
    print("\t\t%s: check quality string" % 'QUALITY-CHECK')
    print("\t-l: read length you expected when '-d QD-DECODER' [Default: %d]" % DEFAULT_READ_LEN)
    print("\t-q: read length you expected when '-d QUALITY-SETTING' [Default: 'F']")
    print("\t-o: Output FASTQ.gz file")
    print("Usage:")
    print("\tpython ./ConnectedReadsDecoder.py -i ../data/NA12878/result-1.0.2-qual-fix-5.fq.gz  "
          "-e %s -d %s -l %d -o ../data/NA12878/result-1.0.2-qual-fix-5.decoded.fq.gz "
          "" % ('QD-3-14', 'QD-DECODER', DEFAULT_READ_LEN))
    print("\tpython ./ConnectedReadsDecoder.py -i ../data/NA12878/result-1.0.2-qual-fix-5.fq.gz  "
          "-e %s -d %s -o ../data/NA12878/result-1.0.2-qual-fix-5.decoded.fq.gz "
          "" % ('QD-3-14', 'DEPTH-ONLY'))
    print("\tpython ./ConnectedReadsDecoder.py -i ../data/NA12878/result-1.0.2-qual-fix-5.fq.gz  "
          "-e %s -d %s -q %c -o ../data/NA12878/result-1.0.2-qual-fix-5.decoded.fq.gz "
          "" % ('QD-3-14', 'QUALITY-SETTING', DEFAULT_QUALITY))
    print("\tpython ./ConnectedReadsDecoder.py -i ../data/NA12878/result-1.0.2-qual-fix-5.fq.gz  "
          "-e %s -d %s -o ../data/NA12878/result-1.0.2-qual-fix-6.fq.gz "
          "" % ('QD-3-14', 'QUALITY-CHECK'))

    return


def qd_decoder(num_layers_quality, qual):
    (decoded_qual, decoded_depth) = ([], [])
    num_layer_depth = PHRAD33_LAYERS / num_layers_quality
    for k in str(qual):
        v = ord(k) - ord(PHRAD33_BASE)
        if v >= PHRAD33_LAYERS:
            v = PHRAD33_LAYERS - 1
        layer_quality = int(v / num_layer_depth)
        layer_depth = v - layer_quality * num_layer_depth + 1
        decoded_qual.append(layer_quality)
        decoded_depth.append(layer_depth)
        # print("%c => %d\t%d" % (k, layer_quality, layer_depth))

    return decoded_qual, decoded_depth


def depth_normalization(decoded_depth):
    normalized = []
    s = 0
    depth_list = list(decoded_depth)
    for i in range(len(depth_list)):
        s += depth_list[i]
    v = int(s / len(depth_list))
    for i in range(len(depth_list)):
        normalized.append(v)
    return normalized


def quality_recovery(decoded_qual, num_layers_depth):
    recovered = ""

    for i in decoded_qual:
        q = chr(int((i + 1) * num_layers_depth - 1 + ord(PHRAD33_BASE)))
        recovered += q
    return recovered


def fully_decoder(num_layers_quality, decoder, read_length, quality, name, sequ, qual):
    result = ""
    if decoder == DECODER_FULLY:
        (decoded_qual, decoded_depth) = qd_decoder(num_layers_quality, qual)

    elif decoder == DECODER_ONLY_DEPTH:
        num_layers_depth = PHRAD33_LAYERS / num_layers_quality
        (decoded_qual, decoded_depth) = qd_decoder(num_layers_quality, qual)
        normalized_depth = depth_normalization(decoded_depth)
        recovered_quality = quality_recovery(decoded_qual, num_layers_depth)
        for i in range(normalized_depth[0]):
            result += "%s-%d\n%s\n+\n%s\n" % (name, i, sequ, recovered_quality)
    return result


def decoder_fun(ifile, ofile, num_layers_quality, decoder, read_length, quality):
    i = 0
    name = ""
    sequ = ""
    qual = ""

    ofd = gzip.open(ofile, 'wb')
    with gzip.open(ifile, 'rb') as gz:
        f = io.BufferedReader(gz)
        for line in f:
            if i % 4 == 3:
                if decoder == DECODER_FULLY or decoder == DECODER_ONLY_DEPTH:
                    # TODO: need to implement the logic
                    qual = line.decode().strip()
                    contig = fully_decoder(num_layers_quality, decoder, read_length, quality, name, sequ, qual)
                elif decoder == DECODER_RESET_QUALITY:
                    qual = DEFAULT_QUALITY * len(sequ)
                    contig = "%s\n%s\n+\n%s\n" % (name, sequ, qual)
                elif decoder == DECODER_QUALITY_CHECK:
                    # idx = 0
                    # print(line.decode().strip())
                    l_qual = []
                    for q in line.decode().strip():
                        v = ord(q) - ord(PHRAD33_BASE)
                        if v >= PHRAD33_LAYERS:
                            # print("q=%s" % q)
                            l_qual.append(PHRAD33_MAX)
                        else:
                            l_qual.append(q)
                        # if v < 28:
                        #     print("%d\t%c" % (idx,q))
                        # idx += 1
                    qual = "".join(l_qual)
                    contig = "%s\n%s\n+\n%s\n" % (name, sequ, qual)
                # print(contig)
                ofd.write(contig.encode())
                name = ""
                sequ = ""
                # break
            elif i % 4 == 0:
                name = line.decode().strip()
            elif i % 4 == 1:
                sequ = line.decode().strip()
            elif line.decode().strip() != "+":
                print("Format ERROR: '%s' should be '+'" % line.decode().strip())
                break
            i += 1

    ofd.close()
    return


def main(argv):
    ifile = ""
    ofile = ""
    num_layers_quality = ENCODERS['QD-3-14']
    decoder = DECODERS['QD-DECODER']
    read_length = DEFAULT_READ_LEN
    quality = DEFAULT_QUALITY

    try:
        opts, args = getopt.getopt(argv, "hi:e:d:l:q:o:")
    except getopt.GetoptError:
        usage()
        sys.exit(1)

    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit()
        elif opt == "-i":
            ifile = arg
        elif opt == "-e":
            if arg not in ENCODERS:
                print("Error: '-e %s' is not valid" % arg)
                usage()
                sys.exit(4)
            num_layers_quality = ENCODERS[arg]
        elif opt == "-d":
            if arg not in DECODERS:
                print("Error: '-d %s' is not valid" % arg)
                usage()
                sys.exit(5)
            decoder = DECODERS[arg]
        elif opt == "-l":
            read_length = int(arg)
        elif opt == "-q":
            quality = arg
        elif opt == "-o":
            ofile = arg

    # error handling for input parameters
    if ifile == "" or ofile == "":
        print("Error: '-i' and '-o' are required")
        usage()
        sys.exit(2)
    elif not os.path.isfile(ifile):
        print("Error: input file(%s) is not existed" % ifile)
        usage()
        sys.exit(3)

    if ofile == "":
        ofile = "%s.unmapped.sam" % ifile

    if decoder == DECODER_FULLY:
        print("[TBD] Coming soon. Please Stay Tune ~~~")
        sys.exit(0)

    # Main Function
    decoder_fun(ifile, ofile, num_layers_quality, decoder, read_length, quality)

    return


if __name__ == '__main__':
    main(sys.argv[1:])
