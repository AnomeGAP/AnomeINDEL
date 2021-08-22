#!/bin/env python
#
# @note Copyright (C) 2020, Atgenomix Incorporated. All Rights Reserved.
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
# @file    PrecisionFDA-TMB-VCF.py
#
# @brief   VCF filtering for PrecisionFDA TMB challenge (Phase II)
#
# @author  Chung-Tsai Su(chungtsai.su@atgenomix.com)
#
# @date    2021/07/29
#
# @version 1.0
#
# @remark
#

import sys
import getopt
import re
import vcf
from collections import defaultdict
import logging
import os


def usage():
    print(
        "PrecisionFDA-TMB.py -i <Input VCF file> -a <VEP Annotation file>")
    print("Argument:")
    print("\t-h: Usage")
    print("\t-i: Input Merged VCF file")
    print("\t-a: VEP Annotation file")
    print("Usage:")
    print("\tpython3 ./PrecisionFDA-TMB-VCF.py -i ~/data/PrecisionFDA/TMB/TMB2/A.vcf -a ~/data/PrecisionFDA/TMB/TMB2-VEP/A")

    return


class Variant(object):
    def __init__(self) -> object:
        self.id = ""
        self.is_frameshift = 0  #False
        self.is_inframe = 0     #False
        self.is_missense = 0    #False
        self.is_synonymous = 0  #False
        self.impact_level = 0   # 0:[LOW|MODIFIER]; 1:MODERATE; 2:HIGH
        self.num_transcripts = 0

        self.is_snp = 0
        self.is_insertion = 0
        self.is_deletion = 0
        self.is_sv = 0
        self.is_qualified = 0


def parse_vep(afn):
    h_var = defaultdict(Variant)
    pre_k = ""
    ifd = open(afn, "r")
    for line in ifd:
        if line.startswith("#"):
            continue

        # logging.debug("%s" % line.strip())
        (k, v) = line.strip().split("\t", 1)
        if k != pre_k:
            h_var[k] = Variant()
        # logging.debug("k=%s" % k)
        # logging.debug("v=%s" % v)

        if v.find("frameshift") > 0:
            setattr(h_var[k], 'is_frameshift', 1)
        elif v.find("inframe") > 0:
            setattr(h_var[k], 'is_inframe', 1)
        elif v.find("missense") > 0:
            setattr(h_var[k], 'is_missense', 1)
        elif v.find("synonymous") > 0:
            setattr(h_var[k], 'is_synonymous', 1)

        impact = getattr(h_var[k], 'impact_level')
        if v.find("IMPACT=HIGH") > 0:
            setattr(h_var[k], 'impact_level', 2)
        elif v.find("IMPACT=MODERATE") > 0 and impact < 2:
            setattr(h_var[k], 'impact_level', 1)

        setattr(h_var[k], 'num_transcripts', getattr(h_var[k], 'num_transcripts') + 1)
        pre_k = k

    num = 0
    for k in h_var:
        num += 1
        logging.debug("%s\t%d:%d:%d:%d:%d:%d" % (k, getattr(h_var[k], 'is_frameshift'), getattr(h_var[k], 'is_inframe'),
                                                 getattr(h_var[k], 'is_missense'), getattr(h_var[k], 'is_synonymous'),
                                                 getattr(h_var[k], 'impact_level'), getattr(h_var[k], 'num_transcripts')))
    logging.debug("%d variants on VEP" % num)
    ifd.close()
    return h_var


def process(ifn, afn):
    ifd = open(ifn, "r")
    vcf_reader = vcf.Reader(ifd)
    h_var = parse_vep(afn)
    num_total = 0
    num_qualified = 0

    for record in vcf_reader:
        # *``Record.CHROM``
        # *``Record.POS``
        # *``Record.ID``
        # *``Record.REF``
        # *``Record.ALT``
        # *``Record.QUAL``
        # *``Record.FILTER``
        # *``Record.INFO``
        # *``Record.FORMAT``
        # *``Record.samples``
        # *``Record.genotype``

        num_total += 1

        logging.debug("REF=%s; ALT=%s" % (record.REF, record.ALT))
        logging.debug("INFO=%s" % record.INFO)
        vtype = 0   #0:SNP; 1:INSERTION; 2:DELETION; 3:SV
        if record.INFO.get("SVTYPE") is not None:
            vtype = 3
            logging.debug("SV")
        elif len(record.REF) < len(record.ALT[0]):
            vtype = 1
            logging.debug("INSERTION %s=>%s" % (record.REF, record.ALT[0]))
        elif len(record.REF) > len(record.ALT[0]):
            vtype = 2
            logging.debug("DELETION")

        #Variant Normalization for VEP
        offset = 1 if vtype != 3 and len(record.REF) != len(record.ALT[0]) else 0
        ref = record.REF
        alt = str(record.ALT[0])
        logging.debug("type(ref) = %s, type(alt) = %s" % (type(ref), type(alt)))

        if vtype == 1:
            alt = alt[len(ref):]
            ref = "-"
        elif vtype == 2:
            ref = ref[len(alt):]
            alt = "-"

        k = "%s_%s_%s/%s" % (record.CHROM, record.POS+offset, ref, alt)
        logging.debug("%s" % k)
        setattr(h_var[k], 'id', "%s_%s" % (record.CHROM, record.POS))

        if vtype == 0:
            setattr(h_var[k], 'is_snp', 1)
        elif vtype == 1:
            setattr(h_var[k], 'is_insertion', 1)
        elif vtype == 2:
            setattr(h_var[k], 'is_deletion', 1)
        elif vtype == 3:
            setattr(h_var[k], 'is_sv', 1)

        # Qualification

        if len(record.FILTER) > 0:
            logging.debug("removed by FILTER")
            continue
        if record.POS == 17296837:
            logging.info("pass FILTER\t%s" % record.FILTER)

        # logging.debug("%s" % record.samples[0])
        if record.FORMAT.find("AF") < 0:
            logging.debug("no AF record")
            continue

        if record.samples[0]['AF'] < 0.05:
            logging.debug("removed by AF")
            continue

        if record.FORMAT.find("AD") < 0:
            logging.debug("no AD record")
            continue

        # logging.debug("AD=%s" % record.samples[0]['AD'])
        if record.samples[0]['AD'][1] < 3:
            logging.debug("AD < 3")
            continue

        if record.samples[0]['AD'][0] + record.samples[0]['AD'][1] < 325:
            logging.debug("t_depth < 25")
            continue

        if k not in h_var:
            logging.debug("%s not in VEP" % k)
            continue

        if h_var[k].is_frameshift + h_var[k].is_inframe + h_var[k].is_missense == 0:
            logging.debug("%s : %s" % (k, h_var[k]))
            continue

        setattr(h_var[k], 'is_qualified', 1)
        logging.debug("PASS: %s" % k)
        # logging.debug("%s\t%s\t%s\t%s\t%s" % (k, record.samples[0]['AD'], record.genotype, record.FORMAT, record.FILTER))
        num_qualified += 1
        # break
    logging.debug("%d/%d" % (num_qualified, num_total))
    ifd.close()

    h_qualified = defaultdict(int)

    num_qualified = 0
    (all_highimpact, all_moderateimpact, all_lowimpact, all_snp, all_insertion, all_deletion, all_indel, all_sv,
     all_frameshift, all_inframe, all_missense, all_synonymous, num_highimpact, num_moderateimpact, num_lowimpact,
     num_snp, num_insertion, num_deletion, num_indel, num_sv, num_frameshift, num_inframe, num_missense, num_synonymous) \
        = (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    for k in h_var:
        if h_qualified[h_var[k].id] == 1:
            continue

        if h_var[k].impact_level == 2:
            all_highimpact += 1
        elif h_var[k].impact_level == 1:
            all_moderateimpact += 1
        else:
            all_lowimpact += 1

        if h_var[k].is_snp:
            all_snp += 1
        elif h_var[k].is_insertion:
            all_insertion += 1
            all_indel += 1
        elif h_var[k].is_deletion:
            all_deletion += 1
            all_indel += 1
        else:
            all_sv += 1

        if h_var[k].is_frameshift:
            all_frameshift += 1
        elif h_var[k].is_inframe:
            all_inframe += 1
        elif h_var[k].is_missense:
            all_missense += 1
        elif h_var[k].is_synonymous:
            all_synonymous += 1

        if h_var[k].is_qualified == 1:
            num_qualified += 1
            if h_var[k].impact_level == 2:
                num_highimpact += 1
            elif h_var[k].impact_level == 1:
                num_moderateimpact += 1
            else:
                num_lowimpact += 1

            if h_var[k].is_snp:
                num_snp += 1
            elif h_var[k].is_insertion:
                num_insertion += 1
                num_indel += 1
            elif h_var[k].is_deletion:
                num_deletion += 1
                num_indel += 1
            else:
                num_sv += 1

            if h_var[k].is_frameshift:
                num_frameshift += 1
            elif h_var[k].is_inframe:
                num_inframe += 1
            elif h_var[k].is_missense:
                num_missense += 1
            elif h_var[k].is_synonymous:
                num_synonymous += 1
            h_qualified[h_var[k].id] = 1
            logging.debug("%s\t%d:%d:%d:%d:%d:%d:%d:%d:%d:%d" % (k, h_var[k].is_snp, h_var[k].is_insertion,
                                                                 h_var[k].is_deletion, h_var[k].is_sv,
                                                                 h_var[k].is_frameshift, h_var[k].is_inframe,
                                                                 h_var[k].is_missense, h_var[k].is_synonymous,
                                                                 h_var[k].num_transcripts, h_var[k].impact_level))
    logging.debug("%d\t%d" % (num_qualified, num_total))
    print(" 1:%d 2:%d 3:%d 4:%d 5:%d 6:%d 7:%d 8:%d 9:%d 10:%d 11:%d 12:%d 13:%d"
          " 14:%d 15:%d 16:%d 17:%d 18:%d 19:%d 20:%d 21:%d 22:%d 23:%d 24:%d 25:%d 26:%d" %
          (num_total, num_qualified, all_highimpact, all_moderateimpact, all_lowimpact, all_snp, all_insertion,
           all_deletion, all_indel, all_sv, all_frameshift, all_inframe, all_missense, all_synonymous,
           num_highimpact, num_moderateimpact, num_lowimpact, num_snp, num_insertion, num_deletion, num_indel,
           num_sv, num_frameshift, num_inframe, num_missense, num_synonymous))


def main(argv):
    ifile = ""
    afile = ""

    logging.basicConfig(level=os.environ.get("LOGLEVEL", "INFO"))
    try:
        opts, args = getopt.getopt(argv, "hi:a:")
    except getopt.GetoptError:
        usage()
        sys.exit(1)

    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit()
        elif opt == "-i":
            ifile = arg
        elif opt == "-a":
            afile = arg

    process(ifile, afile)
    return

if __name__ == '__main__':
    main(sys.argv[1:])
