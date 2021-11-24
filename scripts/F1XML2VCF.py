#!/bin/env python
#
# @note Copyright (C) 2021, Atgenomix Incorporated. All Rights Reserved.
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
# @file    F1XML2VCF.py
#
# @brief   Parse Foundation One CDx xml report and transform to CSV format
#
# @author  Chung-Tsai Su(chungtsai.su@atgenomix.com)
#
# @date    2021/11/09
#
# @version 1.0
#
# @remark
#

import sys
import getopt
import os
import re
from collections import defaultdict
import json
import logging

## refer to https://www.hellocodeclub.com/how-to-convert-xml-to-json-in-python-ultimate-guide/
import xmltodict

# CONSTANT
VCF_HEADER = '''##fileformat=VCFv4.1
##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele frequency based on Flow Evaluator observation counts">
##FORMAT=<ID=AO,Number=A,Type=Integer,Description="Alternate allele observation count">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=LEN,Number=A,Type=Integer,Description="allele length">
##INFO=<ID=END,Number=A,Type=Integer,Description="the end position for the sv">
##INFO=<ID=TYPE,Number=A,Type=String,Description="The type of allele, either snp, mnp, ins, del, dup, or complex.">
##INFO=<ID=HGVS,Number=A,Type=String,Description="the protein effect of the allele">
##INFO=<ID=FUNCTIONAL,Number=A,Type=String,Description="the functional effect of the allele">
##INFO=<ID=TRANSCRIPT,Number=A,Type=String,Description="the transcript of the allele">
##INFO=<ID=GENE,Number=1,Type=String,Description="the gene name of the SV">
##INFO=<ID=EXON,Number=1,Type=String,Description="the exon position of the SV">
##INFO=<ID=COPYNUMBER,Number=1,Type=Integer,Description="the copy number of the SV">
##INFO=<ID=RATIO,Number=1,Type=Float,Description="the copy ratio of the SV">
##INFO=<ID=EQUIVOCAL,Number=1,Type=Bool,Description="the equivocal of the SV">
##INFO=<ID=STATUS,Number=1,Type=String,Description="the status of the SV">
##contig=<ID=chr1,length=249250621,assembly=hg19>
##contig=<ID=chr2,length=243199373,assembly=hg19>
##contig=<ID=chr3,length=198022430,assembly=hg19>
##contig=<ID=chr4,length=191154276,assembly=hg19>
##contig=<ID=chr5,length=180915260,assembly=hg19>
##contig=<ID=chr6,length=171115067,assembly=hg19>
##contig=<ID=chr7,length=159138663,assembly=hg19>
##contig=<ID=chr8,length=146364022,assembly=hg19>
##contig=<ID=chr9,length=141213431,assembly=hg19>
##contig=<ID=chr10,length=135534747,assembly=hg19>
##contig=<ID=chr11,length=135006516,assembly=hg19>
##contig=<ID=chr12,length=133851895,assembly=hg19>
##contig=<ID=chr13,length=115169878,assembly=hg19>
##contig=<ID=chr14,length=107349540,assembly=hg19>
##contig=<ID=chr15,length=102531392,assembly=hg19>
##contig=<ID=chr16,length=90354753,assembly=hg19>
##contig=<ID=chr17,length=81195210,assembly=hg19>
##contig=<ID=chr18,length=78077248,assembly=hg19>
##contig=<ID=chr19,length=59128983,assembly=hg19>
##contig=<ID=chr20,length=63025520,assembly=hg19>
##contig=<ID=chr21,length=48129895,assembly=hg19>
##contig=<ID=chr22,length=51304566,assembly=hg19>
##contig=<ID=chrX,length=155270560,assembly=hg19>
##contig=<ID=chrY,length=59373566,assembly=hg19>
##contig=<ID=chrM,length=16569,assembly=hg19>
##reference=hg19
##source="F1CDx to VCF by ATGENOMIX v1.0"
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	%s
'''
QUAL_DEFAULT = "."
FILTER_DEFAULT = "PASS"
TYPE_SNP = "snp"
TYPE_INS = "ins"
TYPE_DUP = "dup"
TYPE_DEL = "del"


class SNP:
    """
    SNP from F1CDx
    """
    def __init__(self, object):
        (chrom, pos) = object['@position'].split(":")
        s = object['@cds-effect'].find(">")
        if s > 0:
            self.type = TYPE_SNP
            ref = object['@cds-effect'][s-1]
            alt = object['@cds-effect'][s+1]
            self.id = "%s_%s_%s>%s" % (chrom, pos, ref, alt)
        elif object['@cds-effect'].find("ins") > 0:
            self.type = TYPE_INS
            ref = ""
            alt = object['@cds-effect'][object['@cds-effect'].find("ins")+3:]
            self.id = "%s_%s_%dins%s" % (chrom, pos, int(pos)+1, alt)
        elif object['@cds-effect'].find('del'):
            self.type = TYPE_DEL
            ref = object['@cds-effect'][object['@cds-effect'].find("del")+3:]
            alt = ""
            self.id = "%s_%s_%ddel" % (chrom, pos, int(pos)+len(ref))
        else:
            logging.warning("Can't find '>' in %s [%s]" % (object['@position'], object['@cds-effect']))
            return

        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.qual = QUAL_DEFAULT
        self.filter = FILTER_DEFAULT
        self.len = 1
        self.info = "LEN=%d;TYPE=%s;HGVS=%s;FUNCTIONAL=%s;GENE=%s;EQUIVOCAL=%s;STATUS=%s" % \
                    (self.len, self.type, object['@protein-effect'], object['@functional-effect'],
                     object['@gene'], object['@equivocal'], object['@status'])
        self.gt = "0/1"
        self.af = object['@allele-fraction']
        self.dp = object['@depth']
        self.ao = int(float(self.dp)*float(self.af))
        self.format = "GT:AF:AO:DP\t%s:%s:%s:%s" % (self.gt, self.af, self.ao, self.dp)

    def dump_vcf(self):
        result = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (self.chrom, self.pos, self.id, self.ref, self.alt, self.qual,
                                                         self.filter, self.info, self.format)
        return result


class CNV:
    """
    CNV from F1CDx
    """

    def __init__(self, object):
        (self.chrom, pos) = object['@position'].split(":")
        (start, end) = pos.split("-")
        self.alt = "DUP"
        self.type = TYPE_DUP
        if object['@type'] == 'loss':
            self.alt = "DEL"
            self.type = TYPE_DEL

        self.pos = start
        self.ref = "-"
        self.id = "%s_%s_%s%s" % (self.chrom, self.pos, end, self.alt)
        self.qual = QUAL_DEFAULT
        self.filter = FILTER_DEFAULT
        self.len = int(end) - int(start) + 1
        self.info = "LEN=%d;TYPE=%s;END=%s;GENE=%s;EXON=%s;COPYNUMBER=%s;RATIO=%s;EQUIVOCAL=%s;STATUS=%s" \
                    % (self.len, self.type, end, object['@gene'], object['@number-of-exons'], object['@copy-number'],
                       object['@ratio'], object['@equivocal'], object['@status'])
        self.gt = "0/1"
        self.format = "GT\t%s" % self.gt

    def dump_vcf(self):
        result = "%s\t%s\t%s\t%s\t<%s>\t%s\t%s\t%s\t%s" % (self.chrom, self.pos, self.id, self.ref, self.alt, self.qual,
                                                         self.filter, self.info, self.format)
        return result


def usage():
    print("F1XML2VCF.py -i <Input XML file> -o <Output VCF file>")
    print("Argument:")
    print("\t-h: Usage")
    print("\t-i: input file (.xml)")
    print("\t-o: output file (.vcf)")
    print("Usage:")
    print("\tpython ./F1XML2VCF.py -i ~/data/VGH-NGS/data/FoundationOne/S110-99395_\(PF21001\)_deidentified.xml "
          "-o ~/data/VGH-NGS/data/FoundationOne/S110-99395_\(PF21001\)_deidentified.vcf ")

    return


def xml2vcf(ifile, ofile):
    s = ifile.find("(")
    e = ifile.find(")")
    logging.debug("sample name :%s (%d,%d)" % (ifile[s+1:e], s, e))
    if s < 0 or e < 0 or s > e:
        logging.warning("Can't find sample name from %s" % ifile)
        sys.exit(4)

    sample_name = ifile[s+1:e]

    with open(ifile, 'r') as ifd:
        content = re.sub(r"<[\w|\d]{2}>", "", ifd.read())
        # content = ifd.read().replace('<e2>', '').replace('<c3>', '').replace('<85>', '')\
        #     .replace('<c4>', '').replace('<97>', '').replace('<93>', '').replace('<9d>', '')\
        #     .replace('<c5>', '').replace('<84>', '').replace('<c2>', '')

    obj = xmltodict.parse(content)

    with open(ofile, "w") as ofd:
        ofd.write(VCF_HEADER % sample_name)
        num_snv = 0
        #SNV
        for variant in obj["rr:ResultsReport"]["rr:ResultsPayload"]["variant-report"]["short-variants"]["short-variant"]:
            logging.debug("\t%s %s" % (variant['@position'], variant['@cds-effect']))
            var = SNP(variant)
            ofd.write("%s\n" % (var.dump_vcf()))
            num_snv += 1
        logging.info("There are %d short variants in %s" % (num_snv, sample_name))
        # #CNV
        # for variant in obj["rr:ResultsReport"]["rr:ResultsPayload"]["variant-report"]["copy-number-alterations"]["copy-number-alteration"]:
        #     logging.debug("\t%s %s" % (variant['@position'], variant['@copy-number']))
        #     var = CNV(variant)
        #     ofd.write("%s\n" % (var.dump_vcf()))

    return


def main(argv):
    ifile = ""
    ofile = ""

    logging.basicConfig(level=os.environ.get("LOGLEVEL", "INFO"))

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
            ifile = arg
            if ofile == "":
                ofile = "%s.vcf" % ifile
        elif opt == '-o':
            ofile = arg

    if ifile == "":
        print("Error: '-i' is required")
        usage()
        sys.exit(2)
    elif not os.path.isfile(ifile):
        print("Error: input file(%s) is not existed" % ifile)
        usage()
        sys.exit(3)

    xml2vcf(ifile, ofile)

    return


if __name__ == '__main__':
    main(sys.argv[1:])
