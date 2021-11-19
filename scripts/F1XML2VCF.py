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
import elementpath
import xml.etree.ElementTree as ET
from collections import defaultdict
import untangle
## refer to https://www.hellocodeclub.com/how-to-convert-xml-to-json-in-python-ultimate-guide/
import xmltodict
import json

# CONSTANT
VCF_HEADER = ""

def usage():
    print("F1XML2VCF.py -i <Input XML file >")
    print("Argument:")
    print("\t-h: Usage")
    print("\t-i: input file (.xml)")
    print("Usage:")
    print("\tpython ./F1XML2VCF.py -i ~/data/VGH-NGS/data/FoundationOne/S110-99395_\(PF21001\)_deidentified.xml ")

    return


def xml2(ifile):
    with open(ifile, 'r') as ifd:
        content = ifd.read().replace('<e2>', '')

    obj = xmltodict.parse(content)
    #print(json.dumps(obj))

    for variant in obj["rr:ResultsReport"]["rr:ResultsPayload"]["variant-report"]["short-variants"]["short-variant"]:
        #print(json.dumps(variant))
        print("%s %s" % (variant['@position'], variant['@cds-effect']))
    return


def xml1(ifile):

    with open(ifile, 'r') as ifd:
        content = ifd.read().replace('<e2>', '') #.replace('rr:', 'rr_')

    obj = untangle.parse(content)
    # SNP 2 VCF
    # print(obj.rr_ResultsReport.rr_ResultsPayload.variant_report.short_variants)
    for variant in obj.rr_ResultsReport.rr_ResultsPayload.variant_report.short_variants.iter():
        print("＝》%s,%s" % (variant, variant['position']))
    #TODO: CNV 2 VCF

    return


def xml(ifile):

    ## Note: due to invalid xml format from FoundationOne
    #tree = ET.parse(ifile)
    #root = tree.getroot()

    ## Note: xml parsing error due to <rr:*>
    ## error message:
    ##        xml.etree.ElementTree.ParseError: unbound prefix: line 2, column 0

    ## Note: can't find any tags whose name contains '-'
    ## e.g. /rr_ResultsReport/rr_ResultsPayload/variant-report/*

    with open(ifile, 'r') as ifd:
        content = ifd.read().replace('<e2>', '').replace('rr:', 'rr_')#.replace('B2-xyz', 'B2-xyz')
    print(content)
    root = ET.fromstring(content)
    # tree = ET.XMLParser.feed(content)
    # for elem in tree.iter():
    #     print (elem.tag, elem.attrib)

    # for child in root:
    #     print("%s,%s" % (child.tag, child.attrib))

    #selector = elementpath.Selector('/ResultsReport/ResultsPayload/variant-report/short-variants/*')
    #selector = elementpath.Selector('/rr_ResultsReport/rr_ResultsPayload/FinalReport/*')
    selector = elementpath.Selector('/rr_ResultsReport/rr_ResultsPayload/variant-report*/shrot-variants*/*')
    #selector = elementpath.Selector('/rr_A/variant-report*/*')

    print(selector.select(root))
    for item in selector.select(root):
        print("＝》%s,%s" % (item.tag, item.attrib))

    return


def main(argv):
    ifile = ""

    try:
        opts, args = getopt.getopt(argv, "hi:")
    except getopt.GetoptError:
        usage()
        sys.exit(1)

    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit()
        elif opt == '-i':
            ifile = arg

    if ifile == "":
        print("Error: '-i' is required")
        usage()
        sys.exit(2)
    elif not os.path.isfile(ifile):
        print("Error: input file(%s) is not existed" % ifile)
        usage()
        sys.exit(3)

    xml2(ifile)

    return


if __name__ == '__main__':
    main(sys.argv[1:])
