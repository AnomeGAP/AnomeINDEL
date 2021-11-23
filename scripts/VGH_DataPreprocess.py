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
# @file    VGH_DataPreprocess.py
#
# @brief   Parse VCF/XML/ZIP of data folder provided by VGH
#
# @author  Chung-Tsai Su(chungtsai.su@atgenomix.com)
#
# @date    2021/11/23
#
# @version 1.0
#
# @remark
#
import subprocess
import sys
import getopt
import os
import logging
import glob
from enum import Enum

#CONST
STR_ACTOncoPlus = "ACTOnco+"
STR_F1CDx = "F1CDx"
STR_ArcherFusionPlex = "ArcherFusionPlex"
STR_Guardant360 = "Guardant360"
STR_OncomineBRCA = "OncomineBRCA"
STR_OncomineFocus = "OncomineFocus"
STR_OncomineMyeloid = "OncomineMyeloid"
STR_OncomineTMB = "OncomineTMB"

PanelType = {STR_ACTOncoPlus: 1,
             STR_F1CDx: 2,
             STR_ArcherFusionPlex: 3,
             STR_Guardant360: 4,
             STR_OncomineBRCA: 5,
             STR_OncomineFocus: 6,
             STR_OncomineMyeloid: 7,
             STR_OncomineTMB: 8
}
SRC = "/Users/chungtsai_su/src/github/AnomeINDEL/scripts"
CMD_F1CDX = "python " + SRC + "/F1XML2VCF.py -i %s -o %s"


def usage():
    print("VGH_DataPreprocess.py -t <Data type> -i <Input data folder> -o <Output data folder>")
    print("Argument:")
    print("\t-h: Usage")
    s = ""
    for k in PanelType.keys():
        if s == "":
            s = k
        else:
            s += "|" + k
    print("\t-t: data type [%s]" % s)
    print("\t-i: input folder")
    print("\t-o: output folder")
    print("Usage:")
    print("\tpython ./VGH_DataPreprocess.py -t F1CDx -i ~/data/VGH-NGS/data/FoundationOne "
          "-o ~/data/VGH-NGS/vcf/FoundationOne ")

    return

def extract_samplename(file):
    s = file.find("(")
    e = file.find(")")
    logging.debug("sample name :%s (%d,%d)" % (file[s + 1:e], s, e))
    if s < 0 or e < 0 or s > e:
        logging.warning("Can't find sample name from %s" % file)
        sys.exit(5)

    return file[s + 1:e]


def actonco_preprocess(ifolder, ofolder):
    for file in glob.glob(ifolder + "/*.xml"):
        logging.debug("%s" % file)
        sample_name = extract_samplename(file)
        logging.debug("Sample Name = %s" % sample_name)

    return


def f1cdx_preprocess(ifolder, ofolder):
    if not os.path.isdir(ofolder):
        os.mkdir(ofolder)

    for file in glob.glob(ifolder + "/*.xml"):
        if not os.path.isfile(file):
            logging.warning("WARNING: %s can't be found" % file)

        logging.debug("file = %s" % file)
        sample_name = extract_samplename(file)
        logging.debug("Sample Name = %s" % sample_name)
        ofile = "%s/%s.vcf" % (ofolder, sample_name)
        cmd = CMD_F1CDX % (file.replace('(', '\(').replace(')', '\)'), ofile)
        # subprocess.run(cmd)
        try:
            res = subprocess.Popen(cmd, shell=True)
        except IOError:
            logging.ERROR("WARNING: %s " % res.stderr)

    return


def main(argv):
    ifolder = ""
    ofolder = ""
    panel_type = ""

    logging.basicConfig(level=os.environ.get("LOGLEVEL", "DEBUG"))

    try:
        opts, args = getopt.getopt(argv, "ht:i:o:")
    except getopt.GetoptError:
        usage()
        sys.exit(1)

    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit()
        elif opt == '-i':
            ifolder = arg
            if ofolder == "":
                ofolder = "%s.vcf" % ifolder
        elif opt == '-o':
            ofolder = arg
        elif opt == '-t':
            panel_type = arg

    if ifolder == "":
        print("Error: '-i' is required")
        usage()
        sys.exit(2)
    elif not os.path.isdir(ifolder):
        print("Error: input folder(%s) is not existed" % ifolder)
        usage()
        sys.exit(3)

    if panel_type not in PanelType:
        print("Error: %s is not a valid panel name" % panel_type)
        usage()
        sys.exit(4)

    print("%s %d" % (panel_type, PanelType[panel_type]))
    if panel_type == STR_ACTOncoPlus:
        actonco_preprocess(ifolder, ofolder)
    elif panel_type == STR_F1CDx:
        f1cdx_preprocess(ifolder, ofolder)

    # match PanelType[panel_type]:
    #     case 1: #"ACTOnco+"
    #         actonco_preprocess(ifolder, ofolder)
    #     case 2: #"F1CDx"
    #         f1cdx_preprocess(ifolder, ofolder)
    #     case 3: #"ArcherFusionPlex"
    #
    #     case 4: #"Guardant360"
    #
    #     case 5|6|7: #"OncomineBRCA"|"OncomineFocus"|"OncomineMyeloid"
    #
    #     case 8: #"OncomineTMB"
    #
    #     case _:
    #         raise ValueError("Not a valid")

    return


if __name__ == '__main__':
    main(sys.argv[1:])
