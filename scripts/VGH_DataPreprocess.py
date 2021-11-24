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

# CONST
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

TMP_FOLDER = "/tmp/unzip"
CMD_RMTMP = "rm -rf %s" % TMP_FOLDER
SRC = "/Users/chungtsai_su/src/github/AnomeINDEL/scripts"
CMD_F1CDX = "python " + SRC + "/F1XML2VCF.py -i %s -o %s"
CMD_UNZIP = "unzip %s -d " + TMP_FOLDER


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
    print("\tpython ./VGH_DataPreprocess.py -t ACTOnco+ -i ~/data/VGH-NGS/data/行動基因 -o ~/data/VGH-NGS/vcf/ACTG")
    print("\tpython ./VGH_DataPreprocess.py -t F1CDx -i ~/data/VGH-NGS/data/FoundationOne -o ~/data/VGH-NGS/vcf/FoundationOne")
    print("\tpython ./VGH_DataPreprocess.py -t OncomineBRCA -i ~/data/VGH-NGS/data/OncomineBRCA -o ~/data/VGH-NGS/vcf/OncomineBRCA")
    print("\tpython ./VGH_DataPreprocess.py -t OncomineFocus -i ~/data/VGH-NGS/data/OncomineFocus -o ~/data/VGH-NGS/vcf/OncomineFocus")
    print("\tpython ./VGH_DataPreprocess.py -t OncomineMyeloid -i ~/data/VGH-NGS/data/OncomineMyeloid -o ~/data/VGH-NGS/vcf/OncomineMyeloid")
    print("\tpython ./VGH_DataPreprocess.py -t OncomineTMB -i ~/data/VGH-NGS/data/OncomineTMB -o ~/data/VGH-NGS/vcf/OncomineTMB")

    return


def extract_samplename(file):
    fn = os.path.basename(file)
    s = fn.find("(")
    e = fn.find(")")

    if 0 < s < e:
        logging.debug("sample name :%s (%d,%d) from %s" % (fn[s + 1:e], s, e, file))
        return fn[s + 1:e]
    else:
        e = fn.find("_")
        if e > 0:
            logging.debug("Sample name :%s (0,%d) from %s" % (fn[:e], e, file))
            return fn[:e]

    logging.warning("Can't find sample name from %s" % fn)

    return ""


def actonco_preprocess(ifolder, ofolder):
    for file in glob.glob(ifolder + "/*.vcf"):
        logging.debug("%s" % file)
        sample_name = extract_samplename(file)
        logging.debug("Sample Name = %s" % sample_name)

        ofile = "%s/%s.vcf" % (ofolder, sample_name)
        with open(ofile, "w") as ofd:
            ifd = open(file, "r")
            for line in ifd:
                if line.startswith("#CHROM"):
                    items = line.strip().rsplit("\t", 1)
                    ofd.write("%s\t%s\n" % (items[0], sample_name))
                else:
                    ofd.write("%s" % line)
            ifd.close()

    return


def f1cdx_preprocess(ifolder, ofolder):
    for file in glob.glob(ifolder + "/*.xml"):
        if not os.path.isfile(file):
            logging.warning("WARNING: %s can't be found" % file)

        logging.debug("file = %s" % file)
        sample_name = extract_samplename(file)
        logging.debug("Sample Name = %s" % sample_name)
        ofile = "%s/%s.vcf" % (ofolder, sample_name)
        cmd = CMD_F1CDX % (file.replace('(', '\(').replace(')', '\)').replace(' ', '\\ '), ofile)
        # subprocess.run(cmd)
        try:
            res = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
            (output, err) = res.communicate()
            p_status = res.wait()
        except IOError:
            logging.ERROR("WARNING: %s " % res.stderr)

    return


def oncomine_preprocess(ifolder, ofolder):
    num_success = 0
    for file in glob.glob(ifolder + "/*.zip"):
        if not os.path.isfile(file):
            logging.warning("WARNING: %s can't be found" % file)
        try:
            res = subprocess.Popen(CMD_RMTMP, stdout=subprocess.PIPE, shell=True)
            (output, err) = res.communicate()
            p_status = res.wait()
        except IOError:
            logging.ERROR("WARNING: %s " % res.stderr)

        logging.debug("file = %s" % file)
        sample_name = extract_samplename(file)
        logging.debug("Sample Name = %s" % sample_name)
        #unzip
        cmd = CMD_UNZIP % file.replace(' ', '\\ ')
        try:
            logging.debug("%s" % cmd)
            res = subprocess.Popen(r"%s" % cmd, stdout=subprocess.PIPE, shell=True)
            (output, err) = res.communicate()
            p_status = res.wait()
        except IOError:
            logging.ERROR("WARNING: %s " % res.stderr)

        ofile = "%s/%s.vcf" % (ofolder, sample_name)

        flist = glob.glob(TMP_FOLDER + "/Variants/*/*Non-Filtered*.vcf")
        if len(flist) <= 0:
            ## Note: for BR21014
            flist = glob.glob(TMP_FOLDER + "/*/Variants/*/*Non-Filtered*.vcf")
        for vcf in flist:
            num_success += 1
            logging.debug("vcf = %s" % vcf)
            with open(ofile, "w") as ofd:
                ifd = open(vcf, "r")
                count = 0
                for line in ifd:
                    if line.startswith("#CHROM"):
                        items = line.strip().rsplit("\t", 1)
                        ofd.write("%s\t%s\n" % (items[0], sample_name))
                    else:
                        if not line.startswith("#"):
                            count += 1
                        ofd.write("%s" % line)
                ifd.close()
                logging.info("There are %d short variants in %s" % (count, sample_name))

    logging.info("There are %d samples successfully processed." % num_success)

    return


def main(argv):
    ifolder = ""
    ofolder = ""
    panel_type = ""

    logging.basicConfig(level=os.environ.get("LOGLEVEL", "INFO"))

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

    if not os.path.isdir(ofolder):
        os.mkdir(ofolder)

    if panel_type == STR_ACTOncoPlus:
        actonco_preprocess(ifolder, ofolder)
    elif panel_type == STR_F1CDx:
        f1cdx_preprocess(ifolder, ofolder)
    elif panel_type == STR_OncomineBRCA or panel_type == STR_OncomineFocus \
            or panel_type == STR_OncomineMyeloid or panel_type == STR_OncomineTMB:
        oncomine_preprocess(ifolder, ofolder)

    # match PanelType[panel_type]:
    #     case 1: #"ACTOnco+"
    #         actonco_preprocess(ifolder, ofolder)
    #     case 2: #"F1CDx"
    #         f1cdx_preprocess(ifolder, ofolder)
    #     case _:
    #         raise ValueError("Not a valid")

    return


if __name__ == '__main__':
    main(sys.argv[1:])
