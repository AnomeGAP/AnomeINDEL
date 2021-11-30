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
# @file    VGH_ReportPreprocess.py
#
# @brief   Parse report(.txt) of report folder provided by VGH
#
# @author  Chung-Tsai Su(chungtsai.su@atgenomix.com)
#
# @date    2021/11/24
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
import re
from collections import defaultdict

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

STR_NA_TMB = '-1'
STR_NA = "NA"
STR_MSS = "MSS"
STR_MSS2 = "MS-Stable"
STR_MSIH = "MSI-H"
STR_MSIL = "MSI-L"
STR_PASS = "Pass"
STR_SUBOPTIMAL = "Sub-optimal"

COL_SAMPLENAME = "SampleName"
COL_TESTNAME = "TestName"
COL_BLOCKNUMBER = "BlockNumber"
COL_SAMPLETYPE = "SampleType"
COL_TISSUEORIGIN = "TissueOrigin"
COL_PATHOLOGICDIAGNOSIS = "PathologicDiagnosis"
COL_TUMORPERCENTAGE = "TumorPercentage"
COL_BONEMARROWASIRATIONDATE = "BoneMarrowAspirationDate"
COL_QCMEANDEPTH = "MeanDepth"
COL_QCCOVERAGE100X = "Coverage100X"
COL_QCUNIFORMITYOVER80 = "UniformityOver80%"
COL_QCONTARGETOVER90 = "OnTargetOver90%"
COL_QCRNASTARTSITESPERCONTROLGSP2OVER10 = "AverageUniqueRNAStartSitesPerControlGSP2Over10"
COL_MSI = "MSI"
COL_TMB = "TMB"
COL_SNPINDEL = "SNP_INDEL"
COL_CNV = "CNV"
COL_HOMODEL = "HomoDEL"
COL_HETERODEL = "HeteroDEL"
COL_FUSION = "Fusion"


COLUME = [COL_SAMPLENAME, COL_TESTNAME, COL_BLOCKNUMBER, COL_SAMPLETYPE, COL_TISSUEORIGIN, COL_PATHOLOGICDIAGNOSIS,
          COL_TUMORPERCENTAGE, COL_BONEMARROWASIRATIONDATE, COL_QCMEANDEPTH, COL_QCCOVERAGE100X, COL_QCUNIFORMITYOVER80,
          COL_QCONTARGETOVER90, COL_QCRNASTARTSITESPERCONTROLGSP2OVER10, COL_MSI, COL_TMB, COL_SNPINDEL, COL_CNV, COL_HOMODEL,
          COL_HETERODEL, COL_FUSION]

SM_ACTG = ['Test Name',
           'Single Nucleotide',
           'Amplification',
           'Homozygous',
           'Heterozygous',
           'Tumor Mutational Burden(TMB)',
           'Microsatellite',
           'Fusion',
           'Diagnosis',
           'MP No.',
           'Sample Type',
           'Block Number',
           'Tissue Origin',
           'Pathologic Diagnosis',
           'Tumor Percentage',
           'NGS QC parameters',
           'TBD',
           'TBD',
           'TBD'
           ]

SM_F1 = ['Test Name',
         'Blood Tumor Mutation',
         'Microsatellite Status',
         'Tumor Mutation Burden',
         'Tumor Fraction',
         'Genomic Finding',
         'Diagnosis',
         'MP No.',
         'Sample Type',
         'Block Number',
         'Tumor Percentage',
         'NGS QC parameters',
         'TBD',
         'TBD',
         'TBD'
         ]

SM_ARCHER = ['Test Name',
             'Gene fusion:',
             'Diagnosis',
             'MP No.',
             'Sample Type',
             'Block Number',
             'Tumor Percentage',
             'NGS QC parameters',
             'TBD',
             'TBD'
             ]

SM_GUARDANT360 = ['Test Name',
                  'Relevant Biomarkers ',
                  'Diagnosis',
                  'MP No.',
                  'Sample Type',
                  'NGS QC parameters',
                  'TBD',
                  'TBD'
                  ]

SM_ONCOMINE = ['Test Name',
               'Biomarkers',
               'Diagnosis',
               'MP No.',
               'Sample Type',
               'NGS QC parameters',
               'Analytic Interpretation',
               'TBD',
               'TBD'
               ]


def usage():
    print("VGH_ReportPreprocess.py -t <Data type> -i <Input data folder> -o <Output data folder>")
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
    print("\tpython ./VGH_ReportPreprocess.py -t ACTOnco+ -i ~/data/VGH-NGS/report/行動基因 -o ~/data/VGH-NGS/table/ACTG")
    print(
        "\tpython ./VGH_ReportPreprocess.py -t F1CDx -i ~/data/VGH-NGS/report/FoundationOne -o ~/data/VGH-NGS/table/FoundationOne")
    print(
        "\tpython ./VGH_ReportPreprocess.py -t ArcherFusionPlex -i ~/data/VGH-NGS/report/ArcherSarcoma -o ~/data/VGH-NGS/table/ArcherSarcoma")
    print(
        "\tpython ./VGH_ReportPreprocess.py -t Guardant360 -i ~/data/VGH-NGS/report/Guardant360 -o ~/data/VGH-NGS/table/Guardant360")
    print(
        "\tpython ./VGH_ReportPreprocess.py -t OncomineBRCA -i ~/data/VGH-NGS/report/OncomineBRCA -o ~/data/VGH-NGS/table/OncomineBRCA")
    print(
        "\tpython ./VGH_ReportPreprocess.py -t OncomineFocus -i ~/data/VGH-NGS/report/OncomineFocus -o ~/data/VGH-NGS/table/OncomineFocus")
    print(
        "\tpython ./VGH_ReportPreprocess.py -t OncomineMyeloid -i ~/data/VGH-NGS/report/OncomineMyeloid -o ~/data/VGH-NGS/table/OncomineMyeloid")
    print(
        "\tpython ./VGH_ReportPreprocess.py -t OncomineTMB -i ~/data/VGH-NGS/report/OncomineTMB -o ~/data/VGH-NGS/table/OncomineTMB")

    return


def extract_samplename(file):
    fn = os.path.basename(file)
    s = fn.find("(")
    e = fn.find(")")

    if 0 < s < e:
        logging.debug("sample name :%s (%d,%d) from %s" % (fn[s + 1:e], s, e, file))
        return fn[s + 1:e].strip()
    else:
        e = fn.find("_")
        if e > 0:
            logging.debug("Sample name :%s (0,%d) from %s" % (fn[:e], e, file))
            return fn[:e].strip()

    logging.warning("Can't find sample name from %s" % fn)

    return ""


def create_record():
    h = defaultdict(str)

    h[COL_SNPINDEL] = STR_NA
    h[COL_CNV] = STR_NA
    h[COL_HOMODEL] = STR_NA
    h[COL_HETERODEL] = STR_NA
    h[COL_FUSION] = STR_NA
    h[COL_TMB] = STR_NA_TMB
    return h


def actonco_preprocess(ifolder, ofolder):
    for file in glob.glob(ifolder + "/*.txt"):
        logging.debug("%s" % file)
        sample_name = extract_samplename(file)
        logging.debug("Sample Name = %s" % sample_name)
        state = -1
        # h_sample = defaultdict(str)
        h_sample = create_record()
        if len(sample_name) > 5:
            h_sample[COL_SAMPLENAME] = sample_name.strip()

        ifd = open(file, "r", encoding='utf-8', errors='ignore')
        for line in ifd:
            if len(line.strip()) == 0:
                continue
            elif line.startswith(SM_ACTG[state + 1]):
                logging.debug(
                    "STATE %d=>%d due to match (%s) from %s" % (state, state + 1, SM_ACTG[state + 1], line.strip()))
                state += 1
                if state == 0:  # Test Name
                    items = line.strip().split(": ", 1)
                    h_sample[COL_TESTNAME] = items[1]
                    logging.debug("Test Name = %s" % h_sample[COL_TESTNAME])
                elif state == 5:  # TMB
                    if line.find("Cannot be determined") >= 0:
                        h_sample[COL_TMB] = STR_NA_TMB
                    else:
                        items = line.strip().split(": ", 1)
                        if items[1].find("<") >= 0:
                            h_sample[COL_TMB] = "0.5"
                        else:
                            h_sample[COL_TMB] = items[1].split(" ")[0]
                    logging.debug("TMB = %f" % float(h_sample[COL_TMB]))
                elif state == 6:  # MSI
                    h_sample[COL_MSI] = "NA"
                    items = line.split(": ", 1)
                    if items[1].find(STR_MSS) >= 0:
                        h_sample[COL_MSI] = STR_MSS
                    elif items[1].find(STR_MSIL) >= 0:
                        h_sample[COL_MSI] = STR_MSIL
                    elif items[1].find(STR_MSIH) >= 0:
                        h_sample[COL_MSI] = STR_MSIH
                    logging.debug("MSI = %s" % h_sample[COL_MSI])
                elif state == 7:  # Fusion
                    if line.find("Not detected") >= 0 or line.find("Not tested") >=0 or \
                            line.find("The fusion sequencing data did not pass the QC criteria") >= 0:
                        h_sample[COL_FUSION] = "NA"
                    else:
                        items = line.strip().split(": ", 1)
                        h_sample[COL_FUSION] = items[1]
                    logging.debug("Fusion = %s" % h_sample[COL_FUSION])
                elif state == 8:    # Diagnosis
                    items = line.strip().split(": ", 1)
                    h_sample[COL_PATHOLOGICDIAGNOSIS] = items[1]
                    logging.debug("%s = %s" % (SM_ACTG[state], h_sample[COL_PATHOLOGICDIAGNOSIS]))
                elif state == 9:    # MP No.
                    items = line.strip().split(": ", 1)
                    if h_sample[COL_SAMPLENAME] != items[1].strip():
                        logging.warning("Sample Name of %s is not %s" % (items[1], h_sample[COL_SAMPLENAME]))
                    h_sample[COL_SAMPLENAME] = items[1].strip()
                    logging.debug("%s = %s" % (SM_ACTG[state], h_sample[COL_SAMPLENAME]))
                elif state == 10:    # Sample Type
                    items = line.strip().split(": ", 1)
                    h_sample[COL_SAMPLETYPE] = items[1]
                    logging.debug("%s = %s" % (SM_ACTG[state], h_sample[COL_SAMPLETYPE]))
                elif state == 11:    # Block Number
                    items = line.strip().split(": ", 1)
                    h_sample[COL_BLOCKNUMBER] = items[1]
                    logging.debug("%s = %s" % (SM_ACTG[state], h_sample[COL_BLOCKNUMBER]))
                elif state == 12:    # Tissue Origin
                    items = line.strip().split(": ", 1)
                    h_sample[COL_TISSUEORIGIN] = items[1]
                    logging.debug("%s = %s" % (SM_ACTG[state], h_sample[COL_TISSUEORIGIN]))
                elif state == 13:   # Pathologic Diagnosis
                    items = line.strip().split(": ", 1)
                    h_sample[COL_PATHOLOGICDIAGNOSIS] = items[1]
                    logging.debug("%s = %s" % (SM_ACTG[state], h_sample[COL_PATHOLOGICDIAGNOSIS]))
                elif state == 14:    # Tumor Percentage
                    items = line.strip().split(": ", 1)
                    h_sample[COL_TUMORPERCENTAGE] = items[1]
                    logging.debug("%s = %s" % (SM_ACTG[state], h_sample[COL_TUMORPERCENTAGE]))
                elif state == 15:   # NGS QC parameters
                    items = line.strip().split(": ", 2)
                    (h_sample[COL_QCMEANDEPTH], h_sample[COL_QCCOVERAGE100X]) = items[2].split("&", 1)
                    logging.debug("%s = %s, %s" % (SM_ACTG[state], h_sample[COL_QCMEANDEPTH],
                                                   h_sample[COL_QCCOVERAGE100X]))
            elif line.startswith(SM_ACTG[state + 2]):
                logging.debug(
                    "STATE %d=>%d due to match (%s) from %s" % (state, state + 2, SM_ACTG[state + 2], line.strip()))
                state += 2
                if state == 9:    # MP No.
                    items = line.strip().split(": ", 1)
                    if h_sample[COL_SAMPLENAME] != items[1].strip():
                        logging.warning("Sample Name of %s is not %s" % (items[1], h_sample[COL_SAMPLENAME]))
                    h_sample[COL_SAMPLENAME] = items[1].strip()
                    logging.debug("%s = %s" % (SM_ACTG[state], h_sample[COL_SAMPLENAME]))
            elif line.startswith(SM_ACTG[state + 3]):
                logging.debug(
                    "STATE %d=>%d due to match (%s) from %s" % (state, state + 3, SM_ACTG[state + 3], line.strip()))
                state += 3
                if state == 14:    # Tumor Percentage
                    items = line.strip().split(": ", 1)
                    h_sample[COL_TUMORPERCENTAGE] = items[1]
                    logging.debug("%s = %s" % (SM_ACTG[state], h_sample[COL_TUMORPERCENTAGE]))
            elif state == 1:
                if line.startswith("Not detected."):
                    h_sample[COL_SNPINDEL] = "NA"
                elif line.startswith("Gene:"):
                    if COL_SNPINDEL not in h_sample or h_sample[COL_SNPINDEL] == STR_NA:
                        h_sample[COL_SNPINDEL] = line.strip()
                    else:
                        h_sample[COL_SNPINDEL] += "; %s" % line.strip()
            elif state == 2:
                if line.startswith("Not detected.") or line.startswith("Cannot be determined"):
                    h_sample[COL_CNV] = "NA"
                elif line.startswith("Chr:"):
                    if COL_CNV not in h_sample or h_sample[COL_CNV] == STR_NA:
                        h_sample[COL_CNV] = line.strip()
                    else:
                        h_sample[COL_CNV] += "; %s" % line.strip()
            elif state == 3:
                if line.startswith("Not detected.") or line.startswith("Cannot be determined"):
                    h_sample[COL_HOMODEL] = "NA"
                elif line.startswith("Chr:"):
                    if COL_HOMODEL not in h_sample or h_sample[COL_HOMODEL] == STR_NA:
                        h_sample[COL_HOMODEL] = line.strip()
                    else:
                        h_sample[COL_HOMODEL] += "; %s" % line.strip()
            elif state == 4:
                if line.startswith("Not detected.") or line.startswith("Cannot be determined"):
                    h_sample[COL_HETERODEL] = "NA"
                elif line.startswith("Chr:"):
                    if COL_HETERODEL not in h_sample or h_sample[COL_HETERODEL] == STR_NA:
                        h_sample[COL_HETERODEL] = line.strip()
                    else:
                        h_sample[COL_HETERODEL] += "; %s" % line.strip()

        ifd.close()

        ofile = "%s/%s.tsv" % (ofolder, h_sample[COL_SAMPLENAME])
        logging.info("%s => %s" % (file, ofile))
        ofd = open(ofile, "w")
        for idx in range(len(COLUME)):
            if idx == 0:
                ofd.write("%s" % h_sample[COLUME[idx]])
            else:
                ofd.write("\t%s" % h_sample[COLUME[idx]])
        ofd.write("\n")
        logging.debug("SNP = %s" % h_sample[COL_SNPINDEL])
        logging.debug("CNV = %s" % h_sample[COL_CNV])
        logging.debug("Homo DEL = %s" % h_sample[COL_HOMODEL])
        logging.debug("Hetero DEL = %s" % h_sample[COL_HETERODEL])
        logging.debug("Fusion = %s" % h_sample[COL_FUSION])
        ofd.close()
        # break
    return


def f1cdx_preprocess(ifolder, ofolder):
    for file in glob.glob(ifolder + "/*.txt"):
        logging.debug("%s" % file)
        sample_name = extract_samplename(file)
        logging.debug("Sample Name = %s" % sample_name)
        state = -1
        # h_sample = defaultdict(str)
        h_sample = create_record()
        if len(sample_name) > 5:
            h_sample[COL_SAMPLENAME] = sample_name

        if sample_name != "PF21018":
            continue

        ifd = open(file, "r", encoding='utf-8', errors='ignore')
        for line in ifd:
            if len(line.strip()) == 0:
                continue
            elif line.startswith(SM_F1[state + 1]):
                logging.debug(
                    "STATE %d=>%d due to match (%s) from %s" % (state, state + 1, SM_F1[state + 1], line.strip()))
                state += 1
                if state == 0:  # Test Name
                    items = line.strip().split(": ", 1)
                    h_sample[COL_TESTNAME] = items[1]
                    logging.debug("Test Name = %s" % h_sample[COL_TESTNAME])
                elif state == 1:  # Blood Tumor Mutation
                    items = line.strip().split("-", 1)
                    h_sample[COL_TMB] = items[1].strip().split(" ")[0]
                    logging.debug("TMB = %f" % float(h_sample[COL_TMB]))
                elif state == 2:  # MSI
                    h_sample[COL_MSI] = "NA"
                    items = line.split("-", 1)
                    if items[1].find(STR_MSS) >= 0 or items[1].find(STR_MSS2) >= 0:
                        h_sample[COL_MSI] = STR_MSS
                    elif items[1].find(STR_MSIL) >= 0:
                        h_sample[COL_MSI] = STR_MSIL
                    elif items[1].find(STR_MSIH) >= 0:
                        h_sample[COL_MSI] = STR_MSIH
                    logging.debug("MSI = [%s]" % h_sample[COL_MSI])
                elif state == 3:  # Tumor Mutational Burden
                    if line.find("Cannot Be Determined") >= 1:
                        h_sample[COL_TMB] = STR_NA_TMB
                    else:
                        items = line.strip().split("- ", 1)
                        h_sample[COL_TMB] = items[1].split(" ")[0]
                    logging.debug("TMB = %f" % float(h_sample[COL_TMB]))
                elif state == 4:  # Tumor Fraction
                    items = line.strip().split("- ", 1)
                    h_sample[COL_TUMORPERCENTAGE] = items[1]
                    logging.debug("Tumor Percentage = %s" % h_sample[COL_TUMORPERCENTAGE])
                elif state == 6:    # Diagnosis
                    items = line.strip().split(": ", 1)
                    h_sample[COL_PATHOLOGICDIAGNOSIS] = items[1]
                    logging.debug("%s = %s" % (SM_F1[state], h_sample[COL_PATHOLOGICDIAGNOSIS]))
                elif state == 7:    # MP No.
                    items = line.strip().split(": ", 1)
                    if h_sample[COL_SAMPLENAME] != items[1].strip():
                        logging.warning("Sample Name of %s is not %s" % (items[1], h_sample[COL_SAMPLENAME]))
                    h_sample[COL_SAMPLENAME] = items[1].strip()
                    logging.debug("%s = %s" % (SM_F1[state], h_sample[COL_SAMPLENAME]))
                elif state == 8:    # Sample Type
                    items = line.strip().split(": ", 1)
                    h_sample[COL_SAMPLETYPE] = items[1]
                    logging.debug("%s = %s" % (SM_F1[state], h_sample[COL_SAMPLETYPE]))
                elif state == 9:    # Block Number
                    items = line.strip().split(": ", 1)
                    h_sample[COL_BLOCKNUMBER] = items[1]
                    logging.debug("%s = %s" % (SM_F1[state], h_sample[COL_BLOCKNUMBER]))
                elif state == 10:    # Tumor Percentage
                    items = line.strip().split(": ", 1)
                    h_sample[COL_TUMORPERCENTAGE] = items[1]
                    logging.debug("%s = %s" % (SM_F1[state], h_sample[COL_TUMORPERCENTAGE]))
                elif state == 11:   # NGS QC parameters
                    items = line.strip().split(": ", 1)
                    if items[1].find("DNA input mass") >= 0:
                        continue
                    else:
                        h_sample[COL_QCMEANDEPTH] = items[1].split(" ")[-1]
                    logging.debug("%s = %s" % (SM_F1[state], h_sample[COL_QCMEANDEPTH]))
            elif line.startswith(SM_F1[state + 2]):
                logging.debug(
                    "STATE %d=>%d due to match (%s) from %s" % (state, state + 2, SM_F1[state + 2], line.strip()))
                state += 2
                if state == 2:  # MSI
                    h_sample[COL_MSI] = "NA"
                    items = line.split("-", 1)
                    if items[1].find(STR_MSS) >= 0 or items[1].find(STR_MSS2) >= 0:
                        h_sample[COL_MSI] = STR_MSS
                    elif items[1].find(STR_MSIL) >= 0:
                        h_sample[COL_MSI] = STR_MSIL
                    elif items[1].find(STR_MSIH) >= 0:
                        h_sample[COL_MSI] = STR_MSIH
                    logging.debug("MSI = [%s]" % h_sample[COL_MSI])
                elif state == 4:  #Tumor Fraction
                    items = line.strip().split("- ", 1)
                    h_sample[COL_TUMORPERCENTAGE] = items[1]
                    logging.debug("Tumor Percentage = %s" % h_sample[COL_TUMORPERCENTAGE])
                elif state == 7:    # MP No.
                    items = line.strip().split(": ", 1)
                    if h_sample[COL_SAMPLENAME] != items[1].strip():
                        logging.warning("Sample Name of %s is not %s" % (items[1], h_sample[COL_SAMPLENAME]))
                    h_sample[COL_SAMPLENAME] = items[1].strip()
                    logging.debug("%s = %s" % (SM_F1[state], h_sample[COL_SAMPLENAME]))
                elif state == 11:  # NGS QC parameters
                    items = line.strip().split(": ", 1)
                    if items[1].find("DNA input mass") >= 0:
                        continue
                    else:
                        h_sample[COL_QCMEANDEPTH] = items[1].split(" ")[-1]
                    logging.debug("%s = %s" % (SM_F1[state], h_sample[COL_QCMEANDEPTH]))
            elif line.startswith(SM_F1[state + 3]):
                logging.debug(
                        "STATE %d=>%d due to match (%s) from %s" % (state, state + 2, SM_F1[state + 2], line.strip()))
                state += 3
            elif state == 5:
                if line.startswith("For a complete list of ") or line.startswith("Note") \
                        or line.startswith("Please refer to Picture"):
                    continue
                else:
                    if COL_SNPINDEL not in h_sample or h_sample[COL_SNPINDEL] == STR_NA:
                        h_sample[COL_SNPINDEL] = line.strip()
                    else:
                        h_sample[COL_SNPINDEL] += "; %s" % line.strip()

        ifd.close()

        ofile = "%s/%s.tsv" % (ofolder, h_sample[COL_SAMPLENAME])
        logging.info("%s => %s" % (file, ofile))
        ofd = open(ofile, "w")
        for idx in range(len(COLUME)):
            if idx == 0:
                ofd.write("%s" % h_sample[COLUME[idx]])
            else:
                ofd.write("\t%s" % h_sample[COLUME[idx]])
        ofd.write("\n")
        logging.debug("MSI = %s" % h_sample[COL_MSI])
        logging.debug("TMB = %s" % h_sample[COL_TMB])
        logging.debug("SNP = %s" % h_sample[COL_SNPINDEL])
        logging.debug("CNV = %s" % h_sample[COL_CNV])
        logging.debug("Homo DEL = %s" % h_sample[COL_HOMODEL])
        logging.debug("Hetero DEL = %s" % h_sample[COL_HETERODEL])
        logging.debug("Fusion = %s" % h_sample[COL_FUSION])
        ofd.close()
        # break
    return


def archer_preprocess(ifolder, ofolder):
    for file in glob.glob(ifolder + "/*.txt"):
        logging.debug("%s" % file)
        sample_name = extract_samplename(file)
        logging.debug("Sample Name = %s" % sample_name)
        state = -1
        h_sample = create_record()
        if len(sample_name) > 5:
            h_sample[COL_SAMPLENAME] = sample_name

        ifd = open(file, "r", encoding='utf-8', errors='ignore')
        for line in ifd:
            if len(line.strip()) == 0:
                continue
            elif line.startswith(SM_ARCHER[state + 1]):
                logging.debug(
                    "STATE %d=>%d due to match (%s) from %s" % (state, state + 1, SM_ARCHER[state + 1], line.strip()))
                state += 1
                if state == 0:  # Test Name
                    items = line.strip().split(": ", 1)
                    h_sample[COL_TESTNAME] = items[1]
                    logging.debug("Test Name = %s" % h_sample[COL_TESTNAME])
                elif state == 2:  # Diagnosis
                    items = line.strip().split(": ", 1)
                    h_sample[COL_PATHOLOGICDIAGNOSIS] = items[1]
                    logging.debug("%s = %s" % (SM_ARCHER[state], h_sample[COL_PATHOLOGICDIAGNOSIS]))
                elif state == 3:  # MP No.
                    items = line.strip().split(": ", 1)
                    if h_sample[COL_SAMPLENAME] != items[1].strip():
                        logging.warning("Sample Name of %s is not %s" % (items[1], h_sample[COL_SAMPLENAME]))
                    h_sample[COL_SAMPLENAME] = items[1].strip()
                    logging.debug("%s = %s" % (SM_ARCHER[state], h_sample[COL_SAMPLENAME]))
                elif state == 4:  # Sample Type
                    items = line.strip().split(": ", 1)
                    h_sample[COL_SAMPLETYPE] = items[1]
                    logging.debug("%s = %s" % (SM_ARCHER[state], h_sample[COL_SAMPLETYPE]))
                elif state == 5:  # Block Number
                    items = line.strip().split(": ", 1)
                    h_sample[COL_BLOCKNUMBER] = items[1]
                    logging.debug("%s = %s" % (SM_ARCHER[state], h_sample[COL_BLOCKNUMBER]))
                elif state == 6:  # Tumor Percentage
                    items = line.strip().split(": ", 1)
                    h_sample[COL_TUMORPERCENTAGE] = items[1]
                    logging.debug("%s = %s" % (SM_ARCHER[state], h_sample[COL_TUMORPERCENTAGE]))
                elif state == 7:    #NGS QC parameters
                    h_sample[COL_QCRNASTARTSITESPERCONTROLGSP2OVER10] = STR_PASS
            elif state == 1:
                if line.startswith("No fusion variants") or line.startswith("Note") \
                        or line.startswith("Please refer to Picture"):
                    continue
                else:
                    if COL_FUSION not in h_sample or h_sample[COL_FUSION] == STR_NA:
                        h_sample[COL_FUSION] = line.strip()
                    else:
                        h_sample[COL_FUSION] += "; %s" % line.strip()
            elif state == 7:  # NGS QC parameters
                if line.find("[  ] Pass") >= 0:
                    h_sample[COL_QCRNASTARTSITESPERCONTROLGSP2OVER10] = STR_SUBOPTIMAL
                    logging.debug("%s = %s" % (SM_ARCHER[state], h_sample[COL_QCRNASTARTSITESPERCONTROLGSP2OVER10]))

        ifd.close()

        ofile = "%s/%s.tsv" % (ofolder, h_sample[COL_SAMPLENAME])
        logging.info("%s => %s" % (file, ofile))
        ofd = open(ofile, "w")
        for idx in range(len(COLUME)):
            if idx == 0:
                ofd.write("%s" % h_sample[COLUME[idx]])
            else:
                ofd.write("\t%s" % h_sample[COLUME[idx]])
        ofd.write("\n")
        logging.debug("SNP = %s" % h_sample[COL_SNPINDEL])
        logging.debug("CNV = %s" % h_sample[COL_CNV])
        logging.debug("Homo DEL = %s" % h_sample[COL_HOMODEL])
        logging.debug("Hetero DEL = %s" % h_sample[COL_HETERODEL])
        logging.debug("Fusion = %s" % h_sample[COL_FUSION])
        ofd.close()
        # break
    return


def guardant_preprocess(ifolder, ofolder):
    for file in glob.glob(ifolder + "/*.txt"):
        logging.debug("%s" % file)
        sample_name = extract_samplename(file)
        logging.debug("Sample Name = %s" % sample_name)
        state = -1
        sub_state = -1
        h_sample = create_record()
        if len(sample_name) > 5:
            h_sample[COL_SAMPLENAME] = sample_name

        ifd = open(file, "r", encoding='utf-8', errors='ignore')
        for line in ifd:
            if len(line.strip()) == 0:
                continue
            elif line.startswith(SM_GUARDANT360[state + 1]):
                logging.debug(
                    "STATE %d=>%d due to match (%s) from %s" % (state, state + 1, SM_GUARDANT360[state + 1], line.strip()))
                state += 1
                if state == 0:  # Test Name
                    items = line.strip().split(": ", 1)
                    h_sample[COL_TESTNAME] = items[1]
                    logging.debug("Test Name = %s" % h_sample[COL_TESTNAME])
                elif state == 2:  # Diagnosis
                    items = line.strip().split(": ", 1)
                    h_sample[COL_PATHOLOGICDIAGNOSIS] = items[1]
                    logging.debug("%s = %s" % (SM_GUARDANT360[state], h_sample[COL_PATHOLOGICDIAGNOSIS]))
                elif state == 3:  # MP No.
                    items = line.strip().split(": ", 1)
                    if len(items) <= 1:     # MP No.:PG21005
                        items = line.strip().split(":", 1)
                    if h_sample[COL_SAMPLENAME] != items[1].strip():
                        logging.warning("Sample Name of %s is not %s" % (items[1], h_sample[COL_SAMPLENAME]))
                    h_sample[COL_SAMPLENAME] = items[1].strip()
                    logging.debug("%s = %s" % (SM_GUARDANT360[state], h_sample[COL_SAMPLENAME]))
                elif state == 4:  # Sample Type
                    items = line.strip().split(": ", 1)
                    h_sample[COL_SAMPLETYPE] = items[1]
                    logging.debug("%s = %s" % (SM_GUARDANT360[state], h_sample[COL_SAMPLETYPE]))
                elif state == 5:    #NGS QC parameters
                    h_sample[COL_QCMEANDEPTH] = "15,000x"
            elif state == 1:
                logging.debug("state = 1: %s" % line.strip())
                if line.startswith("Note") or line.startswith("Please refer to Picture"):
                    continue
                elif line.find("Gene or Biomarker") >= 0:
                    sub_state = 1
                elif sub_state == 1:
                    if line.find("MSI") >= 0:
                        if line.find("NOT DETECTED"):
                            h_sample[COL_MSI] = STR_MSS
                        else:
                            h_sample[COL_MSI] = STR_MSIH
                    elif line.find("Amplification") >= 0:
                        if COL_CNV not in h_sample or h_sample[COL_CNV] == STR_NA:
                            h_sample[COL_CNV] = line.strip().replace("\t", " ")
                        else:
                            h_sample[COL_CNV] += "; %s" % line.strip().replace("\t", " ")
                    elif line.find("Deletion") >= 0:
                        if COL_HETERODEL not in h_sample or h_sample[COL_HETERODEL] == STR_NA:
                            h_sample[COL_HETERODEL] = line.strip().replace("\t", " ")
                        else:
                            h_sample[COL_HETERODEL] += "; %s" % line.strip().replace("\t", " ")
                    elif COL_SNPINDEL not in h_sample or h_sample[COL_SNPINDEL] == STR_NA:
                        h_sample[COL_SNPINDEL] = line.strip().replace("\t", " ")
                    else:
                        h_sample[COL_SNPINDEL] += "; %s" % line.strip().replace("\t", " ")
                elif sub_state == 2:
                    if line.find("NOT DETECTED") >= 0:
                        h_sample[COL_MSI] = STR_MSS
                    else:
                        h_sample[COL_MSI] = STR_MSIH
                elif line.find("MSI") >= 0:
                    sub_state = 2
                else:
                    if COL_SNPINDEL not in h_sample or h_sample[COL_SNPINDEL] == STR_NA:
                        h_sample[COL_SNPINDEL] = line.strip().replace("\t", " ")
                    elif line.find("Gene") >= 0:
                        h_sample[COL_SNPINDEL] += "; %s" % line.strip().replace("\t", " ")
                    else:
                        h_sample[COL_SNPINDEL] += ", %s" % line.strip().replace("\t", " ")
            elif state == 5:  # NGS QC parameters
                if line.find("Average depth") >= 0:
                    items = line.split(".", 1)
                    h_sample[COL_QCMEANDEPTH] = items[0].split(" ")[-1]
                    logging.debug("%s = %s" % (SM_GUARDANT360[state], h_sample[COL_QCMEANDEPTH]))

        ifd.close()

        ofile = "%s/%s.tsv" % (ofolder, h_sample[COL_SAMPLENAME])
        logging.info("%s => %s" % (file, ofile))
        ofd = open(ofile, "w")
        for idx in range(len(COLUME)):
            if idx == 0:
                ofd.write("%s" % h_sample[COLUME[idx]])
            else:
                ofd.write("\t%s" % h_sample[COLUME[idx]])
        ofd.write("\n")
        logging.debug("SNP = %s" % h_sample[COL_SNPINDEL])
        logging.debug("CNV = %s" % h_sample[COL_CNV])
        logging.debug("Homo DEL = %s" % h_sample[COL_HOMODEL])
        logging.debug("Hetero DEL = %s" % h_sample[COL_HETERODEL])
        logging.debug("Fusion = %s" % h_sample[COL_FUSION])
        ofd.close()
        # break
    return


def oncomine_preprocess(ifolder, ofolder):
    for file in glob.glob(ifolder + "/*.txt"):
        logging.debug("%s" % file)
        sample_name = extract_samplename(file)
        logging.debug("Sample Name = %s" % sample_name)
        state = -1
        sub_state = -1
        h_sample = create_record()
        if len(sample_name) > 5:
            h_sample[COL_SAMPLENAME] = sample_name

        ifd = open(file, "r", encoding='utf-8', errors='ignore')
        for line in ifd:
            if len(line.strip()) == 0:
                continue
            elif line.find(SM_ONCOMINE[state + 1]) >= 0:
                logging.debug(
                    "STATE %d=>%d due to match (%s) from %s" % (state, state + 1, SM_ONCOMINE[state + 1], line.strip()))
                state += 1
                if state == 0:  # Test Name
                    items = line.strip().split(": ", 1)
                    h_sample[COL_TESTNAME] = items[1]
                    logging.debug("Test Name = %s" % h_sample[COL_TESTNAME])
                elif state == 2:  # Diagnosis
                    items = line.strip().split(": ", 1)
                    h_sample[COL_PATHOLOGICDIAGNOSIS] = items[1]
                    logging.debug("%s = %s" % (SM_ONCOMINE[state], h_sample[COL_PATHOLOGICDIAGNOSIS]))
                elif state == 3:  # MP No.
                    items = line.strip().split(": ", 1)
                    if len(items) <= 1:     # MP No.:PG21005
                        items = line.strip().split(":", 1)
                    if h_sample[COL_SAMPLENAME] != items[1].strip().replace("/ ", "_").replace("/", "_"):
                        logging.debug("Sample Name of %s is not %s" % (items[1], h_sample[COL_SAMPLENAME]))
                        if len(h_sample[COL_SAMPLENAME]) <= len(items[1].strip().replace("/ ", "_").replace("/", "_")):
                            h_sample[COL_SAMPLENAME] = items[1].strip().replace("/ ", "_").replace("/", "_")

                    # NOTE： fix Wrong Sample Name (e.g. F1901 ==> F19001)
                    no = re.search(r'(\D+)(\d+)', h_sample[COL_SAMPLENAME], re.M| re.I)
                    if len(no.group(2)) == 4:
                        logging.debug("SampleName[%s] has %d digits [%s], %d chars [%s]" % (
                        h_sample[COL_SAMPLENAME], len(no.group(2)), no.group(2), len(no.group(1)), no.group(1)))
                        h_sample[COL_SAMPLENAME] = "%s%s0%s" % (no.group(1), no.group(2)[:2], no.group(2)[2:])
                    logging.debug("%s = %s" % (SM_ONCOMINE[state], h_sample[COL_SAMPLENAME]))
                elif state == 4:  # Sample Type
                    items = line.strip().split(": ", 1)
                    h_sample[COL_SAMPLETYPE] = items[1]
                    logging.debug("%s = %s" % (SM_ONCOMINE[state], h_sample[COL_SAMPLETYPE]))
                elif state == 5:    #NGS QC parameters
                    h_sample[COL_QCMEANDEPTH] = "500x"
                    h_sample[COL_QCUNIFORMITYOVER80] = STR_PASS
                    h_sample[COL_QCONTARGETOVER90] = STR_PASS
                elif state == 6:    # Analytic Interpretation
                    break
            elif state == 1:
                logging.debug("state = 1 (sub_state=%d): %s" % (sub_state, line.strip()))
                if line.startswith("(Gene") or line.startswith("Note") or line.startswith("--")\
                        or line.startswith("Please refer to Picture") or line.startswith("Relevant Biomarkers")\
                        or line.startswith('No relevant biomarkers found in this sample.'):
                    sub_state = -1
                    continue
                elif line.find("Tumor Mutational Burden") >= 0:
                    sub_state = -1
                    items = line.strip().split(":")[1].strip().split(" ")
                    h_sample[COL_TMB] = items[0].strip()
                    logging.debug("TMB=[%f]" % float(h_sample[COL_TMB]))
                elif line.find("deletion") >= 0:
                    sub_state = 1
                    if COL_HETERODEL not in h_sample or h_sample[COL_HETERODEL] == STR_NA:
                        h_sample[COL_HETERODEL] = line.strip()
                    else:
                        h_sample[COL_HETERODEL] += "; %s" % line.strip()
                elif line.find("amplification") >= 0:
                    sub_state = 2
                    if COL_CNV not in h_sample or h_sample[COL_CNV] == STR_NA:
                        h_sample[COL_CNV] = line.strip()
                    else:
                        h_sample[COL_CNV] += "; %s" % line.strip()
                elif line.find("fusion") >= 0:
                    sub_state = 3
                    if COL_FUSION not in h_sample or h_sample[COL_FUSION] == STR_NA:
                        h_sample[COL_FUSION] = line.strip()
                    else:
                        h_sample[COL_FUSION] += "; %s" % line.strip()
                elif line.find("Gene:") >= 0:
                    sub_state = -1
                    if COL_SNPINDEL not in h_sample or h_sample[COL_SNPINDEL] == STR_NA:
                        h_sample[COL_SNPINDEL] = line.strip()
                    elif line.find("Gene") >= 0:
                        h_sample[COL_SNPINDEL] += "; %s" % line.strip()
                    else:
                        h_sample[COL_SNPINDEL] += ", %s" % line.strip()
                elif sub_state == 1:
                    if COL_HETERODEL not in h_sample or h_sample[COL_HETERODEL] == STR_NA:
                        h_sample[COL_HETERODEL] = line.strip()
                    else:
                        h_sample[COL_HETERODEL] += "; %s" % line.strip()
                elif sub_state == 2:
                    if COL_CNV not in h_sample or h_sample[COL_CNV] == STR_NA:
                        h_sample[COL_CNV] = line.strip()
                    else:
                        h_sample[COL_CNV] += "; %s" % line.strip()
                elif sub_state == 3:
                    if COL_FUSION not in h_sample or h_sample[COL_FUSION] == STR_NA:
                        h_sample[COL_FUSION] = line.strip()
                    else:
                        h_sample[COL_FUSION] += "; %s" % line.strip()
                else:
                    sub_state = -1
                    if COL_SNPINDEL not in h_sample or h_sample[COL_SNPINDEL] == STR_NA:
                        h_sample[COL_SNPINDEL] = line.strip()
                    elif line.find("Gene") >= 0:
                        h_sample[COL_SNPINDEL] += "; %s" % line.strip()
                    else:
                        h_sample[COL_SNPINDEL] += ", %s" % line.strip()
            elif state == 5:  # NGS QC parameters
                if line.find("[  ] Pass") >= 0:
                    if line.find("Read Depth") >= 0:
                        h_sample[COL_QCMEANDEPTH] = "100x"
                    elif line.find("Uniformity over") >= 0:
                        h_sample[COL_QCUNIFORMITYOVER80] = STR_SUBOPTIMAL
                    elif line.find("On target over") >= 0:
                        h_sample[COL_QCONTARGETOVER90] = STR_SUBOPTIMAL

        ifd.close()

        ofile = "%s/%s.tsv" % (ofolder, h_sample[COL_SAMPLENAME])
        logging.info("%s => %s" % (file, ofile))
        ofd = open(ofile, "w")
        for idx in range(len(COLUME)):
            if idx == 0:
                ofd.write("%s" % h_sample[COLUME[idx]])
            else:
                ofd.write("\t%s" % h_sample[COLUME[idx]])
        ofd.write("\n")
        logging.debug("SNP = %s" % h_sample[COL_SNPINDEL])
        logging.debug("CNV = %s" % h_sample[COL_CNV])
        logging.debug("Homo DEL = %s" % h_sample[COL_HOMODEL])
        logging.debug("Hetero DEL = %s" % h_sample[COL_HETERODEL])
        logging.debug("Fusion = %s" % h_sample[COL_FUSION])
        logging.debug("TMB = %s" % h_sample[COL_TMB])
        ofd.close()
        # break
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

    if not os.path.isdir(ofolder):
        os.makedirs(ofolder, exist_ok=True)

    if panel_type == STR_ACTOncoPlus:
        actonco_preprocess(ifolder, ofolder)
    elif panel_type == STR_F1CDx:
        f1cdx_preprocess(ifolder, ofolder)
    elif panel_type == STR_ArcherFusionPlex:
        archer_preprocess(ifolder, ofolder)
    elif panel_type == STR_Guardant360:
        guardant_preprocess(ifolder, ofolder)
    elif panel_type == STR_OncomineBRCA or panel_type == STR_OncomineFocus \
            or panel_type == STR_OncomineMyeloid or panel_type == STR_OncomineTMB:
        oncomine_preprocess(ifolder, ofolder)

    return


if __name__ == '__main__':
    main(sys.argv[1:])
