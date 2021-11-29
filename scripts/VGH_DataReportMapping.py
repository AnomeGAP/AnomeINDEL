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
import glob
import subprocess

import re
from collections import defaultdict
import logging

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
##INFO=<ID=EQUIVOCAL,Number=.,Type=String,Description="the equivocal of the SV">
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

PANEL_FOLDER = ['ACTG',
                'FoundationOne',
                'ArcherSarcoma',
                'Guardant360',
                'OncomineBRCA',
                'OncomineFocus',
                'OncomineMyeloid',
                'OncomineTMB'
                ]



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



def usage():
    print("VGH_DataReportMapping.py -i <Input Folder> -o <Output VCF and Table>")
    print("Argument:")
    print("\t-h: Usage")
    print("\t-i: input folder")
    print("\t-o: output file (.vcf and Table)")
    print("Usage:")
    print("\tpython ./VGH_DataReportMapping.py -i ~/data/VGH-NGS -o ~/data/VGH-NGS/output ")

    return


def mapping(ifolder, ofolder):
    os.makedirs("%s/vcf" % ofolder, exist_ok=True)
    merge_cmd = "bcftools merge -O z -o %s/vgh.vcf.gz " % ofolder
    ofd = open("%s/vgh.tsv" % ofolder, "w")
    ofd.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %
              (COL_SAMPLENAME, COL_TESTNAME, COL_BLOCKNUMBER, COL_SAMPLETYPE, COL_TISSUEORIGIN,
               COL_PATHOLOGICDIAGNOSIS, COL_TUMORPERCENTAGE, COL_BONEMARROWASIRATIONDATE, COL_QCMEANDEPTH, COL_QCCOVERAGE100X,
               COL_QCUNIFORMITYOVER80, COL_QCONTARGETOVER90, COL_QCRNASTARTSITESPERCONTROLGSP2OVER10, COL_MSI, COL_TMB,
               COL_SNPINDEL, COL_CNV, COL_HOMODEL, COL_HETERODEL, COL_FUSION))
    for panel in PANEL_FOLDER:
        count = 0
        for file in glob.glob("%s/table/%s/*.tsv" % (ifolder, panel)):
            fn = os.path.basename(file)
            if panel == "OncomineTMB" and fn.find("_") >= 0:
                sample_name = fn.split(".", 1)[0].split("_")[1]
            else:
                sample_name = fn.split(".", 1)[0].split("_")[0]
            logging.debug("\t%s" % sample_name)
            vcf = "%s/vcf/%s/%s.vcf" % (ifolder, panel, sample_name)
            os.makedirs("%s/gz/%s" % (ofolder, panel), exist_ok=True)
            gz = "%s/gz/%s/%s.vcf.gz" % (ofolder, panel, sample_name)
            if not os.path.exists(vcf):
                logging.info("%s is not existed from %s" % (vcf, file))
                if panel == "ArcherSarcoma" or panel == 'Guardant360':
                    # create an empty vcf
                    os.makedirs("%s/vcf/%s" % (ifolder, panel), exist_ok=True)
                    ovcf = open(vcf, "w")
                    ovcf.write(VCF_HEADER % sample_name)
                    ovcf.close()
                else:
                    continue
            count += 1

            os.makedirs("%s/vcf/%s" % (ofolder, panel), exist_ok=True)
            new_vcf = "%s/vcf/%s/%s" % (ofolder, panel, os.path.basename(vcf))
            new_vfd = open(new_vcf, "w")
            vfd = open(vcf, "r")
            vname = ""
            for line in vfd:
                if line.startswith("#CHROM"):
                    items = line.strip().split("\t")
                    vname = items[-1]
                elif line.startswith("##INFO=<ID=SUBSET,") and panel == "OncomineFocus":
                    logging.warning("%s **************************************  SOURCE correct" % fn)
                    new_vfd.write("##INFO=<ID=SUBSET,Number=A,Type=String,Description=\"1-based index in ALT list of genotyped allele(s) that are a strict superset\">\n")
                    continue
                elif line.startswith("##INFO=<ID=FR,Number=A") and panel == "OncomineTMB":
                    new_vfd.write("##INFO=<ID=FR,Number=.,Type=String,Description=\"Reason why the variant was filtered.\">\n")
                    continue
                elif line.startswith("##INFO=<ID=LEN,") and panel == "OncomineBRCA":
                    new_vfd.write("##INFO=<ID=LOW_QUALITY_SCORE,Number=.,Type=String,Description=\"Missing in the original vcf(added by ATGNEOMIX)\">\n")
                elif line.startswith("##INFO=<ID=LEN,Number=A,") and panel == "ACTG":
                    new_vfd.write("##INFO=<ID=SOURCE,Number=1,Type=String,Description=\"Missing in the original vcf(added by ATGNEOMIX)\">\n")
                new_vfd.write("%s" % line)
            new_vfd.close()
            vfd.close()
            tfd = open(file, "r")
            tname = ""
            for line in tfd:
                if not line.startswith("#"):
                    items = line.strip().split("\t", 1)
                    tname = items[0]
                    ofd.write("%s\t%s\n" % (vname, items[1]))
                    break
            tfd.close()
            if vname != tname:
                logging.info("Name Inconsistency: [%s] from report but [%s] from data" % (tname, vname))

            cmd = "bcftools sort -O z -o %s %s" % (gz, new_vcf)
            logging.debug("cmd = %s" % cmd)
            try:
                res = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
                (output, err) = res.communicate()
                p_status = res.wait()
            except IOError:
                logging.ERROR("WARNING: %s " % res.stderr)
            cmd = "bcftools index %s " % gz
            logging.debug("cmd = %s" % cmd)
            try:
                res = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
                (output, err) = res.communicate()
                p_status = res.wait()
            except IOError:
                logging.ERROR("WARNING: %s " % res.stderr)
            merge_cmd += " %s" % gz
        # break
        logging.info("%s has %d samples" % (panel, count))
    ofd.close()
    try:
        logging.debug("cmd = %s" % merge_cmd)
        res = subprocess.Popen(merge_cmd, stdout=subprocess.PIPE, shell=True)
        (output, err) = res.communicate()
        p_status = res.wait()
    except IOError:
        logging.ERROR("WARNING: %s " % res.stderr)

    return


def main(argv):
    ifolder = ""
    ofolder = ""

    logging.basicConfig(level=os.environ.get("LOGLEVEL", "WARN"))

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
            ifolder = arg
            if ofolder == "":
                ofolder = "%s/output" % ifolder
        elif opt == '-o':
            ofolder = arg

    if ifolder == "":
        print("Error: '-i' is required")
        usage()
        sys.exit(2)
    elif not os.path.isdir(ifolder):
        print("Error: input folder(%s) is not existed" % ifolder)
        usage()
        sys.exit(3)

    mapping(ifolder, ofolder)

    return


if __name__ == '__main__':
    main(sys.argv[1:])
