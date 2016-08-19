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
# @file    AGCFilter.py
#
# @brief   A program to statistic Anome Gene Curated Database
#
# @author  A-Tsai Su(chungtsai_su@anome.com.tw)
#
# @date    2016/07/01
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
    """ Description: Program Usage
        Argument:    NONE
        Return:	     NONE
    """
    print(
        "AGCFilter.py -i <AGC database file> -o <Output file> -t <minimal total supports> -p <minimal supports for pathogenic> -b <minimal supports for benign")
    print("Argument:")
    print("\t-h: Usage")
    print("\t-i: <AGC database file>")
    print("\t-o: <Output file>")
    print("\t-t: <minimal total supports")
    print("\t-p: <minimal supports for pathogenic>")
    print("\t-b: <minimal supports for benign>")

    print("Usage:")
    print("\tpython ~/src/github/AnomeINDEL/scripts/ACGFilter.py -i ~/gene.tsv -o ~/gene.output -p 10 -b 100")
    print("\tpython ./AGCFilter.py -i ~/gene_005_1000.tsv -t 10 -p 10 -b 100 -o /Users/chungtsai_su/data/AGCD/output.tsv")

    return


MAX_COUNT = 100
MIN_TOTAL = 10
MIN_PATHOGENIC = 10
MIN_BENIGN = 100
MAX_PSV1_BENIGN = 3
MIN_PSV1_PATHO = 5
MIN_TOTAL_MISSENSE = 10
MAX_PP2_THRESHOLD = 0.1
MIN_BP1_THRESHOLD = 0.9

# 0-8  gene_name	pathogenic	uncertain	benign	amr_common	afr_common	eur_common	eas_common	sas_common
# 9-21 synonymous	missense	nonsense	frameshift	start_loss	splicing_site	long_deletion	sv_in_splicing  stop_loss	del	ins	snp	mnp
# 22-34    pathogenic_synonymous	pathogenic_missense	pathogenic_nonsense	pathogenic_frameshift	pathogenic_start_loss	pathogenic_splicing_site	pathogenic_long_deletion	pathogenic_sv_in_splicing	pathogenic_stop_loss	pathogenic_del	pathogenic_ins	pathogenic_snp	pathogenic_mnp
# 35-47    uncertain_synonymous	uncertain_missense	uncertain_nonsense	uncertain_frameshift	uncertain_start_loss	uncertain_splicing_site	uncertain_long_deletion	uncertain_sv_in_splicing	uncertain_stop_loss	uncertain_del	uncertain_ins	uncertain_snp	uncertain_mnp
# 48-60    benign_synonymous	benign_missense	benign_nonsense	benign_frameshift	benign_start_loss	benign_splicing_site	benign_long_deletion	benign_sv_in_splicing	benign_stop_loss	benign_del	benign_ins	benign_snp	benign_mnp
# 61-73    amr_common_synonymous	amr_common_missense	amr_common_nonsense	amr_common_frameshift	amr_common_start_loss	amr_common_splicing_site	amr_common_long_deletion	amr_common_sv_in_splicing	amr_common_stop_loss	amr_common_del	amr_common_ins	amr_common_snp	amr_common_mnp
# 74-86    afr_common_synonymous	afr_common_missense	afr_common_nonsense	afr_common_frameshift	afr_common_start_loss	afr_common_splicing_site	afr_common_long_deletion	afr_common_sv_in_splicing	afr_common_stop_loss	afr_common_del	afr_common_ins	afr_common_snp	afr_common_mnp
# 87-99    eur_common_synonymous	eur_common_missense	eur_common_nonsense	eur_common_frameshift	eur_common_start_loss	eur_common_splicing_site	eur_common_long_deletion	eur_common_sv_in_splicing	eur_common_stop_loss	eur_common_del	eur_common_ins	eur_common_snp	eur_common_mnp
# 100-112   eas_common_synonymous	eas_common_missense	eas_common_nonsense	eas_common_frameshift	eas_common_start_loss	eas_common_splicing_site	eas_common_long_deletion	eas_common_sv_in_splicing	eas_common_stop_loss	eas_common_del	eas_common_ins	eas_common_snp	eas_common_mnp
# 113-125  sas_common_synonymous	sas_common_missense	sas_common_nonsense	sas_common_frameshift	sas_common_start_loss	sas_common_splicing_site	sas_common_long_deletion	sas_common_sv_in_splicing	sas_common_stop_loss	sas_common_del	sas_common_ins	sas_common_snp	sas_common_mnp
def enum(*sequential, **named):
    enums = dict(zip(sequential, range(len(sequential))), **named)
    return type('Enum', (), enums)


COLUMN = enum('GENE_NAME', 'PATHOGENIC', 'UNCERTAIN', 'BENIGN', 'AMR_COMMON', 'AFR_COMMON', 'EUR_COMMON', 'EAS_COMMON',
              'SAS_COMMON',
              'SYNONYMOUS', 'MISSENSE', 'NONSENSE', 'FRAMESHIFT', 'START_LOSS', 'SPLICING_SITE', 'LONG_DELETION',
              'SV_IN_SPLICING', 'STOP_LOSS', 'DEL', 'INS', 'SNP', 'MNP',
              'PATHOGENIC_SYNONYMOUS', 'PATHOGENIC_MISSENSE', 'PATHOGENIC_NONSENSE', 'PATHOGENIC_FRAMESHIFT',
              'PATHOGENIC_START_LOSS', 'PATHOGENIC_SPLICING_SITE', 'PATHOGENIC_LONG_DELETION', 'PATHOGENIC_SV_IN_SPLICING', 'PATHOGENIC_STOP_LOSS',
              'PATHOGENIC_DEL', 'PATHOGENIC_INS', 'PATHOGENIC_SNP', 'PATHOGENIC_MNP',
              'UNCERTAIN_SYNONYMOUS', 'UNCERTAIN_MISSENSE', 'UNCERTAIN_NONSENSE', 'UNCERTAIN_FRAMESHIFT',
              'UNCERTAIN_START_LOSS', 'UNCERTAIN_SPLICING_SITE', 'UNCERTAIN_LONG_DELETION', 'UNCERTAIN_SV_IN_SPLICING', 'UNCERTAIN_STOP_LOSS',
              'UNCERTAIN_DEL', 'UNCERTAIN_INS', 'UNCERTAIN_SNP', 'UNCERTAIN_MNP',
              'BENIGN_SYNONYMOUS', 'BENIGN_MISSENSE', 'BENIGN_NONSENSE', 'BENIGN_FRAMESHIFT', 'BENIGN_START_LOSS',
              'BENIGN_SPLICING_SITE', 'BENIGN_LONG_DELETION', 'BENIGN_SV_IN_SPLICING', 'BENIGN_STOP_LOSS', 'BENIGN_DEL', 'BENIGN_INS',
              'BENIGN_SNP', 'BENIGN_MNP',
              'AMR_COMMON_SYNONYMOUS', 'AMR_COMMON_MISSENSE', 'AMR_COMMON_NONSENSE', 'AMR_COMMON_FRAMESHIFT',
              'AMR_COMMON_START_LOSS', 'AMR_COMMON_SPLICING_SITE', 'AMR_COMMON_LONG_DELETION', 'AMR_COMMON_SV_IN_SPLICING', 'AMR_COMMON_STOP_LOSS',
              'AMR_COMMON_DEL', 'AMR_COMMON_INS', 'AMR_COMMON_SNP', 'AMR_COMMON_MNP',
              'AFR_COMMON_SYNONYMOUS', 'AFR_COMMON_MISSENSE', 'AFR_COMMON_NONSENSE', 'AFR_COMMON_FRAMESHIFT',
              'AFR_COMMON_START_LOSS', 'AFR_COMMON_SPLICING_SITE', 'AFR_COMMON_LONG_DELETION', 'AFR_COMMON_SV_IN_SPLICING', 'AFR_COMMON_STOP_LOSS',
              'AFR_COMMON_DEL', 'AFR_COMMON_INS', 'AFR_COMMON_SNP', 'AFR_COMMON_MNP',
              'EUR_COMMON_SYNONYMOUS', 'EUR_COMMON_MISSENSE', 'EUR_COMMON_NONSENSE', 'EUR_COMMON_FRAMESHIFT',
              'EUR_COMMON_START_LOSS', 'EUR_COMMON_SPLICING_SITE', 'EUR_COMMON_LONG_DELETION', 'EUR_COMMON_SV_IN_SPLICING', 'EUR_COMMON_STOP_LOSS',
              'EUR_COMMON_DEL', 'EUR_COMMON_INS', 'EUR_COMMON_SNP', 'EUR_COMMON_MNP',
              'EAS_COMMON_SYNONYMOUS', 'EAS_COMMON_MISSENSE', 'EAS_COMMON_NONSENSE', 'EAS_COMMON_FRAMESHIFT',
              'EAS_COMMON_START_LOSS', 'EAS_COMMON_SPLICING_SITE', 'EAS_COMMON_LONG_DELETION', 'EAS_COMMON_SV_IN_SPLICING', 'EAS_COMMON_STOP_LOSS',
              'EAS_COMMON_DEL', 'EAS_COMMON_INS', 'EAS_COMMON_SNP', 'EAS_COMMON_MNP',
              'SAS_COMMON_SYNONYMOUS', 'SAS_COMMON_MISSENSE', 'SAS_COMMON_NONSENSE', 'SAS_COMMON_FRAMESHIFT',
              'SAS_COMMON_START_LOSS', 'SAS_COMMON_SPLICING_SITE', 'SAS_COMMON_LONG_DELETION', 'SAS_COMMON_SV_IN_SPLICING', 'SAS_COMMON_STOP_LOSS',
              'SAS_COMMON_DEL', 'SAS_COMMON_INS', 'SAS_COMMON_SNP', 'SAS_COMMON_MNP',
              )

HEADER = "GENE_NAME\tPATHOGENIC\tUNCERTAIN\tBENIGN\tAMR_COMMON\tAFR_COMMON\tEUR_COMMON\tEAS_COMMON\tSAS_COMMON\t" \
                "SYNONYMOUS\tMISSENSE\tNONSENSE\tFRAMESHIFT\t" \
                "START_LOSS\tSPLICING_SITE\tLONG_DELETION\tSV_IN_SPLICING\tSTOP_LOSS\t" \
                "DEL\tINS\tSNP\tMNP\t" \
                "PATHOGENIC_SYNONYMOUS\tPATHOGENIC_MISSENSE\tPATHOGENIC_NONSENSE\tPATHOGENIC_FRAMESHIFT\t" \
                "PATHOGENIC_START_LOSS\tPATHOGENIC_SPLICING_SITE\tPATHOGENIC_LONG_DELETION\tPATHOGENIC_SV_IN_SPLICING\tPATHOGENIC_STOP_LOSS\t" \
                "PATHOGENIC_DEL\tPATHOGENIC_INS\tPATHOGENIC_SNP\tPATHOGENIC_MNP\t" \
                "UNCERTAIN_SYNONYMOUS\tUNCERTAIN_MISSENSE\tUNCERTAIN_NONSENSE\tUNCERTAIN_FRAMESHIFT\tUNCERTAIN_START_LOSS\t" \
                "UNCERTAIN_SPLICING_SITE\tUNCERTAIN_LONG_DELETION\tUNCERTAIN_SV_IN_SPLICING\tUNCERTAIN_STOP_LOSS\t" \
                "UNCERTAIN_DEL\tUNCERTAIN_INS\tUNCERTAIN_SNP\tUNCERTAIN_MNP\t" \
                "BENIGN_SYNONYMOUS\tBENIGN_MISSENSE\tBENIGN_NONSENSE\tBENIGN_FRAMESHIFT\tBENIGN_START_LOSS\t" \
                "BENIGN_SPLICING_SITE\tBENIGN_LONG_DELETION\tBENIGN_SV_IN_SPLICING\tBENIGN_STOP_LOSS\t" \
                "BENIGN_DEL\tBENIGN_INS\tBENIGN_SNP\tBENIGN_MNP\t" \
                "AMR_COMMON_SYNONYMOUS\tAMR_COMMON_MISSENSE\tAMR_COMMON_NONSENSE\tAMR_COMMON_FRAMESHIFT\t" \
                "AMR_COMMON_START_LOSS\tAMR_COMMON_SPLICING_SITE\tAMR_COMMON_LONG_DELETION\tAMR_COMMON_SV_IN_SPLICING\tAMR_COMMON_STOP_LOSS\t" \
                "AMR_COMMON_DEL\tAMR_COMMON_INS\tAMR_COMMON_SNP\tAMR_COMMON_MNP\t" \
                "AFR_COMMON_SYNONYMOUS\tAFR_COMMON_MISSENSE\tAFR_COMMON_NONSENSE\tAFR_COMMON_FRAMESHIFT\t" \
                "AFR_COMMON_START_LOSS\tAFR_COMMON_SPLICING_SITE\tAFR_COMMON_LONG_DELETION\tAFR_COMMON_SV_IN_SPLICING\tAFR_COMMON_STOP_LOSS\t" \
                "AFR_COMMON_DEL\tAFR_COMMON_INS\tAFR_COMMON_SNP\tAFR_COMMON_MNP\t" \
                "EUR_COMMON_SYNONYMOUS\tEUR_COMMON_MISSENSE\tEUR_COMMON_NONSENSE\tEUR_COMMON_FRAMESHIFT\t" \
                "EUR_COMMON_START_LOSS\tEUR_COMMON_SPLICING_SITE\tEUR_COMMON_LONG_DELETION\tEUR_COMMON_SV_IN_SPLICING\tEUR_COMMON_STOP_LOSS\t" \
                "EUR_COMMON_DEL\tEUR_COMMON_INS\tEUR_COMMON_SNP\tEUR_COMMON_MNP\t" \
                "EAS_COMMON_SYNONYMOUS\tEAS_COMMON_MISSENSE\tEAS_COMMON_NONSENSE\tEAS_COMMON_FRAMESHIFT\t" \
                "EAS_COMMON_START_LOSS\tEAS_COMMON_SPLICING_SITE\tEAS_COMMON_LONG_DELETION\tEAS_COMMON_SV_IN_SPLICING\tEAS_COMMON_STOP_LOSS\t" \
                "EAS_COMMON_DEL\tEAS_COMMON_INS\tEAS_COMMON_SNP\tEAS_COMMON_MNP\t" \
                "SAS_COMMON_SYNONYMOUS\tSAS_COMMON_MISSENSE\tSAS_COMMON_NONSENSE\tSAS_COMMON_FRAMESHIFT\t" \
                "SAS_COMMON_START_LOSS\tSAS_COMMON_SPLICING_SITE\tSAS_COMMON_LONG_DELETION\tSAS_COMMON_SV_IN_SPLICING\tSAS_COMMON_STOP_LOSS\t" \
                "SAS_COMMON_DEL\tSAS_COMMON_INS\tSAS_COMMON_SNP\tSAS_COMMON_MNP"

OUTPUTHEADER = "GENE_NAME\tPATHOGENIC\tUNCERTAIN\tBENIGN\t" \
                "PATHOGENIC_SYNONYMOUS\tPATHOGENIC_MISSENSE\tPATHOGENIC_NONSENSE\tPATHOGENIC_FRAMESHIFT\t" \
                "PATHOGENIC_START_LOSS\tPATHOGENIC_SPLICING_SITE\tPATHOGENIC_LONG_DELETION\tPATHOGENIC_SV_IN_SPLICING\tPATHOGENIC_STOP_LOSS\t" \
                "PATHOGENIC_DEL\tPATHOGENIC_INS\tPATHOGENIC_SNP\tPATHOGENIC_MNP\t" \
                "BENIGN_SYNONYMOUS\tBENIGN_MISSENSE\tBENIGN_NONSENSE\tBENIGN_FRAMESHIFT\tBENIGN_START_LOSS\t" \
                "BENIGN_SPLICING_SITE\tBENIGN_LONG_DELETION\tBENIGN_SV_IN_SPLICING\tBENIGN_STOP_LOSS\t" \
                "BENIGN_DEL\tBENIGN_INS\tBENIGN_SNP\tBENIGN_MNP\t" \
                "EAS_COMMON\tEAS_COMMON_MISSENSE\n"

def AGCFilter(ifn, ofn):
    num_gene = 0
    cPathogenic = [0 for x in range(MAX_COUNT + 1)]
    cAMR_MISSENSE = [0 for x in range(MAX_COUNT + 1)]
    cAFR_MISSENSE = [0 for x in range(MAX_COUNT + 1)]
    cEUR_MISSENSE = [0 for x in range(MAX_COUNT + 1)]
    cEAS_MISSENSE = [0 for x in range(MAX_COUNT + 1)]
    cSAS_MISSENSE = [0 for x in range(MAX_COUNT + 1)]

    ofd = open(ofn, "w")
    ifd = open(ifn, "r")
    ofd.write(OUTPUTHEADER)

    loffd = open("/Users/chungtsai_su/data/AGCD/psv1.tsv", "w")
    loffd.write("PATHOGENIC_LOF\tUNCERTAIN_LOF\tBENIGN_LOF\tAMR_LOF\tAFR_LOF\tEUR_LOF\tEAS_LOF\tSAS_LOF\t%s\n" % HEADER)

    pp2fd = open("/Users/chungtsai_su/data/AGCD/pp2.tsv", "w")
    pp2fd.write("TOTAL_MISSENSE\tPATHOGENIC_MISSENSE\tBENIGN_MISSENSE\tCOMMON_MISSENSE\t%s\n" % HEADER)
    bp1fd = open("/Users/chungtsai_su/data/AGCD/bp1.tsv", "w")
    bp1fd.write("TOTAL_MISSENSE\tPATHOGENIC_MISSENSE\tBENIGN_MISSENSE\tCOMMON_MISSENSE\t%s\n" % HEADER)

    for line in ifd:
        if re.match("^#", line):
            continue

        line = line.strip()
        items = line.split('\t')
        count = int(items[COLUMN.PATHOGENIC]) + int(items[COLUMN.UNCERTAIN]) + int(items[COLUMN.BENIGN]) + int(
            items[COLUMN.EAS_COMMON]);
        if count < MIN_TOTAL:
            continue
        # if ((int(items[COLUMN.PATHOGENIC]) < MIN_PATHOGENIC) or (int(items[COLUMN.BENIGN]) < MIN_BENIGN)):
        #    continue

        ## LOF mechanism (PSV1)
        iPathogenic_LOF = int(items[COLUMN.PATHOGENIC_NONSENSE]) + int(items[COLUMN.PATHOGENIC_FRAMESHIFT]) + \
                          int(items[COLUMN.PATHOGENIC_START_LOSS]) + int(items[COLUMN.PATHOGENIC_SPLICING_SITE]) + \
                          int(items[COLUMN.PATHOGENIC_LONG_DELETION]) + int(items[COLUMN.PATHOGENIC_SV_IN_SPLICING])
        iUncertain_LOF  = int(items[COLUMN.UNCERTAIN_NONSENSE]) + int(items[COLUMN.UNCERTAIN_FRAMESHIFT]) + \
                          int(items[COLUMN.UNCERTAIN_START_LOSS]) + int(items[COLUMN.UNCERTAIN_SPLICING_SITE]) + \
                          int(items[COLUMN.UNCERTAIN_LONG_DELETION]) + int(items[COLUMN.UNCERTAIN_SV_IN_SPLICING])
        iBenign_LOF     = int(items[COLUMN.BENIGN_NONSENSE]) + int(items[COLUMN.BENIGN_FRAMESHIFT]) + \
                          int(items[COLUMN.BENIGN_START_LOSS]) + int(items[COLUMN.BENIGN_SPLICING_SITE]) + \
                          int(items[COLUMN.BENIGN_LONG_DELETION]) + int(items[COLUMN.BENIGN_SV_IN_SPLICING])
        iAMR_LOF        = int(items[COLUMN.AMR_COMMON_NONSENSE]) + int(items[COLUMN.AMR_COMMON_FRAMESHIFT]) + \
                          int(items[COLUMN.AMR_COMMON_START_LOSS]) + int(items[COLUMN.AMR_COMMON_SPLICING_SITE]) + \
                          int(items[COLUMN.AMR_COMMON_LONG_DELETION]) + int(items[COLUMN.AMR_COMMON_SV_IN_SPLICING])
        iAFR_LOF        = int(items[COLUMN.AFR_COMMON_NONSENSE]) + int(items[COLUMN.AFR_COMMON_FRAMESHIFT]) + \
                          int(items[COLUMN.AFR_COMMON_START_LOSS]) + int(items[COLUMN.AFR_COMMON_SPLICING_SITE]) + \
                          int(items[COLUMN.AFR_COMMON_LONG_DELETION]) + int(items[COLUMN.AFR_COMMON_SV_IN_SPLICING])
        iEUR_LOF        = int(items[COLUMN.EUR_COMMON_NONSENSE]) + int(items[COLUMN.EUR_COMMON_FRAMESHIFT]) + \
                          int(items[COLUMN.EUR_COMMON_START_LOSS]) + int(items[COLUMN.EUR_COMMON_SPLICING_SITE]) + \
                          int(items[COLUMN.EUR_COMMON_LONG_DELETION]) + int(items[COLUMN.EUR_COMMON_SV_IN_SPLICING])
        iEAS_LOF        = int(items[COLUMN.EAS_COMMON_NONSENSE]) + int(items[COLUMN.EAS_COMMON_FRAMESHIFT]) + \
                          int(items[COLUMN.EAS_COMMON_START_LOSS]) + int(items[COLUMN.EAS_COMMON_SPLICING_SITE]) + \
                          int(items[COLUMN.EAS_COMMON_LONG_DELETION]) + int(items[COLUMN.EAS_COMMON_SV_IN_SPLICING])
        iSAS_LOF        = int(items[COLUMN.SAS_COMMON_NONSENSE]) + int(items[COLUMN.SAS_COMMON_FRAMESHIFT]) + \
                          int(items[COLUMN.SAS_COMMON_START_LOSS]) + int(items[COLUMN.SAS_COMMON_SPLICING_SITE]) + \
                          int(items[COLUMN.SAS_COMMON_LONG_DELETION]) + int(items[COLUMN.SAS_COMMON_SV_IN_SPLICING])

        # if (iBenign_LOF + max(iAMR_LOF, iAFR_LOF, iEUR_LOF, iEAS_LOF, iSAS_LOF)) < MAX_PSV1_BENIGN and \
        #    (iPathogenic_LOF > MIN_PSV1_PATHO):
        if iPathogenic_LOF > 0:
            loffd.write("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n" % (iPathogenic_LOF, iUncertain_LOF, iBenign_LOF,
                                                                  iAMR_LOF, iAFR_LOF, iEUR_LOF, iEAS_LOF, iSAS_LOF, line))
        ##TODO: debug usage
        if items[0] == 'TSC2':
            print("%s\t%d\t%d\t%d\t%d\t%d\t%d" % (items[0], iBenign_LOF, iAMR_LOF, iAFR_LOF, iEUR_LOF, iEAS_LOF,
                                                  iSAS_LOF))

        ##Missense
        benign_missense = int(items[COLUMN.BENIGN_MISSENSE])
        common_missense = max(int(items[COLUMN.AMR_COMMON_MISSENSE]), int(items[COLUMN.AFR_COMMON_MISSENSE]), int(items[COLUMN.EUR_COMMON_MISSENSE]), int(items[COLUMN.EAS_COMMON_MISSENSE]), int(items[COLUMN.SAS_COMMON_MISSENSE]))
        total_missense = int(items[COLUMN.PATHOGENIC_MISSENSE]) + benign_missense + common_missense
        ## pathogenic benign (PP2)
        if (total_missense >= MIN_TOTAL_MISSENSE) and ((benign_missense+common_missense) <= MAX_PP2_THRESHOLD * total_missense):
            pp2fd.write("%d\t%d\t%d\t%d\t%s\n" % (total_missense, int(items[COLUMN.PATHOGENIC_MISSENSE]), benign_missense, common_missense, line))
        ## commom benign (BP1)
        if (total_missense >= MIN_TOTAL_MISSENSE) and ((benign_missense+common_missense) >= MIN_BP1_THRESHOLD * total_missense):
            bp1fd.write("%d\t%d\t%d\t%d\t%s\n" % (total_missense, int(items[COLUMN.PATHOGENIC_MISSENSE]), benign_missense, common_missense, line))

        ## Statistical profiling
        count = MAX_COUNT if (int(items[COLUMN.PATHOGENIC]) >= MAX_COUNT) else int(items[COLUMN.PATHOGENIC]);
        cPathogenic[count] += 1
        cAMR_MISSENSE[MAX_COUNT if (int(items[COLUMN.AMR_COMMON_MISSENSE]) >= MAX_COUNT) else int(
            items[COLUMN.AMR_COMMON_MISSENSE])] += 1
        cAFR_MISSENSE[MAX_COUNT if (int(items[COLUMN.AFR_COMMON_MISSENSE]) >= MAX_COUNT) else int(
            items[COLUMN.AFR_COMMON_MISSENSE])] += 1
        cEUR_MISSENSE[MAX_COUNT if (int(items[COLUMN.EUR_COMMON_MISSENSE]) >= MAX_COUNT) else int(
            items[COLUMN.EUR_COMMON_MISSENSE])] += 1
        cEAS_MISSENSE[MAX_COUNT if (int(items[COLUMN.EAS_COMMON_MISSENSE]) >= MAX_COUNT) else int(
            items[COLUMN.EAS_COMMON_MISSENSE])] += 1
        cSAS_MISSENSE[MAX_COUNT if (int(items[COLUMN.SAS_COMMON_MISSENSE]) >= MAX_COUNT) else int(
            items[COLUMN.SAS_COMMON_MISSENSE])] += 1

        #ofd.write("%s\n" % line)

        ofd.write("%s\t%d\t%d\t%d"
                  "\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s"
                  "\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s"
                  "\t%s\t%s\n" %
                  (items[COLUMN.GENE_NAME], int(items[COLUMN.PATHOGENIC]), int(items[COLUMN.UNCERTAIN]),int(items[COLUMN.BENIGN]),
                   items[COLUMN.PATHOGENIC_SYNONYMOUS], items[COLUMN.PATHOGENIC_MISSENSE], items[COLUMN.PATHOGENIC_NONSENSE],
                   items[COLUMN.PATHOGENIC_FRAMESHIFT], items[COLUMN.PATHOGENIC_START_LOSS], items[COLUMN.PATHOGENIC_SPLICING_SITE],
                   items[COLUMN.PATHOGENIC_LONG_DELETION], items[COLUMN.SV_IN_SPLICING], items[COLUMN.PATHOGENIC_STOP_LOSS], items[COLUMN.PATHOGENIC_DEL],
                   items[COLUMN.PATHOGENIC_INS], items[COLUMN.PATHOGENIC_SNP], items[COLUMN.PATHOGENIC_MNP],
                   items[COLUMN.BENIGN_SYNONYMOUS], items[COLUMN.BENIGN_MISSENSE], items[COLUMN.BENIGN_NONSENSE],
                   items[COLUMN.BENIGN_FRAMESHIFT], items[COLUMN.BENIGN_START_LOSS], items[COLUMN.BENIGN_SPLICING_SITE],
                   items[COLUMN.BENIGN_LONG_DELETION], items[COLUMN.BENIGN_SV_IN_SPLICING], items[COLUMN.BENIGN_STOP_LOSS], items[COLUMN.BENIGN_DEL],
                   items[COLUMN.BENIGN_INS], items[COLUMN.BENIGN_SNP], items[COLUMN.BENIGN_MNP],
                   items[COLUMN.EAS_COMMON], items[COLUMN.EAS_COMMON_MISSENSE]))
        num_gene += 1

    loffd.close()
    pp2fd.close()
    bp1fd.close()

    ## Profiling output
    tfd = open("/Users/chungtsai_su/data/AGCD/pathogenic-distribution.tsv", "w")
    for i in range(MAX_COUNT + 1):
        tfd.write("%d\t%d\n" % (i, cPathogenic[i]))
    tfd.close()
    tfd = open("/Users/chungtsai_su/data/AGCD/common-missense-distribution.tsv", "w")
    for i in range(MAX_COUNT + 1):
        tfd.write("%d\t%d\t%d\t%d\t%d\t%d\n" % (
        i, cAMR_MISSENSE[i], cAFR_MISSENSE[i], cEUR_MISSENSE[i], cEAS_MISSENSE[i], cSAS_MISSENSE[i]))
    tfd.close()

    ifd.close()
    ofd.close()
    print("Total %d genes" % num_gene)

    return


def main(argv):
    infile = ""
    outfile = ""
    global MIN_TOTAL, MIN_PATHOGENIC, MIN_BENIGN

    try:
        opts, args = getopt.getopt(argv, "hi:o:t:p:b:")
    except getopt.GetoptError:
        Usage()
        sys.exit(1)

    for opt, arg in opts:
        if opt == '-h':
            Usage()
            sys.exit()
        elif opt in "-i":
            infile = arg
            outfile = infile + ".out"
        elif opt in "-o":
            outfile = arg
        elif opt in "-t":
            MIN_TOTAL = int(arg)
        elif opt in "-p":
            MIN_PATHOGENIC = int(arg)
        elif opt in "-b":
            MIN_BENIGN = int(arg)
    if not infile:
        Usage()
        sys.exit(2)
    elif not os.path.isfile(infile):
        print("Error: input file(%s) is not existed" % (infile))
        Usage()
        sys.exit(3)

    # Main Function
    AGCFilter(infile, outfile)

    return


if __name__ == '__main__':
    main(sys.argv[1:])
