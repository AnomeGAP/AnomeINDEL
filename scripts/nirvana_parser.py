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
# @file    nirvana_parser.py
#
# @brief   Parsing Nirvana output JSON file
#
# @author  Chung-Tsai Su(chungtsai_su@atgenomix.com)
#
# @date    2019/04/02
#
# @version 1.0
#
# @remark
#

import sys
import getopt
import json
from pprint import pprint
from collections import defaultdict

# CONSTANT
DAMAGE_LEVEL = {'stop_lost': 4,
                'start_lost': 4,
                'stop_gained': 4,
                'frameshift_variant': 4,
                'inframe_deletion': 3,
                'inframe_insertion': 3,
                'stop_retained_variant': 3,
                'start_retained_variant': 3,
                'incomplete_terminal_codon_variant': 3,
                'mature_miRNA_variant': 3,
                'splice_acceptor_variant': 3,
                'splice_donor_variant': 3,
                'coding_sequence_variant': 2,
                'missense_variant': 2,
                'splice_region_variant': 2,
                'upstream_gene_variant': 1,
                'downstream_gene_variant': 1,
                'synonymous_variant': 1,
                '5_prime_UTR_variant': 1,
                '3_prime_UTR_variant': 1,
                'intron_variant': 1,
                'non_coding_transcript_variant': 1,
                'non_coding_transcript_exon_variant': 1,
                'NMD_transcript_variant': 1}


# Mutation Type
h_mutation = defaultdict(lambda: 0)


def usage():
    print("zcat <Nirvana JSON Gzip file> | nirvana_parser.py -f <COSMIC Cancer Gene Census> ")
    print("Argument:")
    print("\t-h: Usage")
    print("\t-f: COSMIC Cancer Gene Census")
    print("Usage:")
    print("\tzcat ./KSTF_NIRVANA/KSTF_0601.json.gz | python ./nirvana_parser.py -f ./Census_all.tsv")
    print("\tgzcat /Users/chungtsai_su/KSTF/KSTF_NIRVANA/KSTF_0601.json.gz | python ./nirvana_parser.py -f /Users/chungtsai_su/data/COSMIC/Census_all.tsv")

    return


def check_tier(vars, gene_list):
    (tier, level) = (3, 0)
    gene_name = ""
    consequence = ""
    for j in range(len(vars['variants'])):
        # print("%d, %s, %s, %s, %s" % (j, vars['variants'][j]['begin'], vars['variants'][j]['end'],
        #                                          vars['variants'][j]['refAllele'],
        #                                          vars['variants'][j]['altAllele']))
        # print(json.dumps(vars['variants'][j], indent=4))
        # 'cosmic' / 'gene'
        type = "N/A"
        var = vars['variants'][j]
        if 'cosmic' in var:
            for k in range(len(var['cosmic'])):
                if 'gene' in var['cosmic'][k]:
                    temp = gene_list[var['cosmic'][k]['gene']]
                    if temp < tier:
                        tier = temp
                        gene_name = var['cosmic'][k]['gene']
        # 'transcripts' / 'ensembl' / 'hgnc'
        if 'transcripts' in var:
            # print(json.dumps(var['transcripts'], indent=4))
            if 'ensembl' in var['transcripts']:
                for k in range(len(var['transcripts']['ensembl'])):
                    if 'hgnc' in var['transcripts']['ensembl'][k]:
                        # print("*** %s ***" % var['transcripts']['ensembl'][k]['hgnc'])
                        temp = gene_list[var['transcripts']['ensembl'][k]['hgnc']]
                        if temp < tier:
                            tier = temp
                            gene_name = var['transcripts']['ensembl'][k]['hgnc']
                        # print("%s" % var['transcripts']['ensembl'][k]['consequence'])
                        for l in range(len(var['transcripts']['ensembl'][k]['consequence'])):
                            if level < DAMAGE_LEVEL[var['transcripts']['ensembl'][k]['consequence'][l]]:
                                level = DAMAGE_LEVEL[var['transcripts']['ensembl'][k]['consequence'][l]]
                                consequence = var['transcripts']['ensembl'][k]['consequence'][l]
                                # h_level[var['transcripts']['ensembl'][k]['consequence'][l]] += 1

    # if tier < 3 and level >= 3:
    #     print("%s:%s-%s\t%s\t%s\t%s\t%s" % (
    #         vars['variants'][j]['chromosome'], vars['variants'][j]['begin'], vars['variants'][j]['end'],
    #         vars['variants'][j]['refAllele'], vars['variants'][j]['altAllele'], gene_name, consequence))

    return tier, level


def nirvana_parser(ifile):
    h_gene = defaultdict(lambda: 3)
    tier_count = [0] * (3 + 1)
    multiallele_count = [0] * (3 + 1)
    level_count = [0] * (4 + 1)

    ifd = open(ifile, "r")
    # https://cancer.sanger.ac.uk/cosmic/census?tier=all
    # Gene Symbol,Name,Entrez GeneId,Genome Location,Tier,Hallmark,Chr Band,Somatic,Germline,Tumour Types(Somatic),
    # Tumour Types(Germline),Cancer Syndrome,Tissue Type,Molecular Genetics,Role in Cancer,Mutation Types,Translocation Partner,Other Germline Mut,Other Syndrome,Synonyms

    for line in ifd:
        items = line.strip().split("\t")
        # print("%s\t%s\t%s" % (items[0], items[4], items[14]))
        if items[4] != 'Tier':
            h_gene[items[0]] = int(items[4])
    ifd.close()

    data = json.load(sys.stdin)
    # pprint(data)
    pos = data['positions']
    # pprint('There are %d variants' % len(pos))
    for i in range(len(pos)):
        # print("i=%d" % i)
        # print(json.dumps(pos[i], indent=4))
        (tier, level) = check_tier(pos[i], h_gene)
        tier_count[tier] += 1
        if len(pos[i]['variants']) > 1:
            multiallele_count[tier] += 1
        level_count[level] += 1

    # print("Total\tTier1\tTier2\tNon-Oncogenes\tTier1[MultiAllele]\tTier2[MultiAllele]\tNon-Oncogenes[MultiAllele]\tLevel0\tLevel1\tLevel2\tLevel3\tLevel4")
    print("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d" %
          (len(pos), tier_count[1], tier_count[2], tier_count[3],
           multiallele_count[1], multiallele_count[2], multiallele_count[3],
           level_count[0], level_count[1], level_count[2], level_count[3], level_count[4]))

    # for k in h_level.keys():
    #     print("%s\t%d" % (k, h_level[k]))

    return


def main(argv):
    ifile = ""

    try:
        opts, args = getopt.getopt(argv, "hf:")
    except getopt.GetoptError:
        usage()
        sys.exit(1)

    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit()
        elif opt == "-f":
            ifile = arg

    nirvana_parser(ifile)

    return


if __name__ == '__main__':
    main(sys.argv[1:])
