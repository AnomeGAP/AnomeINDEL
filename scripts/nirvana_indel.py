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
# @file    nirvana_indel.py
#
# @brief   Parsing Nirvana output JSON file and extract indel information
#
# @author  Chung-Tsai Su(chungtsai_su@atgenomix.com)
#
# @date    2019/06/19
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


# Mutation Type
h_mutation = defaultdict(lambda: 0)


def usage():
    print("zcat <Nirvana JSON Gzip file> | nirvana_indel.py -f <COSMIC Cancer Gene Census> -o <output file>")
    print("Argument:")
    print("\t-h: Usage")
    print("\t-f: COSMIC Cancer Gene Census")
    print("\t-f: Output file")
    print("Usage:")
    print("\tzcat ./KSTF_NIRVANA/KSTF_0601.json.gz | python ./nirvana_indel.py -f ./Census_all.tsv -o ./KSTF_NIRVANA/KSTF_0601.indel")
    print("\tgzcat /Users/chungtsai_su/KSTF/KSTF_NIRVANA/KSTF_0601.json.gz | python ./nirvana_indel.py -f /Users/chungtsai_su/data/COSMIC/Census_all.tsv -o /Users/chungtsai_su/KSTF/KSTF_NIRVANA/KSTF_0601.indel")

    return


def check_tier(vars, gene_list):
    is_indel = False
    is_insertion = False
    tier = 3
    refAllele = ""
    altAllele = ""
    gene_name = ""

    # print(json.dumps(vars, indent=4))

    for j in range(len(vars['variants'])):
        # print("%d, %s, %s, %s, %s" % (j, vars['variants'][j]['begin'], vars['variants'][j]['end'],
        #                                          vars['variants'][j]['refAllele'],
        #                                          vars['variants'][j]['altAllele']))
        # print(json.dumps(vars['variants'][j], indent=4))

        # 'cosmic' / 'gene'
        type = "N/A"
        var = vars['variants'][j]
        if len(var['altAllele']) == len(var['refAllele']) and var['altAllele'] != "-" and var['refAllele'] != "-":
            continue
        elif len(var['altAllele']) > len(var['refAllele']) or var['refAllele'] == '-':
            is_insertion = True
        else:
            is_insertion = False

        refAllele = var['refAllele']
        altAllele = var['altAllele']
        is_indel = True
        if 'cosmic' in var:
            for k in range(len(var['cosmic'])):
                if 'gene' in var['cosmic'][k]:
                    # print("*** %s ***" % var['cosmic'][k]['gene'])
                    temp = gene_list[var['cosmic'][k]['gene']]
                    if temp <= tier:
                        tier = temp
                        gene_name = var['cosmic'][k]['gene']
                        # print("tier=%d, genename = %s" % (tier, gene_name))
        # 'transcripts' / 'ensembl' / 'hgnc'
        if 'transcripts' in var:
            # print(json.dumps(var['transcripts'], indent=4))
            if 'ensembl' in var['transcripts']:
                for k in range(len(var['transcripts']['ensembl'])):
                    if 'hgnc' in var['transcripts']['ensembl'][k]:
                        # print("*** %s ***" % var['transcripts']['ensembl'][k]['hgnc'])
                        temp = gene_list[var['transcripts']['ensembl'][k]['hgnc']]
                        if temp <= tier:
                            tier = temp
                            gene_name = var['transcripts']['ensembl'][k]['hgnc']
                            # print("tier=%d, genename = %s" % (tier, gene_name))

    return is_indel, is_insertion, tier, refAllele, altAllele, gene_name


def nirvana_indel(ifile, ofile):
    h_gene = defaultdict(lambda: 3)
    tier_count = [0] * (3 + 1)
    num_insertion = num_deletion = 0

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
    ofn = open(ofile, "w")
    pos = data['positions']
    # pprint('There are %d variants' % len(pos))
    for i in range(len(pos)):
        # print("i=%d" % i)
        # print(json.dumps(pos[i], indent=4))
        len_ref = pos[i]["refAllele"]
        (is_indel, is_insertion, tier, refAllele, altAllele, gene_name) = check_tier(pos[i], h_gene)
        if is_indel:
            # print(gene_name)
            # print(json.dumps(pos[i], indent=4))
            tier_count[tier] += 1
            if is_insertion:
                num_insertion += 1
            else:
                num_deletion += 1
            # print("%s\t%s\t%s\t%s\t%s\t%d\n" % ( pos[i]["chromosome"], pos[i]["position"], refAllele, altAllele, gene_name, tier))
            ofn.write("%s\t%s\t%s\t%s\t%s\t%d\n" % ( pos[i]["chromosome"], pos[i]["position"], refAllele, altAllele, gene_name, tier))
    ofn.close()
    # print("Total\tIndel\tInsertion\tDeletion\tTeir 1\tTeir 2\n")
    print("%d\t%d\t%d\t%d\t%d\t%d\n" % (len(pos), tier_count[3]+tier_count[1]+tier_count[2], num_insertion, num_deletion, tier_count[1], tier_count[2]))
    return


def main(argv):
    ifile = ""
    ofile = ""

    try:
        opts, args = getopt.getopt(argv, "hf:o:")
    except getopt.GetoptError:
        usage()
        sys.exit(1)

    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit()
        elif opt == "-f":
            ifile = arg
        elif opt == "-o":
            ofile = arg

    nirvana_indel(ifile, ofile)

    return


if __name__ == '__main__':
    main(sys.argv[1:])
