#!/bin/env python
#
# @note Copyright (C) 2015, Anome Incorporated. All Rights Reserved.
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
# @file    es_query.py
#
# @brief   An example to query Anome DB by chr_start_end_ref_alt via Elasticsearch
#
# @author  Chung-Tsai Su(chungtsai_su@anome.com)
#
# @date    2015/08/11
#
# @version 1.0
#
# @remark
#

import sys
import getopt
from elasticsearch import Elasticsearch

ES_SERVER_IP="140.112.183.200"
ES_SERVER_PORT=9200

def Usage():
    ''' Description: Program Usage
        Argument:    NONE
        Return:	     NONE
    '''
    print("es_search.py -i <CHROM_START_END_REF_ALT> ")
    print("Argument:")
    print("\t-h: Usage")
    print("\t-i: <CHROM_START_END_REF_ALT>")
    print("Usage:")
    print("\t./es_search.py -i 9_135804266_135804266_G_A")

    return

def main(argv):
    key = ""

    try:
        opts, args = getopt.getopt(argv,"hi:")
    except getopt.GetoptError:
        Usage()
        sys.exit(1)

    for opt, arg in opts:
        if opt == '-h':
            Usage()
            sys.exit()
        elif opt in ("-i"):
            key = arg
    if not key:
        Usage()
        sys.exit(2)

    es = Elasticsearch(hosts=[{"host": ES_SERVER_IP,"port": ES_SERVER_PORT}])

    ##TODO: add your query here
    result = es.get(index="anome_base", doc_type='anome_table', id=key)
    #print(result['_source'])

    print("_id=%(_id)s" % result)
    print("\tdbSNP id: %(dbsnp_rs)s" % result["_source"])
    print("\tid: %(id)s" % result["_source"])
    if "1000genome_anome_tot_maf" in result["_source"]:
        print("\t1000genome_anome_tot_maf=%(1000genome_anome_tot_maf).4f" % result["_source"])
        print("\t1000genome_anome_afr_maf=%(1000genome_anome_afr_maf).4f" % result["_source"])
        print("\t1000genome_anome_amr_maf=%(1000genome_anome_amr_maf).4f" % result["_source"])
        print("\t1000genome_anome_eur_maf=%(1000genome_anome_eur_maf).4f" % result["_source"])
        print("\t1000genome_anome_eas_maf=%(1000genome_anome_eas_maf).4f" % result["_source"])
        print("\t1000genome_anome_sas_maf=%(1000genome_anome_sas_maf).4f" % result["_source"])
    else:
         print("\tnot in the 1000 Genome Project")

    return

if __name__ == '__main__':
    main(sys.argv[1:])
