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
# @file    es_search.py
#
# @brief   An example to query Anome DB by dbSNP RS ID via Elasticsearch
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
    print("es_search.py -q <RS id> ")
    print("Argument:")
    print("\t-h: Usage")
    print("\t-q: dbSNP id")
    print("Usage:")
    print("\t./es_search.py -q 62621221")

    return

def main(argv):
    rs = ""

    try:
        opts, args = getopt.getopt(argv,"hq:")
    except getopt.GetoptError:
        Usage()
        sys.exit(1)

    for opt, arg in opts:
        if opt == '-h':
            Usage()
            sys.exit()
        elif opt in ("-q"):
            rs = arg
    if not rs:
        Usage()
        sys.exit(2)

    es = Elasticsearch(hosts=[{"host": ES_SERVER_IP,"port": ES_SERVER_PORT}])

    ##TODO: add your query here
    result = es.search(index="anome_base",doc_type="anome_table",
                        body={"query": {"match": {"dbsnp_rs":rs}}})

    #print("Query Result:\n")
    #print("%s\n" % (result))

    print("There are %d documents matched by the query" % (result["hits"]["total"]))

    for item in result["hits"]["hits"]:
        print("\t_id=%(_id)s" % item)
        print("\tdbSNP id: %(dbsnp_rs)s" % item["_source"])
        print("\tid: %(id)s" % item["_source"])
        if "1000genome_anome_tot_maf" in item["_source"]:
            print("\t1000genome_anome_tot_maf=%(1000genome_anome_tot_maf).4f" % item["_source"])
            print("\t1000genome_anome_afr_maf=%(1000genome_anome_afr_maf).4f" % item["_source"])
            print("\t1000genome_anome_amr_maf=%(1000genome_anome_amr_maf).4f" % item["_source"])
            print("\t1000genome_anome_eur_maf=%(1000genome_anome_eur_maf).4f" % item["_source"])
            print("\t1000genome_anome_eas_maf=%(1000genome_anome_eas_maf).4f" % item["_source"])
            print("\t1000genome_anome_sas_maf=%(1000genome_anome_sas_maf).4f" % item["_source"])
        else:
            print("\tnot in the 1000 Genome Project")
        print("")

    return

if __name__ == '__main__':
    main(sys.argv[1:])
