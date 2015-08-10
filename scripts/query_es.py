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
# @file    query_es.py
#
# @brief   An example to query Anome DB via Elasticsearch
#
# @author  Chung-Tsai Su(chungtsai_su@anome.com)
#
# @date    2015/08/11
#
# @version 1.0
#
# @remark
#

__author__ = 'chungtsai_su'

import sys
import getopt
from elasticsearch import Elasticsearch

def Usage():
    ''' Description: Program Usage
        Argument:    NONE
        Return:	     NONE
    '''
    print("query_es.py -q <RS id> ")
    print("Argument:")
    print("\t-h: Usage")
    print("\t-q: dbSNP id")
    print("Usage:")
    print("\t./query_es.py -q 62621221")

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

    es = Elasticsearch(hosts=[{"host": "140.112.183.200","port": 9200}])

    ##TODO: add your query here
    result = es.search(index="anome_base",doc_type="anome_table",
                        body={"query": {"match": {"dbsnp_rs":rs}}})

    print("Query Result:\n")
    print("%s" % (result))
    
    return

if __name__ == '__main__':
    main(sys.argv[1:])
