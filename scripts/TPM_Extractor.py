import sys
import gzip
from collections import defaultdict

#INPUT_FILE="./RNA-BAM/%s_genes.out"
#OUTPUT_FILE="./RNA-BAM/matrix.tsv"
#files = ["sample01", "sample02", "sample03", "sample04", "sample05", "sample06", "sample07", "sample08", "sample09", "sample10"]
##files = ["sample01", "sample02", "sample03", "sample04", "sample05", "sample06", "sample07", "sample08", "sample09", "sample10", "4T1",  "MC38", "CT26.WT", "A-8", "K14-L", "K14-H"]
##files = ["sample01", "sample02", "sample07", "sample09", "4T1",  "MC38", "CT26.WT", "A-8", "K14-L", "K14-H"]

INPUT_FILE="/Users/chungtsai_su/data/GA816_RNA-Seq/7.TPM/%s_genes.out"
OUTPUT_FILE="/Users/chungtsai_su/data/GA816_RNA-Seq/7.TPM/matrix.tsv"
files = ["YCC1", "YCC2", "YCC3", "YCC4", "YCC5", "YCC6", "GA816-NC_S7"]
#files = ["YCC1", "YCC2", "YCC3", "YCC4", "GA816-NC_S7"]

MIN_SINGAL = 2
num_singals = MIN_SINGAL

h_gene = defaultdict(str)
h_tpm = defaultdict(str)
h_found = defaultdict(int)

for k in files:
    #print("%s" % k)
    fd = open(INPUT_FILE % k, "r")
    idx = 0
    h_found.clear()
    sum_tpm = 0.0
    cnt = 0
    for line in fd:
        items = line.strip().split("\t")
        if idx > 0:
            if not items[0] in h_found:
                h_found[items[0]] = 1
                # print("%s\t%s" % (items[0], items[6]))
                if h_gene[idx] == "":
                    h_gene[idx] = items[0]
                    h_tpm[items[0]] = items[6]
                elif h_gene[idx] == items[0]:
                    h_tpm[items[0]] = "%s\t%s" % (h_tpm[items[0]], items[6])
                else:
                    print("ERROR")
                sum_tpm += float(items[6])
                if float(items[6]) > 0:
                    cnt += 1
        idx += 1
    # print("%s has %d genes and its TPM summary = %f" % (k,idx-1, sum_tpm))
    print("%s has %d/%d genes with non-zero TPM" % (k, cnt, idx-1))

out = open(OUTPUT_FILE, "w")
for i in range(1, idx):
    if h_tpm[h_gene[i]].find(".") >= 0:
        items = h_tpm[h_gene[i]].split("\t")
        #num_singals = 0
        #for j in items:
        #    if float(j) > 0:
        #        num_singals += 1
        if len(items) == len(files):
            if num_singals >= MIN_SINGAL:
                #print("%s" % h_tpm[h_gene[i]])
                out.write("%s\n" % h_tpm[h_gene[i]])
                #out.write("%s\t%s\n" % (h_gene[i],h_tpm[h_gene[i]]))
            else:
                print("Remove the weak signal - %s, %s" % (h_gene[i], h_tpm[h_gene[i]]))
        else:
            print("ERROR at %s(%d) has %d elements %s" % (h_gene[i], i, len(items), h_tpm[h_gene[i]]))
out.close()