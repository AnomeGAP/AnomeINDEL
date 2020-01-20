import sys
import gzip
from collections import defaultdict

files = ["sample01", "sample02", "sample03", "sample04", "sample05", "sample06", "sample07", "sample08", "sample09", "sample10", "4T1",  "MC38", "CT26.WT", "A-8", "K14-L", "K14-H"]
#files = ["sample01", "sample02", "sample07", "sample09", "4T1",  "MC38", "CT26.WT", "A-8", "K14-L", "K14-H"]

MIN_SINGAL = 3
num_singals = MIN_SINGAL

h_gene = defaultdict(str)
h_tpm = defaultdict(str)
h_found = defaultdict(int)

for k in files:
    print("%s" % k)
    fd = open("./RNA-BAM/%s_genes.out" % k, "r")
    idx = 0
    h_found.clear()
    sum_tpm = 0.0
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
        idx += 1
    print("%s has %d genes and its TPM summary = %f" % (k,idx-1, sum_tpm))

out = open("./RNA-BAM/matrix.tsv", "w")
for i in range(1, idx):
    if h_tpm[h_gene[i]].find(".") >= 0:
        items = h_tpm[h_gene[i]].split("\t")
        # num_singals = 0
        # for j in items:
        #     if float(j) > 0:
        #         num_singals += 1
        if len(items) == len(files):
            if num_singals >= MIN_SINGAL:
                #print("%s" % h_tpm[h_gene[i]])
                out.write("%s\n" % h_tpm[h_gene[i]])
            else:
                print("Remove the weak signal - %s, %s" % (h_gene[i], h_tpm[h_gene[i]]))
        else:
            print("ERROR at %s(%d) has %d elements %s" % (h_gene[i], i, len(items), h_tpm[h_gene[i]]))
out.close()
