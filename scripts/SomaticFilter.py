
import sys
import gzip
from collections import defaultdict

# CONSTANT
MIN_DP = 0
TEMPLATE="Y%02d_Y01.BALB_cJ.contam.filtered2.d100.vcf.gz"
MIN_SUPPORT = 3


def somatic():
    h_primary = defaultdict(int)
    h_ctc = defaultdict(int)
    h_secondary = defaultdict(int)
    num_primary = 0
    num_ctc = 0
    num_secondary = 0

    for i in range(2, 5):
        #print("%d" % i)
        ifile = TEMPLATE % i
        ifd = gzip.open(ifile, "rt")

        for line in ifd:
            if not line.startswith("#"):
                # 1       4351880 .       TA      T       .       artifact_in_normal;str_contraction      CONTQ=93;DP=141;ECNT=1;GERMQ=265;MBQ=30,30;MFRL=257,277;MMQ=60,60;MPOS=20;NALOD=-1.178e+00;NLOD=11.66;POPAF=6.00;RPA=11,10;RU=A;SAAF=0.051,0.051,0.065;SAPP=5.753e-03,0.021,0.973;STR;TLOD=5.48  GT:AD:AF:DP:F1R2:F2R1:OBAM:OBAMRC       0/0:56,2:0.049:58:28,1:28,1:false:false 0/1:72,5:0.075:77:25,3:47,2:false:false
                items = line.strip().split("\t")
                if items[6] == "PASS":
                    variant = "%s:%s:%s:%s" % (items[0], items[1], items[3], items[4])
                    h_primary[variant] += 1
                    #print("%s" % line)
                    if h_primary[variant] == 1:
                        num_primary += 1
    num_common = 0
    for k in h_primary:
        if h_primary[k] >= MIN_SUPPORT:
            print("%s" % k)
            num_common += 1
    print ("primary has %d common variants" % num_common)
    print("%d variants are primary" % num_primary)

    for i in range(8, 11):
        #print("%d" % i)
        ifile = TEMPLATE % i
        ifd = gzip.open(ifile, "rt")

        for line in ifd:
            if not line.startswith("#"):
                # 1       4351880 .       TA      T       .       artifact_in_normal;str_contraction      CONTQ=93;DP=141;ECNT=1;GERMQ=265;MBQ=30,30;MFRL=257,277;MMQ=60,60;MPOS=20;NALOD=-1.178e+00;NLOD=11.66;POPAF=6.00;RPA=11,10;RU=A;SAAF=0.051,0.051,0.065;SAPP=5.753e-03,0.021,0.973;STR;TLOD=5.48  GT:AD:AF:DP:F1R2:F2R1:OBAM:OBAMRC       0/0:56,2:0.049:58:28,1:28,1:false:false 0/1:72,5:0.075:77:25,3:47,2:false:false
                items = line.strip().split("\t")
                if items[6] == "PASS":
                    dp = int(items[10].split(":")[3])
                    if dp >= MIN_DP:
                        variant = "%s:%s:%s:%s" % (items[0], items[1], items[3], items[4])
                        h_ctc[variant] += 1
                        #print("%s" % line)
                        if h_ctc[variant] == 1:
                            num_ctc += 1
    num_common = 0
    for k in h_ctc:
        if h_ctc[k] >= MIN_SUPPORT:
            print("%s" % k)
            num_common += 1
    print("ctc has %d common variants" % num_common)

    num_pass = 0
    for k in h_ctc:
        if h_ctc[k] >= MIN_SUPPORT and h_primary[k] < 1:
            print("%s" % k)
            num_pass += 1
    print("%d/%d qualified variants in ctc group are not found in primary" % (num_pass, num_ctc))

    for i in range(5, 8):
        #print("%d" % i)
        ifile = TEMPLATE % i
        ifd = gzip.open(ifile, "rt")

        for line in ifd:
            if not line.startswith("#"):
                # 1       4351880 .       TA      T       .       artifact_in_normal;str_contraction      CONTQ=93;DP=141;ECNT=1;GERMQ=265;MBQ=30,30;MFRL=257,277;MMQ=60,60;MPOS=20;NALOD=-1.178e+00;NLOD=11.66;POPAF=6.00;RPA=11,10;RU=A;SAAF=0.051,0.051,0.065;SAPP=5.753e-03,0.021,0.973;STR;TLOD=5.48  GT:AD:AF:DP:F1R2:F2R1:OBAM:OBAMRC       0/0:56,2:0.049:58:28,1:28,1:false:false 0/1:72,5:0.075:77:25,3:47,2:false:false
                items = line.strip().split("\t")
                if items[6] == "PASS":
                    dp = int(items[10].split(":")[3])
                    if dp >= MIN_DP:
                        variant = "%s:%s:%s>%s" % (items[0], items[1], items[3], items[4])
                        h_secondary[variant] += 1
                        if h_secondary[variant] == 1:
                            num_secondary += 1
    num_common = 0
    for k in h_ctc:
        if h_ctc[k] >= MIN_SUPPORT:
            print("%s" % k)
            num_common += 1
    print("secondary has %d common variants" % num_common)

    num_pass = 0
    for k in h_secondary:
        if h_secondary[k] >= MIN_SUPPORT and h_ctc[k] >= MIN_SUPPORT and h_primary[k] < 1:
            print("%s" % k)
            num_pass += 1
    print("%d/%d qualified variants in Secondary group are found in ctc" % (num_pass, num_secondary))

    num_pass = 0
    for k in h_secondary:
        if h_secondary[k] >= MIN_SUPPORT and h_primary[k] < 1 and h_ctc[k] < 1:
            print("%s" % k)
            num_pass += 1
    print("%d/%d qualified variants in Secondary group are not found in primary and ctc" % (num_pass, num_secondary))

    return


def main(argv):

    somatic()

    return


if __name__ == '__main__':
    main(sys.argv[1:])
