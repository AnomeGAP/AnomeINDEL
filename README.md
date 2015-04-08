# AnomeINDEL
  There are some useful scripts for manipulating BAM files.

## Requirement

- [Pysam](https://github.com/pysam-developers/pysam): Python module for BAM/SAM parser

## Scripts

```
ExtractUnMappedReads.py -i <Input BAM> -o <Output SAM>
Argument:
	-h: Usage
	-i: Input BAM file
	-o: Output SAM file which contains all unmapped reads[Default=<inputfile>.unmapped.sam]
Usage:
	./ExtractUnMappedReads.py -i ../example/example.bam -o ../example/example.unmapped.sam
```

```
FindLongDeletion.py -i <Input BAM> -o <Output SAM> -m <Minimal Segment Length> -M <Maximal Segment Length> -q <Minimal Quality> -c <Minimal Coverage> -s <Minimal Number of Supports>
Argument:
	-h: Usage
	-i: Input BAM file
	-o: Output SAM file which contains all read pairs within specific segment length [Default=<inputfile>.ld.sam]
	-m: Minimal segment length [Default=500]
	-M: Maximal segment length [Default=50000]
	-q: Minimal read quality [Default=10]
	-c: Minimal ratio of overlapping reads [Default=0.900000]
	-s: Minimal number of supports [Default=5]
Usage:
	./FindLongDeletion.py -i ../example/example.bam -o ../example/example.ld.sam -m 500 -M 50000 -q 10 -c 0.9 -s 5
```

See also https://github.com/AnomeGAP/AnomeINDEL

last modified at 2015/04/08

