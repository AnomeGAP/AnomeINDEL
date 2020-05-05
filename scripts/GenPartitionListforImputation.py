#!/bin/env python3

import sys

#CONST
PARTITION_SIZE = 10000
num_partitions = 0
start_list = []
end_list = []
idx = 0
previous = ""
for line in sys.stdin:
    previous = line.strip()
    if idx % PARTITION_SIZE == 1:
        start_list.append(previous)
        num_partitions += 1
    idx += 1

buf = "START=("
for item in start_list:
    buf = buf + " " + item

buf += ")"
print(buf)

is_skip = True
buf = "END=("
for item in start_list:
    if is_skip:
        is_skip = False
    else:
        buf = buf + " " + str(int(item)-1)

buf = buf + " " + previous + ")"
print(buf)

print("Total: %s partitions" % num_partitions)
