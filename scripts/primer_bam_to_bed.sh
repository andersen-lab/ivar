#!/bin/bash

# TODO uncomment if appropriate. See https://sipb.mit.edu/doc/safe-shell/
# set -euf -o pipefail

c=$1
p=$2
r=$3

cut -f 1,2 -d ',' "$c" | sed 's/^/>/g' | tr ',' '\n' > "$p".fa
~/Documents/bwa/bwa mem -k 5 -T 20 "$r" "$p".fa | samtools view -bS -F 4 -o "$p".bam # Tune -T parameter.
bedtools bamtobed -i "$p".bam > "$p".bed
