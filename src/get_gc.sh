#!/bin/bash
export PATH=/programs/bioawk:$PATH

echo -e "filename\tgc_content" > gc_content.tsv
for f in /workdir/hdd29/theRefseqening/*.fa*; do
    avg_gc=$(bioawk -c fastx '{gc+=gsub(/[GCgc]/,"",$seq); len+=length($seq)} END {if (len>0) print (gc/len)*100; else print 0}' "$f")
    echo -e "$(basename "$f")\t$avg_gc" >> gc_content.tsv
done
