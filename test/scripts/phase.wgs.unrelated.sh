#!/bin/bash

CRD=info/chunks.coordinates.txt
LST=tmp/chunks.files.txt
rm $LST

#PHASE SCAFFOLD SITES WITH MAF >=0.1%
../phase_common/bin/phase_common --input wgs/target.unrelated.bcf --filter-maf 0.001 --region 1 --map info/chr1.gmap.gz --output tmp/target.scaffold.bcf --thread 8

#PHASE RARE SITES ONTO SCAFFOLD
while read LINE; do
	CHK=$(echo $LINE | awk '{ print $1; }')
	SRG=$(echo $LINE | awk '{ print $3; }')
	IRG=$(echo $LINE | awk '{ print $4; }')

	OUT=tmp/target.phased.chunk$CHK\.bcf
	../phase_rare/bin/phase_rare --input wgs/target.unrelated.bcf --scaffold tmp/target.scaffold.bcf --map info/chr1.gmap.gz --input-region $IRG --scaffold-region $SRG --output $OUT  --thread 8

	echo $OUT >> $LST
done < $CRD
bcftools concat -n -Ob -o target.phased.bcf -f $LST
bcftools index target.phased.bcf


