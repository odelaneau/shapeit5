#!/bin/bash

../phase_common/bin/phase_common --input array/target.haploid.bcf --haploids info/target.haploid.txt --region 1 --map info/chr1.gmap.gz --output target.phased.bcf --thread 8

