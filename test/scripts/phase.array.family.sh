#!/bin/bash

../phase_common/bin/phase_common --input array/target.family.bcf --pedigree info/target.family.fam --region 1 --map info/chr1.gmap.gz --output target.phased.bcf --thread 8

