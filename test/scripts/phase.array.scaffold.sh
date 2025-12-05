#!/bin/bash

../phase_common/bin/phase_common --input array/target.unrelated.bcf --scaffold array/target.scaffold.bcf --region 1 --map info/chr1.gmap.gz --output target.phased.bcf --thread 8

