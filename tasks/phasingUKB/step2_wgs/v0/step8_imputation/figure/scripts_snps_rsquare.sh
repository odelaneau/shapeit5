#!/bin/bash

for MET in beagle5.4  shapeit5; do
	zcat /home/srubinac/fuse/PhasingUKB/Phasing/PhasingWGS/step8_imputation/concordance_test/imputed_1k_rp_${MET}_chr20.error.grp.txt.gz | grep "^GCsSAF" | awk '{print $2,$3,$4,$19,$19}' | gzip > /home/srubinac/Dropbox/workspace/shapeit5/tasks/phasingUKB/step2_wgs/step8_imputation/figure/data/imputed_1k_rp_${MET}_chr20.rsquare.grp.txt.gz
done
