#!/bin/bash

for MET in beagle5.4  shapeit5; do
	zcat /home/srubinac/fuse/PhasingUKB/UKB_PHASING_EXOME_ARRAY/step6_imputation/imputation/concordance_test/imputed_1k_rp_${MET}_allchrs.error.grp.txt.gz | grep "^GCsSAF" | awk '{print $2,$3,$4,$19,$19}' | gzip > /home/srubinac/Dropbox/workspace/shapeit5/tasks/phasingUKB/step3_wes/step6_imputation/figure/data/imputed_1k_rp_${MET}_allchrs.rsquare.grp.txt.gz
done
