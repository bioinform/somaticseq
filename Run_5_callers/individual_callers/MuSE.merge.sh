#!/bin/bash

cp 1.MuSE.txt ALL.MuSE.txt
cat 2.MuSE.txt 3.MuSE.txt 4.MuSE.txt 5.MuSE.txt 6.MuSE.txt 7.MuSE.txt 8.MuSE.txt 9.MuSE.txt 10.MuSE.txt 11.MuSE.txt 12.MuSE.txt 13.MuSE.txt 14.MuSE.txt 15.MuSE.txt 16.MuSE.txt 17.MuSE.txt 18.MuSE.txt 19.MuSE.txt 20.MuSE.txt 21.MuSE.txt 22.MuSE.txt X.MuSE.txt Y.MuSE.txt MT.MuSE.txt | egrep -v '^#' >> ALL.MuSE.txt
/home/ltfang/apps/downloads/MuSE/MuSEv1.0rc_submission_c039ffa sump -I ALL.MuSE.txt -E -O variants.E.vcf -D /home/ltfang/references/dbsnp/dbsnp144.GRCh37.vcf.gz
