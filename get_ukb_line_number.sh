#!/bin/bash

chr=$1

if test $chr -ne 21 && test $chr -ne 22;then
    grep -xnf tmp_dir/chr${chr}_match.txt ukb_snp_list/chr${chr}_ukb.txt |cut -d ':' -f1|sort -n > tmp_ukb_hap/chr${chr}_line_num.txt
else
    grep -xnf tmp_dir/chr${chr}_match.txt  ukb_snp_list/chr${chr}_ukb.txt |cut -d ':' -f1|sort -n >tmp_ukb_hap/chr${chr}_line_num_match.txt
    grep -xnf <(awk 'BEGIN{OFS="\t"}{$1=$1;print $1,$2,$3,$5,$4}' tmp_dir/chr${chr}_flip_match.txt) ukb_snp_list/chr${chr}_ukb.txt |cut -d ':' -f1|sort -n >tmp_ukb_hap/chr${chr}_line_num_flip_match.txt
fi
