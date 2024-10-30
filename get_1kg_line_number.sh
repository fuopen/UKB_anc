#!/bin/bash

chr=$1

if test $chr -ne 21 && test $chr -ne 22;then
    grep -xnf tmp_dir/chr${chr}_match.txt 1kg_snp_list/chr${chr}_1kg.txt |cut -d ':' -f1|sort -n > tmp_1kg_hap/chr${chr}_line_num.txt
else
    grep -xnf <(cat tmp_dir/chr${chr}_match.txt tmp_dir/chr${chr}_flip_match.txt) 1kg_snp_list/chr${chr}_1kg.txt |cut -d ':' -f1|sort -n >tmp_1kg_hap/chr${chr}_line_num.txt
fi
