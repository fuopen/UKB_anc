#!/bin/bash

BGENIEPATH=
bgenfile=
threads=
Snpfile=
phenofile=
covarfile=
outfile=

$BGENIEPATH \
--bgen "$bgenfile" \
--thread $threads \
--include_rsids $Snpfile \
--pheno "$phenofile" \
--covar "$covarfile" \
--pvals \
--out "$outfile"
