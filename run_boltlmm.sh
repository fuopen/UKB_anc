#!/bin/bash -x

BOLT= # path to the executable bolt software
chipdir= #directory at which the plink model builds files are
modelbuildprefix= #plink bed/fam/bim prefix
removelistfile= #sample ids to be removed
impdir= #directory where the (imputed) bgen files are,bgen files are named as "ukb_imp_chr2_v3.bgen" and sample files are named as "ukb27960_imp_chr2_v3_s487324.sample"
dataodir= # output directory for bolt-lmm
phenoPath= # path to the phenotype file
covarPath= # path to the covariate file
ldscorePath= # path to the LDSCORE files
genmapPath= # path to the genetic map files

if ! test -d $dataodir;then
  mkdir -p $dataodir
fi

run_bolt () {
chr=$1
pheno=$2

BOLT \
    --bfile=$chipdir/$modelbuildprefix \
    --remove=$chipdir/$removelistfile \
    --bgenFile=$impdir/ukb_imp_chr${chr}_v3.bgen \
    --sampleFile=$impdir/ukb27960_imp_chr22_v3_s487324.sample \
    --phenoFile=$phenoPath \
    --phenoCol=$pheno \
    --covarFile=$covarPath \
    --qCovarCol=age \
    --qCovarCol=age2 \
    --qCovarCol=agesex \
    --covarCol=sex \
    --covarCol=age2sex \
    --LDscoresFile=$ldscorePath \
    --geneticMapFile=$genmapPath \
    --lmmForceNonInf \
    --maxMissingPerSnp=1 \
    --numThreads=16 \
    --LDscoresMatchBp \
    --statsFile=$dataodir/bolt_${pheno}_chr${chr}.txt \
    --statsFileBgenSnps=$dataodir/bolt_imp_stats_${pheno}_chr${chr}.txt \
    --covarMaxLevels 300 \
    --verboseStats \
    --bgenMinMAF=5e-3 \
    --bgenMinINFO=0.3
#    --qCovarCol=PC{1:10} \
#    --qCovarCol=batch.{1:11} \
#    --qCovarCol=batch{1:94} \
}
for chr in {1..22};do
  for phe in birth_north dt;do
      run_bolt ${chr} $phe
  done
done
