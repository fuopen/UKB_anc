#!/bin/bash

batch_dir=$1

chr=$2

if test -d $batch_dir;then
    cd $batch_dir
    cmd="perl /data/muscovy/not-backed-up/shu/hapmix/bin/runHapmix.pl chr${chr}.par"
    $cmd
fi
