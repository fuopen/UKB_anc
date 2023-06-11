#!/bin/bash

Help () {
	echo "you need to provide the following arguments to run the EM based ancestry specific allele frequency estimation software"
	echo "-a, --ancestry"
	echo "ancestry matrix file, should be NxA where N is the number of samples and A is the number of ancestry columns,this option is mandatory!"

	echo "-g, --genotype"
	echo "genotype matrix file, should be NxM where N is the number of samples and M is the number of SNPs/markers, the order of sample in ancestry matrix and genotype matrix should match each other,this option is mandatory!"

	echo "-s, --snplist"
	echo "SNP annotation file, should be M rows and there must be a column called "SNP", which is the ID column with respect to the column of genotype matrix, this option is mandatory"

	echo "-O, --outdir"
	echo "Output directory. The software will generate two files (.freq.gz) and (.rds) files with prefix specified by -s (snplist) option"

	echo "-h, --help"
	echo "help information"

	echo "-v, --verbose"
	echo "verbose setting, enabled by default"
}

if test $# -eq 0;then
	echo "No arguments provided, quit!"
	Help
	exit 1
fi

sft=1
vb=
ancfile=
genofile=
snpfile=
outDir=
while ! test $# -eq 0;do
	case $1 in
	-h | --help)
		Help
		exit 0
	;;
	-a | --snplist)
		ancfile=$2
		if ! test -f $ancfile;then
			echo "error: ancestry file doesn't exist!"
			exit 1
		fi
		sft=2
	;;
	-g | --genotype)
		genofile=$2
		if ! test -f $genofile;then
			echo "error: genotype file doesn't exist!"
			exit 1
		fi
		sft=2
	;;
	-s | --snplist)
		snpfile=$2
		if ! test -f $snpfile;then
			echo "error: snp annotation file doesn't exist!"
			exit 1
		fi
		sft=2
	;;
	-O | --outdir)
		outDir=$2
		if ! test -d $outDir;then
			echo " the output directory doesn't exist, will try to create one (please bear in mind it might fail due to the file permission issue"
			mkdir $outDir
		fi
		sft=2
	;;
	-v | --verbose)
		vb=1
		sft=1
	;;
	*)
		echo "error:unrecongnised arguments,quit"
		Help
		exit 1
	;;
	esac
	shift $sft
done

if test -z $ancfile;then
	echo "error: ancestry file is not given!"
	exit 1
fi

if test -z $genofile;then
	echo "error: genotype file is not given!"
	exit 1
fi
if test -z $snpfile;then
	echo "error: snp annotation file is not given!"
	exit 1
fi

if test -z $outDir;then
	echo " the output directory is not given, will try to create one (please bear in mind it might fail due to the file permission issue) at your current directory!"
	outDir="$(pwd)/outAF"
	mkdir $outDir
fi

CURDIR=$(pwd)

if $vb -eq 1;then
	echo "verbose mode is turned on!"
	begin_time=$(date +%s)
	echo "Start to run the analysis!"
	
fi
Rscript "${CURDIR}/EM_allele_frequency_estimation.R" $ancfile $genofile $snpfile $outDir

if $vb -eq 1;then
	end_time=$(date +%s)
	echo "The analysis is completed in $((end_time - begin_time)) seconds"
fi
