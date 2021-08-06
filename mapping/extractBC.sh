#!/bin/bash

if [ $# -ne 4 ]
then
    echo "Please, give: "
    echo "1) input root to fastq files"
    echo "2) protocol [celseq1, celseq2, vasaplate, vasadrop, 10x]"
    echo "3) path to concatenator.py"
    echo "4) output folder for fastq files"
    exit
fi

outfq=$1
protocol=$2
path2scripts=$3
outdir=$4
bcfile=$5 # ${path2scripts}/bc_celseq1.tsv, ${path2scripts}/bc_celseq2.tsv

mkdir -p $outdir

if [ $protocol == 'celseq1' ]
then
    python3 ${path2scripts}/concatenator.py --fqf ${outfq} --cbcfile ${path2scripts}/bc_celseq1.tsv --cbchd 0 --lenumi 4 --outdir ${outdir}
elif [ $protocol == 'celseq2' ]
then
    python3 ${path2scripts}/concatenator.py --fqf ${outfq} --cbcfile ${path2scripts}/bc_celseq2.tsv --cbchd 0 --lenumi 6 --umifirst --outdir ${outdir}
elif [ $protocol == 'vasaplate' ]
then
    python3 ${path2scripts}/concatenator.py --fqf ${outfq} --cbcfile ${path2scripts}/bc_celseq2.tsv --cbchd 0 --lenumi 6 --umifirst --demux --outdir ${outdir}
elif [ $protocol == 'vasadrop' ]
then
    python3 ${path2scripts}/concatenator.py --fqf ${outfq} --cbcfile $bcfile --cbchd 0 --lencbc 16 --lenumi 6 --umifirst --demux --outdir ${outdir}
elif [ $protocol == '10x' ]
then
    python3 ${path2scripts}/concatenator.py --fqf ${outfq} --cbcfile $bcfile --cbchd 0 --lencbc 16 --lenumi 12 --demux --outdir ${outdir}
else
    echo "unknown protocol [celseq1, celseq2, vasaplate, vasadrop]"
    exit
fi
