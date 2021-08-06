#!/bin/bash


if [ $# -ne 5 ]
then
    echo "Please, give:"
    echo "1) input bam file"
    echo "2) bed file for introns, exons and tRNA"
    echo "3) stranded protocol (n/y)"
    echo "4) path to samtools"
    echo "5) path to bedtools"
    exit
fi

inbam=$1
refBED=$2
stranded=$3
p2samtools=$4
p2bedtools=$5
#refBED=/hpc/hub_oudenaarden/aalemany/vasaseq/ref_seqs/Mus_musculus.GRCm38.99.homemade_IntronExonTrna.bed

${p2samtools}/samtools view -h $inbam | awk 'BEGIN{OFS="\t"} {
    if ($1 ~ /^@/) {print $0} 
    else if ($0 ~ /NH:i:1\tHI:i:1\t/) {
        for (i=1; i<=NF; i++) {
            if ($i ~ /nM:i:[0-9]/) {
                col=i; nm=substr($col, 6, length($col))
            }
        };
        $1=$1";CG:"$6";nM:"nm; print $0
    }
}' | ${p2samtools}/samtools view -Sb > ${inbam%.bam}.singlemappers.bam

${p2bedtools}/bedtools bamtobed -i ${inbam%.bam}.singlemappers.bam | ${p2bedtools}/bedtools sort > ${inbam%.bam}.singlemappers.bed

if [ $stranded == "y" ]
then
    ${p2bedtools}/bedtools intersect -a ${inbam%.bam}.singlemappers.bed -b $refBED -wa -wb | awk 'BEGIN {OFS="\t"; w="T"} {
        chr=$1; readstart=$2; readend=$3; readname=$4; readstrand=$6; refstart=$8; refend=$9; refstrand=$10; refname=$11; genelen=$12; genestart=$13; geneend=$14;
        sx=match(readname, /;CG:/); rn=substr(readname, 0, sx-1); rq=substr(readname,sx+1,length(readname))
        if (readstrand==refstrand) {
            if ((readstart >= refstart) && (readend <= refend)) {
                readname=readname";jS:IN"; w="T"
            } else if ((readstart < refstart) && (readend > refend)) {
                readname=readname";jS:OUT"; w="F";
            } else if ( ((readstart < refstart)&&(readend <= refend)) || ((readstart <= refstart)&&(readend < refend)) ) {
                if (readstrand="+") {readname=readname";jS:5"} else {readname=readname";jS:3"}; w="T"
                } else if ( ((readstart > refstart) && (readend >= refend)) || ((readstart >= refstart) && (readend > refend)) ) {
                if (readstrand="+") {readname=readname";jS:3"} else {readname=readname";jS:5"}; w="T"
            } else {print $0 > "checkme.txt"}
            if (readstrand=="+") {x=1-(geneend-readend)/genelen} else {x=1-(readstart-genestart)/genelen}
            sx=match(readname, /;CG:/); rn=substr(readname, 0, sx-1); rq=substr(readname,sx+1,length(readname))
            if (w=="T") {print chr, readstart, readend, rn, readstrand, refname, rq, refend-refstart, x}
        }
    }' > ${inbam%.bam}.singlemappers_genes.bed

elif [ $stranded == 'n' ]
then
    ${p2bedtools}/bedtools intersect -a ${inbam%.bam}.singlemappers.bed -b $refBED -wa -wb | awk 'BEGIN {OFS="\t"; w="T"} {
        chr=$1; readstart=$2; readend=$3; readname=$4; readstrand=$6; refstart=$8; refend=$9; refstrand=$10; refname=$11; genelen=$12; genestart=$13; geneend=$14;
        sx=match(readname, /;CG:/); rn=substr(readname, 0, sx-1); rq=substr(readname,sx+1,length(readname))
        if ((readstart >= refstart) && (readend <= refend)) {
            readname=readname";jS:IN"; w="T"
        } else if ((readstart < refstart) && (readend > refend)) {
            readname=readname";jS:OUT"; w="F";
        } else if ( ((readstart < refstart)&&(readend <= refend)) || ((readstart <= refstart)&&(readend < refend)) ) {
            if (readstrand="+") {readname=readname";jS:5"} else {readname=readname";jS:3"}; w="T"
            } else if ( ((readstart > refstart) && (readend >= refend)) || ((readstart >= refstart) && (readend > refend)) ) {
            if (readstrand="+") {readname=readname";jS:3"} else {readname=readname";jS:5"}; w="T"
        } else {print $0 > "checkme.txt"}
        if (readstrand=="+") {x=1-(geneend-readend)/genelen} else {x=1-(readstart-genestart)/genelen}
        sx=match(readname, /;CG:/); rn=substr(readname, 0, sx-1); rq=substr(readname,sx+1,length(readname))
        if (w=="T") {print chr, readstart, readend, rn, readstrand, refname, rq, refend-refstart, x}
    }' > ${inbam%.bam}.singlemappers_genes.bed
fi

gzip ${inbam%.bam}.singlemappers_genes.bed

#sort -k4 ${inbam%.bam}.singlemappers_genes.bed > ${inbam%.bam}.nsorted.singlemappers_genes.bed # unnecessary

rm ${inbam%.bam}.singlemappers.bam ${inbam%.bam}.singlemappers.bed
