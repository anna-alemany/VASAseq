#!/bin/bash

### input paths (to modify by user)
p2s=/exports/ana-scarlab/aalemany/bin/vasaplate_split   # path to mapping scripts in your computer/HPC
p2trimgalore=/exports/ana-scarlab/bin/TrimGalore-0.6.6  # path to TrimGalore
p2cutadapt=/home/aalemany/anaconda3/bin/                # path to cutadapt
p2bwa=/exports/ana-scarlab/bin/bwa-0.7.17               # path to BWA
p2samtools=/exports/ana-scarlab/bin/samtools-1.11       #  path to samtools
p2star=/exports/ana-scarlab/bin/STAR-2.7.7a/bin/Linux_x86_64    # path to STAR
p2bedtools=/exports/ana-scarlab/bin/bedtools2/bin          # path to bedtools
email=a.alemany@lumc.nl                                 # email

### check input parameters
if [ $# -ne 7 ]
then
    echo "Please, give:"
    echo "1) library name (prefix of the fastq files, name before _R1.fastq.gz and _R2.fastq.gz)"
    echo "2) genome: MOUSE /  HUMAN / MIXED "
    echo "3) read length (for MOUSE: 59, 74, 96, 246; for HUMAN: 74, 136; for MIXED: 74, 91, 135)"
    echo "4) prefix for output files"
    echo "5) folder for output files"
    echo "6) fastqfile extraction (y/n)"
    echo "7) cellID from filename or readname (f/r)"
    exit
fi

lib=$1
ref=$2
n=$3
out=$4
folder=$5
fqext=$6
cellidori=$7

### check existence of input fastq files
r1=$(ls ${lib}*_R1*.fastq.gz)
r2=$(ls ${lib}*_R2*.fastq.gz)
echo $r1 $r2
if [ ${#r1} == 0 ]
then
    echo "R1 fastq files not found"
    exit
fi
if [ ${#r2} == 0 ]
then
    echo "R2 fastq files not found"
    exit
fi

### check python version (we want version 3)
v=$(python -c 'import sys; print(".".join(map(str, sys.version_info[:3])))' | awk -F "." '{print $1}')
if [ $v -ne "3" ]
then
    echo "python needs to be 3"
    exit
fi

### set references
if [[ $ref == "MOUSE" ]]
then
    riboref=/hpc/hub_oudenaarden/aalemany/vasaseq/ref_seqs/rRNA_mouse_Rn45S.fa
    genome=/hpc/hub_oudenaarden/group_references/ensembl/99/mus_musculus/star_v273a_NOMASK_NOERCC_index_$n
    refBED=/hpc/hub_oudenaarden/aalemany/vasaseq/ref_seqs/Mus_musculus.GRCm38.99.homemade_IntronExonTrna.bed
elif [[ $ref == "HUMAN" ]]
then
    riboref=/exports/ana-scarlab/group_references/ensembl/human/99/unique_rRNA_human.fa
    genome=/exports/ana-scarlab/group_references/ensembl/human/99/star_v277a_index_$n
    refBED=/exports/ana-scarlab/group_references/ensembl/human/99/Homosapines_ensemble99.homemade_IntronExonTrna.bed
fi 

if [ ! -d $genome ]
then
    echo "genome not found"
    exit
fi

### extract cell barcodes
jcbc=1
if [ $fqext == "y" ]
then
    jobid=extract-${lib}
    jcbc=$(sbatch --export=All -c 1 -N 1 -J ${jobid} -e ${jobid}.err -o ${jobid}.out -t 48:00:00 --mem=10G --mail-type=END --mail-user=${email} --wrap="${p2s}/extractBC.sh ${lib} vasaplate ${p2s} ${folder}")
    jcbc=$(echo $jcbc | awk '{print $NF}')
    exit
fi

jbeds=(); lane=0
for file in ${folder}/${lib}*cbc.fastq.gz
do
    lib=${file%_cbc.fastq.gz}
    lib=${lib#*/}
    echo $lib
    ### trim
    jobid=trim-${lib}
    jtrim=2
    jtrim=$(sbatch --export=All -N 1 -J ${jobid} -e ${folder}/${jobid}.err -o ${folder}/${jobid}.out --dependency=afterany:$jcbc -t 05:00:00 --mem=10G --wrap="${p2s}/trim.sh ${folder}/${lib}_cbc.fastq.gz ${folder} ${p2trimgalore} ${p2cutadapt}")
    jtrim=$(echo $jtrim | awk '{print $NF}')

    ### ribo-map
    jobid=ribo-${lib}
    jribo=3
    jribo=$(sbatch --export=All -c 8 -J $jobid -o ${folder}/${jobid}.err -t 10:00:00 --mem=40G --dependency=afterany:$jtrim --wrap="${p2s}/ribo-bwamem.sh $riboref ${folder}/${lib}_cbc_trimmed_homoATCG.fq.gz ${folder}/${lib}_cbc_trimmed_homoATCG $p2bwa $p2samtools y $p2s")
    jribo=$(echo $jribo | awk '{print $NF}')

    ### map to genome
    jobid=gmap-$lib
    jgmap=4
    jgmap=$(sbatch --export=All -c 8 -J $jobid -o ${folder}/${jobid}.err -t 10:00:00 --mem=40G --dependency=afterany:$jribo --wrap="${p2s}/map_star.sh ${p2star} ${p2samtools} ${genome} ${folder}/${lib}_cbc_trimmed_homoATCG.nonRibo.fastq.gz ${folder}/${lib}_cbc_trimmed_homoATCG.nonRibo_E99_")
    jgmap=$(echo $jgmap | awk '{print $NF}')

    jobid=b2bs-$lib
    jb2bs=5
    jb2bs=$(sbatch --export=All -c 1 -N 1 -J $jobid -o ${folder}/${jobid}.err -t 12:00:00 --mem=40G --dependency=afterany:$jgmap --wrap="${p2s}/deal_with_singlemappers.sh ${folder}/${lib}_cbc_trimmed_homoATCG.nonRibo_E99_Aligned.out.bam ${refBED} y ${p2samtools} ${p2bedtools}")
    jb2bs=$(echo $jb2bs | awk '{print $NF}')
    lane=$((lane+1))
    jbeds[$lane]=$jb2bs

    jobid=b2bm-$lib
    jb2bm=6
    jb2bm=$(sbatch --export=All -c 1 -N 1 -J $jobid -o ${folder}/${jobid}.err -t 10:00:00 --mem=40G --dependency=afterany:$jgmap --wrap="${p2s}/deal_with_multimappers.sh ${folder}/${lib}_cbc_trimmed_homoATCG.nonRibo_E99_Aligned.out.bam ${refBED} y ${p2samtools} ${p2bedtools}")
    jb2bm=$(echo $jb2bm | awk '{print $NF}')
    lane=$((lane+1))
    jbeds[$lane]=$jb2bm
done

### count table
jobid=cout-$lib
j=''
for k in ${jbeds[@]}
do
    j=${j}:$k
done
jcout=7
jcout=$(sbatch --export=All -c 8 -t 60:00:00 --mem=160G --dependency=afterany${j} -J cnt${folder} -e cnt${folder}.err -o cnt${folder}.out --mail-type=END --mail-user=${email} --wrap="${p2s}/countTables_2pickle_cellsSpliced.py ${folder} ${folder} vasa $cellidori";)
jcout=$(echo $jcout | awk '{print $NF}')

jobid=pick-$lib
jpick=8
jpick=$(sbatch --export=All -c 1 -t 15:00:00 --mem=160G --dependency=afterany:$jcout -J pick${folder} -e pick${folder}.err -o pick${folder}.out --mail-type=END --mail-user=${email} --wrap="${p2s}/countTables_fromPickle.py ${folder}.pickle.gz ${folder} vasa y") 



