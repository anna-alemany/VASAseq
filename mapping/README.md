# VASAseq Mapping Pipeline

This same pipeline was used to map 10X and SmartSeq data, changing the proper parameters (such as UMI length, cell-barcode lenght, or protocol strandness). 

### Contents 
- `bc_celseq2.tsv`: list of barcodes used in vasa-plate. Feel free to change this list if you work with a different list of barcodes (which will be the case when working with vasa-drop)
- `concatenator.py`: python script that takes as an input two fastq files (R1 and R2), extracts the UFI and the barcode information (provided in `bc_celseq2.tsv`) and generates filtered an annotated fastq files for each barcode. Barcode information is present in the file name. In addition, barcode and UFI information (sequence and phred-score) are appended under each read name after the tags 'SS' (sequenced cell barcode),'CB' (assigned cell barcode sequence),'QT' (phred-score for cell barcode,'RX' (UFI sequence),'RQ' (UFI pref-score), and 'SM' (cellID assignment provided in `bc_cellseq2.tsv`). This script has several customizable options and it will work for different barcode and UFI lengths, and if they are present in R1 or R2, etc. For now, default input parameters are set to vasa-plate requirements.
- `countTables_2pickle_cellsSpliced.py`: python script that processes all the bed files produced by `deal_with_multimappers.sh` and `deal_with_singlemappers.sh` in order to produce a dictionary. That dictionary contains the number of observed UFIs for each cell and each annotation.
- `countTables_fromPickle.py`: python script that reads the dictionary produced by `countTables_2pickle_cellsSpliced.py` and produces several count tables (spliced, unspliced, total, tRNA transcripts).
- `deal_with_multimappers.sh`: bash script that assigns annotation to multiple-mappers and produces an output bed file.
- `deal_with_singlemappers.sh`: bash script that assigns annotations to single-mappers and produces an output bed file.
- `extractBC.sh`: bash wrapper that calls the `concatenator.py` script with different input options.
- `map_star.sh`: bash wrapper to map a fastq file to a reference genome using STAR.
- `resque_unmapped_fromBAM.sh`: bash script, not included in the main pipeline, used to recover unmapped reads form an input BAM file and generate a new fastq file.
- `ribo-bwamem.sh`: wrapper script that performs _in silico_ ribosomal depletion, by first maping an input fastq file to a ribosomal reference and then filtering reads that do not map.
- `riboread-selection.py`: python script that selects reads from an input BAMfile that are never mapped (used for the _in silico_ ribosomal depletion). Multi-mappers are also filtered out. This script runs for both stranded and not-stranded data.
- `submit_vasaplate_map.sh`
- `trim.sh`: bash script that takes an input fastq file and produces an output fastq file in which illumina adaptors and polyG, polyA, polyT and polyC stretches are trimmed out from the 3' end of the read.

### Quick run of the pipeline for VASA-plate
If you are dealing with **VASA-plate**, save all the scripts in the same folder. Then, open the `submit_vasaplate_map.sh` script and provide all paths to important software used by the pipeline and your email. Then run: 
```{bash}
path_2_script=set_path_to_mapping_scripts
${path_2_script}/submit_vasaplate_map.sh ${prefix_input_fastqfiles} ${reference_genome} ${read_length} ${output_directory} ${prefix_output_count_tables} y r
```
Once this is done (you should get an email), run: 
```{bash}
${path_2_script}/submit_vasaplate_map.sh ${prefix_input_fastqfiles} ${reference_genome} ${read_length} ${output_directory} ${prefix_output_count_tables} n r
```
(notice, that only one argument changes from `y` to `n`)

### Steps and scripts involved in the pipeline
1. Fastq file pre-processing

The pipeline creates one fastq for each barcode, in which the barcode and UFI sequence and phred-scores of each read are stored in the read name. This is done by the script `concatenator.py`, as follows:
```{bash}
python3 concatenator.py --fqf ${prefix_input_fastq_files} --cbcfile bc_celseq2.tsv --cbchd 0 --lenumi 6 --umifirst --demux --outdir ${name_output_directory}

# --fqf: mandatory argument that provides the input fastq files. All fastq files that start their name as ${prefix_input_fastq_files} will be used to generate the output fastq files. 
# --cbcfile: mandatory argument that provides the input file with the barcodes. For VASA-plate, this file is bc_celseq2.tsv. For VASA-drop, you need to produce the barcode file for each library. 
# --cbchd: optional argument to set number of accepted mismatches between sequenced barcode and reference barcode
# --lenumi: UFI length
# --demux: when given, this optional parameter makes sure that each barcode gets its own fastq file. When not provided, all reads are saved in the same fastq file. 
# --umifirst: logical argument that tells the script whether the UFI goes before the cell barcode in the read
# --outdir: argument that sets the outptu directory for the output fastq files
```
Alternatively, the user can call the wrapper `extracth.sh`:
```{bash}
extractBC.sh ${prefix_input_fastq_file} vasaplate ${path_2_script} ${name_output_directory}
```
where `path_2_script` in the path to the folder with the script (scripts can be saved in a different folder from the fastq files).

2. Fastq file trimming

Due to technical reasons, the VASA protocol can generate polyA tails at the end of the read. Therefore, we trim homopolymers away. We use `trim_galore` to trim Illumina adatpros and `cutadapt` to trim homopolymers. These 2 steps are included in the bash wrapper script `trim.sh`, which needs to be called for each of the fastq files generated in the previous step:
```{bash}
trim.sh ${input_cbcfastq_file} ${output_directory} ${path_2_trimgalore} ${path_2_cutadapt}

# input_cbcfastq_file: this is the input fastq file, that was generated in the previous step
```

3. _In silico_ ribosomal depletion

Like any other total RNA-seq protocol, VASA includes an experimental ribosomal-depletion step. However, this is not always perfect and a lot of ribosomal reads survive and get sequenced. It is good to remove them before mapping the library to the genome, since they generate mapping artifacts that can be confused by false small non-coding RNA.

To deplete ribosomal reads _in silico_ we map our trimmed fastq files to a fasta file containing the reference sequences for the ribosomal RNA (_e.g._ Rn45s, Rn6s, 12s, 16s, 47s for mouse, 12S, 16S, 45SN1-5, 45S9, 45S1-17 for human). To do so, we use `bwa`. Because ribosomal reads span a lot of read lengths, we found that both `bwa aln` and `bwa mem` need to be used for a successful ribosomal depletion.  This step results in two BAM files (one for aln, another for mem), which are merged into one that subsequently is sorted by read name. Finally, only reads that remain unmapped with the two approaches are kept for downstream processing (script `riboread-selection.py`). Therefore, both single and multi-mappers are filtered out. These steps are all performed in the bash wrapper script `ribo-bwamem.sh`, which can be called as:
```{bash}
ribo-bwamem.sh $riboref ${input_cbc_trimmed_fastqfile} ${prefix_output_file} ${path_2_bwa} ${path_2_samtools} ${standed_protocol} ${path_2_scripts}
```
The parameter `${standed_protocol}` should be set to "y" when dealing with VASA data. It indicates whether the protocol is stranded.

4. Mapping to the genome

We map the ribosomal-depleted-trimmed-annotated fastq files to a reference genome using STAR:
```{bash}
STAR --runThreadN 8 --genomeDir ${genome} --readFilesIn ${input_fastq_file} --readFilesCommand zcat --outFilterMultimapNmax 20 --outSAMunmapped Within --outSAMtype BAM Unsorted --outSAMattributes All  --outFileNamePrefix ${outprefix}
```
For this, we have a bash wrapper `map_star.sh`.

5. Count tables generation

The generation of count tables highly relies on a complete `gtf` file that contains annotated protein-coding genes, but also the non-coding transcripts. such a `gtf` file does not exist and needs to be prepared by the user.
Here, we integrated `gtf` files containing protein-coding and some long-non-coding transcripts from ENSEMBL with other specialized gtf file specific for small non-coding biotypes such as miRNA, snoRNA, etc, and also tRNA. We reduced as much as possible multi-annotations detected in identical genomic regions.
Next, we converted such `gtf` file into a `bed`  file. For protein-coding and long-non-coding, we had different entries for exonic and intronic fragments, and included biotype information, total gene length, and initial and final position of the whole gene. This reads as:
```{bash}
1       3073253 3074322 +       ENSMUSG00000102693_4933401J01Rik_TEC_exon       1069    3073253 3074322
1       3102016 3102125 +       ENSMUSG00000064842_Gm26206_snRNA_exon   109     3102016 3102125
1       3205901 3207317 -       ENSMUSG00000051951_Xkr4_ProteinCoding_exon      465597  3205901 3671498
1       3207318 3213438 -       ENSMUSG00000051951_Xkr4_ProteinCoding_intron    465597  3205901 3671498
1       3213439 3216968 -       ENSMUSG00000051951_Xkr4_ProteinCoding_exon      465597  3205901 3671498
1       3216969 3421701 -       ENSMUSG00000051951_Xkr4_ProteinCoding_intron    465597  3205901 3671498
... (etc) ...
```

In this step of the pipeline, we go throught the reads of the BAM file, and we check if they fall inside one or more annotated genomic region. In most of the cases, each read is assigned to a unique annotation. However, sometimes a single-mapper read falls into a genomic region where we have multiple annotations (_e.g._  a miRNA located within an intron of a protein-coding gene). In addition, multiple-mappers also get several annotations, from the different regions they map to. In several occasions, we found that multiple annotations for the same read reveal transcriptonal regulatroy mechanisms (_e.g._ we find that tRNA reads are multimappers, since the same tRNA can be transcribed from multiple regions in the genome). Taken together, we established the following hierarchy to assign multi-annotated reads:

* Only multi-annotations with best mapping scores and less mismatches are considered.
* if a read fully falls within a small non-coding RNA region, other annotations are ignored.
* If a read fully falls within more than one small non-coding RNA region, the shortest is given priority.
* Purely exonic annotations are given priority over purely intronic or intron-exon junctions.
* Identical annotations are summarized to a singe one, whose name is the union of all the small non-coding annotations (sorted alphabetically).
* At the end, UFI information is combined in order to properly annotate spliced and unspliced transcripts.

The resulting table might contain some entries that look like `geneA_geneB`. If required, when we find either geneA or geneB expressed alone, the entries attribued to the gene combination are assigned to the uniquely expressed gene. If both geneA and geneB are also found to be expressed alone, then the entry is left as combination of two genes.

All these steps are performed by the script wrappers `deal_with_singlemappers.sh`, `deal_with_multimappers.sh` and `countTables_2pickle_cellsSpliced.py` and `countTables_fromPickle.py`. The first two (`deal_with_*mappers.sh`) run for each cell separatedly:
```{bash}
deal_with_singlemappers.sh ${input_bam_file} ${refBED} ${stranded}
deal_with_multimappers.sh ${input_bam_file} ${refBED} ${stranded}
```
where ${stranded} is set to "y" in vasa data.

The second two (`countTables*py`) reads all the bed files obtained for all the cells and produces the final count tables. First:
```{bash}
countTables_2pickle_cellsSpliced.py ${input_folder} ${prefix_output} ${protocol}$ ${cellIDorigin}"
```
here, for vasa-plate and vasa-drop, ${protocol} should be set to "vasa". ${cellIDorigin} has the options "r" and "f", standing for "read" and "file", respectively. That is, if "r" is given, the cell ID will be extracted from the read name. If "f" is given, the cell ID will be extracted from the file name. This is done in case the user wants to combine several experiments at once in this step already.
The second script can be called as:
```{bash}
countTables_fromPickle.py ${prefix_output}.pickle.gz ${prefix_output} ${protocol}$ ${filter_genes}
```
It takes the `pickle.gz` file generated by the previous script. As before, ${protocol} should be set to "vasa". The last label, ${filter_genes}, refers to the possibility of colapsing multiple gene anntotations to single ones based on expression patterns (explained above).
