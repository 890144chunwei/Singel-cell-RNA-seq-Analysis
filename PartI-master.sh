#!/usr/bin/bash
#SBATCH --mem 32768
#SBATCH -p short
#SBATCH --output=/storage/goodell/projects/chunweic/slurm_out/220901_demux_%j.out
#SBATCH -e /storage/goodell/projects/chunweic/slurm_out/220901_demux_%j.err # Standard output and error log
#SBATCH --mail-type=ALL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=chunweic@bcm.edu # Email to which notifications will be sent

pwd; hostname; date

PROJECTDIR="/storage/goodell/projects/chunweic/220910_scRNAseq"
cellranger mkfastq --run=/storage/goodell/bcl/220405_chunwei_scRNAseq/Files/ --csv=$PROJECTDIR/SampleSheet_scRNAseq_220831.csv --output-dir=$PROJECTDIR/Fastq2
fastqc $PROJECTDIR/Fastq2/AAAMHM7HV/*.gz
gunzip $PROJECTDIR/Fastq2/AAAMHM7HV/*gz
multiqc $PROJECTDIR/Fastq2/AAAMHM7HV/*zip
cellranger count --id=WT_LSK_1 --transcriptome=/storage/goodell/home/chunweic/cellranger/refdata-gex-mm10-2020-A --fastqs=$PROJECTDIR/Fastq2/AAAMHM7HV/ --sample=WT_LSK_1 --localcores=8 --localmem=64 --r1-length=28
cellranger count --id=WT_LSK_2 --transcriptome=/storage/goodell/home/chunweic/cellranger/refdata-gex-mm10-2020-A --fastqs=$PROJECTDIR/Fastq2/AAAMHM7HV/ --sample=WT_LSK_2 --localcores=8 --localmem=64 --r1-length=28
cellranger count --id=Mut_LSK_1 --transcriptome=/storage/goodell/home/chunweic/cellranger/refdata-gex-mm10-2020-A --fastqs=$PROJECTDIR/Fastq2/AAAMHM7HV/ --sample=Mut_LSK_1 --localcores=8 --localmem=64 --r1-length=28
cellranger count --id=Mut_LSK_2 --transcriptome=/storage/goodell/home/chunweic/cellranger/refdata-gex-mm10-2020-A --fastqs=$PROJECTDIR/Fastq2/AAAMHM7HV/ --sample=Mut_LSK_2 --localcores=8 --localmem=64 --r1-length=28

cd $PROJECTDIR/Align/
cellranger aggr --id=Combined_LSK --csv=$PROJECTDIR/Cellranger_aggr_220831.csv
