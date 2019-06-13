## data description
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE101086
#
# SRR5815193 SRX2993641 GSM2699060 WT_RNAseq_rep1
# SRR5815194 SRX2993642 GSM2699061 WT_RNAseq_rep2
# SRR5815195 SRX2993643 GSM2699062 WT_RNAseq_rep3
# SRR5815199 SRX2993647 GSM2699066 pRnhA_RNAseq_rep1
# SRR5815200 SRX2993648 GSM2699067 pRnhA_RNAseq_rep2
# SRR5815201 SRX2993649 GSM2699068 pRnhA_RNAseq_rep3

## download fastq files
mkdir -p ~/projects/datashare/GSE101086/raw
cd ~/projects/datashare/GSE101086/raw
parallel-fastq-dump --threads 4 --tmpdir ./ --gzip --split-files --outdir ./ --sra-id SRR5815193
parallel-fastq-dump --threads 4 --tmpdir ./ --gzip --split-files --outdir ./ --sra-id SRR5815194
parallel-fastq-dump --threads 4 --tmpdir ./ --gzip --split-files --outdir ./ --sra-id SRR5815195
parallel-fastq-dump --threads 4 --tmpdir ./ --gzip --split-files --outdir ./ --sra-id SRR5815199
parallel-fastq-dump --threads 4 --tmpdir ./ --gzip --split-files --outdir ./ --sra-id SRR5815200
parallel-fastq-dump --threads 4 --tmpdir ./ --gzip --split-files --outdir ./ --sra-id SRR5815201
# SR or PE?
ls -lha ~/projects/datashare/GSE101086/raw

# QC
fastqc SRR5815193_1.fastq.gz

# index genome
mkdir -p ~/projects/datashare/genomes/Schizosaccharomyces_pombe/Ensembl/ASM294v2/Sequence/StarIndex
STAR \
  --runThreadN 4 \
  --runMode genomeGenerate \
  --genomeDir ~/projects/datashare/genomes/Schizosaccharomyces_pombe/Ensembl/ASM294v2/Sequence/StarIndex \
  --genomeFastaFiles  ~/projects/datashare/genomes/Schizosaccharomyces_pombe/Ensembl/ASM294v2/Sequence/WholeGenomeFasta/genome.fa \
  --sjdbGTFfile ~/projects/datashare/genomes/Schizosaccharomyces_pombe/Ensembl/ASM294v2/Annotation/Genes/genes.gtf \
  --sjdbOverhang 100


# align reads
cd ~/projects/datashare/GSE101086/
STAR \
  --runThreadN 4 \
  --genomeDir  ~/projects/datashare/genomes/Schizosaccharomyces_pombe/Ensembl/ASM294v2/Sequence/StarIndex \
  --sjdbGTFfile ~/projects/datashare/genomes/Schizosaccharomyces_pombe/Ensembl/ASM294v2/Annotation/Genes/genes.gtf \
  --readFilesCommand gunzip -c \
  --readFilesIn raw/SRR5815193_1.fastq.gz \
  --outFileNamePrefix GSM2699060_notrim_star_Schizosaccharomyces_pombe_ASM294v2_ \
  --outReadsUnmapped Fastx \
  --outSAMtype BAM SortedByCoordinate
samtools index 	_notrim_star_Schizosaccharomyces_pombe_ASM294v2_Aligned.sortedByCoord.out.bam



## qc align count
cd ~/projects/datashare/GSE101086/
echo raw/SRR5815193_1.fastq.gz > GSM2699060_notrim_fqgz.info
echo raw/SRR5815194_1.fastq.gz > GSM2699061_notrim_fqgz.info
echo raw/SRR5815195_1.fastq.gz > GSM2699062_notrim_fqgz.info
echo raw/SRR5815199_1.fastq.gz > GSM2699066_notrim_fqgz.info
echo raw/SRR5815200_1.fastq.gz > GSM2699067_notrim_fqgz.info
echo raw/SRR5815201_1.fastq.gz > GSM2699068_notrim_fqgz.info

# put wf on luke and luachn 

rsync -auvP ~/projects/practicle_sessions/data/GSE101086/ luke:~/projects/practicle_sessions/data/GSE101086/
snakemake -s ~/projects/practicle_sessions/data/GSE101086/wf.py --cores 8 -pn


## get results
mkdir -p ~/projects/datashare/GSE101086/raw/
rsync -auvP luke:~/projects/datashare/GSE101086/raw/*.html ~/projects/datashare/GSE101086/raw/
rsync -auvP luke:~/projects/datashare/GSE101086/*.txt ~/projects/datashare/GSE101086/
rsync -auvP luke:~/projects/practicle_sessions/data/GSE101086/multiqc_notrim* ~/projects/practicle_sessions/data/GSE101086/.