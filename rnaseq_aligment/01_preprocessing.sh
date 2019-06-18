# # Data description
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE101086

# ## Design
# SRR5815193 SRX2993641 GSM2699060 WT_RNAseq_rep1
# SRR5815194 SRX2993642 GSM2699061 WT_RNAseq_rep2
# SRR5815195 SRX2993643 GSM2699062 WT_RNAseq_rep3
# SRR5815199 SRX2993647 GSM2699066 pRnhA_RNAseq_rep1
# SRR5815200 SRX2993648 GSM2699067 pRnhA_RNAseq_rep2
# SRR5815201 SRX2993649 GSM2699068 pRnhA_RNAseq_rep3

## Download fastq files from NCBI Geo
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


# # Pipeline "by hand""
# ## QC
# Calculation of quality of the reads  
fastqc SRR5815193_1.fastq.gz

# ## Index genome
# Create a new StarIndex directory that will group the different indexes created  
# Launch of STAR software 
# Set the number of cores according to the server used
# Command to generate the genome index  
# Specify the direction in which the directory is located which will contain the saved indexes 
# Recovery of the reference genome 
# Recovery of annotations 
# Specifies the length of the annotated genomic sequence
mkdir -p ~/projects/datashare/genomes/Schizosaccharomyces_pombe/Ensembl/ASM294v2/Sequence/StarIndex
STAR \
  --runThreadN 4 \
  --runMode genomeGenerate \
  --genomeDir ~/projects/datashare/genomes/Schizosaccharomyces_pombe/Ensembl/ASM294v2/Sequence/StarIndex \
  --genomeFastaFiles  ~/projects/datashare/genomes/Schizosaccharomyces_pombe/Ensembl/ASM294v2/Sequence/WholeGenomeFasta/genome.fa \
  --sjdbGTFfile ~/projects/datashare/genomes/Schizosaccharomyces_pombe/Ensembl/ASM294v2/Annotation/Genes/genes.gtf \
  --sjdbOverhang 100

# ## Align reads
# To give the direction of STAR software 
# Launch of STAR software 
# Set the number of cores according to the server used 
# Recovery of annotations 
# Decompression of file containing  the sequence to map  
# Reading of file with sequences to map 
# To define the name of new file including the alignment 
# To indicate the exit of no map reads of the alignment file 
# Files BAM sorted thanks to their coordinates 
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
samtools index GSM2699060_notrim_star_Schizosaccharomyces_pombe_ASM294v2_Aligned.sortedByCoord.out.bam

# ## Count
# Count of the reads wich have positive sense for each position range of the gene  
# Recovery of file with reads alignment 
# Recovery of file with annotations 
# Retranscription of the number of reads for each position range 
htseq-count -t exon -f bam -r pos --stranded=yes -m intersection-strict --nonunique none ~/projects/datashare/GSE101086/GSM2699060_notrim_star_Schizosaccharomyces_pombe_ASM294v2_Aligned.sortedByCoord.out.bam ~/projects/datashare/genomes/Schizosaccharomyces_pombe/Ensembl/ASM294v2/Annotation/Genes/genes.gtf > ~/projects/datashare/GSE101086/GSM2699060_notrim_star_Schizosaccharomyces_pombe_ASM294v2_genes_strandedyes_classiccounts.txt

htseq-count -t exon -f bam -r pos --stranded=no -m intersection-strict --nonunique none ~/projects/datashare/GSE101086/GSM2699060_notrim_star_Schizosaccharomyces_pombe_ASM294v2_Aligned.sortedByCoord.out.bam ~/projects/datashare/genomes/Schizosaccharomyces_pombe/Ensembl/ASM294v2/Annotation/Genes/genes.gtf > ~/projects/datashare/GSE101086/GSM2699060_notrim_star_Schizosaccharomyces_pombe_ASM294v2_genes_strandedno_classiccounts.txt

htseq-count -t exon -f bam -r pos --stranded=reverse -m intersection-strict --nonunique none ~/projects/datashare/GSE101086/GSM2699060_notrim_star_Schizosaccharomyces_pombe_ASM294v2_Aligned.sortedByCoord.out.bam ~/projects/datashare/genomes/Schizosaccharomyces_pombe/Ensembl/ASM294v2/Annotation/Genes/genes.gtf > ~/projects/datashare/GSE101086/GSM2699060_notrim_star_Schizosaccharomyces_pombe_ASM294v2_genes_strandedreverse_classiccounts.txt
# Stranded or not?



# # Pipeline using snakemake
## Linking sample and raw files
cd ~/projects/datashare/GSE101086/
echo raw/SRR5815193_1.fastq.gz > GSM2699060_notrim_fqgz.info
echo raw/SRR5815194_1.fastq.gz > GSM2699061_notrim_fqgz.info
echo raw/SRR5815195_1.fastq.gz > GSM2699062_notrim_fqgz.info
echo raw/SRR5815199_1.fastq.gz > GSM2699066_notrim_fqgz.info
echo raw/SRR5815200_1.fastq.gz > GSM2699067_notrim_fqgz.info
echo raw/SRR5815201_1.fastq.gz > GSM2699068_notrim_fqgz.info

## snakemake on local machine 
cd ~/projects/practicle_sessions/data/GSE101086/
snakemake -s wf.py --dag 
snakemake -s wf.py --cores 2 -p

# Put wf on luke and launch it
rsync -auvP ~/projects/practicle_sessions/data/GSE101086/ luke:~/projects/practicle_sessions/data/GSE101086/
snakemake -s ~/projects/practicle_sessions/data/GSE101086/wf.py --cores 8 -pn
snakemake -s ~/projects/practicle_sessions/data/GSE101086/wf.py --cores 49 --cluster "oarsub --project epimed -l nodes=1/core={threads},walltime=6:00:00 " -pn

## get results
mkdir -p ~/projects/datashare/GSE101086/raw/
rsync -auvP luke:~/projects/datashare/GSE101086/raw/*.html ~/projects/datashare/GSE101086/raw/
rsync -auvP luke:~/projects/datashare/GSE101086/*.txt ~/projects/datashare/GSE101086/
rsync -auvP luke:~/projects/practicle_sessions/data/GSE101086/multiqc_notrim* ~/projects/practicle_sessions/data/GSE101086/.