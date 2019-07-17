## Preprocessing of RNA-seq data 
## Celine Mandier, Florent Chuffart (EpiMed/IAB/INSERM)
## 16 july 2019

## Data description
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE45332

## Design
# SRR787296 GSM1102624 	HCT116 total RNA 5h DSN
# SRR787297 GSM1102625 	HCT116 total RNA 6h DSN
# SRR787298 GSM1102626 	HCT116 total RNA 7h DSN
# SRR787299 GSM1102627 	HCT116 total RNA 8h DSN
# SRR787300 GSM1102628 	HCT116 total RNA 9h DSN
# SRR787301 GSM1102629 	HCT116 total RNA 10h DSN
# SRR787302 GSM1102630 	HCT116 total RNA 11h DSN
# SRR787303 GSM1102631 	HCT116 polyA_mRNA
# SRR787304 GSM1102632 	DKO total RNA
# SRR787305 GSM1102633 	DKO polyA mRNA

## Download fastq files from NCBI Geo
mkdir -p ~/projects/datashare/GSE45332/raw
cd ~/projects/datashare/GSE45332/raw
parallel-fastq-dump --threads 8 --tmpdir ./ --gzip --split-files --outdir ./ --sra-id SRR787296
parallel-fastq-dump --threads 8 --tmpdir ./ --gzip --split-files --outdir ./ --sra-id SRR787297
parallel-fastq-dump --threads 8 --tmpdir ./ --gzip --split-files --outdir ./ --sra-id SRR787298
parallel-fastq-dump --threads 8 --tmpdir ./ --gzip --split-files --outdir ./ --sra-id SRR787299
parallel-fastq-dump --threads 8 --tmpdir ./ --gzip --split-files --outdir ./ --sra-id SRR787300
parallel-fastq-dump --threads 8 --tmpdir ./ --gzip --split-files --outdir ./ --sra-id SRR787301
parallel-fastq-dump --threads 8 --tmpdir ./ --gzip --split-files --outdir ./ --sra-id SRR787302
parallel-fastq-dump --threads 8 --tmpdir ./ --gzip --split-files --outdir ./ --sra-id SRR787303
parallel-fastq-dump --threads 8 --tmpdir ./ --gzip --split-files --outdir ./ --sra-id SRR787304
parallel-fastq-dump --threads 8 --tmpdir ./ --gzip --split-files --outdir ./ --sra-id SRR787305

## SR or PE?
ls -lha ~/projects/datashare/GSE45332/raw
sequencing_read_type=SR


## Pipeline "by hand""
### QC
# Calculation of quality of the reads  
cd ~/projects/datashare/GSE45332/raw
fastqc SRR787296_1.fastq.gz

### Index genome
# Create a new StarIndex directory that will group the different indexes created  
# Launch of STAR software 
# Set the number of cores according to the server used
# Command to generate the genome index  
# Specify the direction in which the directory is located which will contain the saved indexes 
# Recovery of the reference genome 
# Recovery of annotations 
# Specifies the length of the annotated genomic sequence
mkdir -p ~/projects/datashare/genomes/Homo_sapiens/UCSC/hg19/Sequence/StarIndex
STAR \
  --runThreadN 8 \
  --runMode genomeGenerate \
  --genomeDir ~/projects/datashare/genomes/Homo_sapiens/UCSC/hg19/Sequence/StarIndex \
  --genomeFastaFiles  ~/projects/datashare/genomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa \
  --sjdbGTFfile ~/projects/datashare/genomes/Homo_sapiens/UCSC/hg19/Annotation/Genes/geneswchrm.gtf \
  --sjdbOverhang 100

### Align reads
# To give the direction of STAR software 
# Launch of STAR software 
# Set the number of cores according to the server used 
# Recovery of annotations 
# Decompression of file containing  the sequence to map  
# Reading of file with sequences to map 
# To define the name of new file including the alignment 
# To indicate the exit of no map reads of the alignment file 
# Files BAM sorted thanks to their coordinates 
cd ~/projects/datashare/GSE45332/
STAR \
  --runThreadN 8 \
  --genomeDir  ~/projects/datashare/genomes/Homo_sapiens/UCSC/hg19/Sequence/StarIndex \
  --sjdbGTFfile ~/projects/datashare/genomes/Homo_sapiens/UCSC/hg19/Annotation/Genes/geneswchrm.gtf \
  --readFilesCommand gunzip -c \
  --readFilesIn raw/SRR787296_1.fastq.gz \
  --outFileNamePrefix GSM1102624_notrim_star_Homo_sapiens_hg19_ \
  --outReadsUnmapped Fastx \
  --outSAMtype BAM SortedByCoordinate
samtools index GSM1102624_notrim_star_Homo_sapiens_hg19_Aligned.sortedByCoord.out.bam

### Count
# Count of the reads wich have positive strand for each position range of the gene  
# Recovery of file with reads alignment 
# Recovery of file with annotations 
# Retranscription of the number of reads for each position range 
htseq-count -t exon -f bam -r pos --stranded=yes -m intersection-strict --nonunique none ~/projects/datashare/GSE45332/GSM1102624_notrim_star_Homo_sapiens_hg19_Aligned.sortedByCoord.out.bam ~/projects/datashare/genomes/Homo_sapiens/UCSC/hg19/Annotation/Genes/geneswchrm.gtf > ~/projects/datashare/GSE45332/GSM1102624_notrim_star_Homo_sapiens_hg19_geneswchrm_strandedyes_classiccounts.txt




## Pipeline using snakemake
### Linking sample and raw files
cd ~/projects/datashare/GSE45332/
echo raw/SRR787296_1.fastq.gz > GSM1102624_notrim_fqgz.info
echo raw/SRR787297_1.fastq.gz > GSM1102625_notrim_fqgz.info
echo raw/SRR787298_1.fastq.gz > GSM1102626_notrim_fqgz.info
echo raw/SRR787299_1.fastq.gz > GSM1102627_notrim_fqgz.info
echo raw/SRR787300_1.fastq.gz > GSM1102628_notrim_fqgz.info
echo raw/SRR787301_1.fastq.gz > GSM1102629_notrim_fqgz.info
echo raw/SRR787302_1.fastq.gz > GSM1102630_notrim_fqgz.info
echo raw/SRR787303_1.fastq.gz > GSM1102631_notrim_fqgz.info
echo raw/SRR787304_1.fastq.gz > GSM1102632_notrim_fqgz.info
echo raw/SRR787305_1.fastq.gz > GSM1102633_notrim_fqgz.info

### snakemake on local machine (it needs at least 32Go ok memory)
cd ~/projects/practicle_sessions/polya_signature/
snakemake -s wf.py --dag | dot -Tpng > dag.png
snakemake -s wf.py --cores 2 -pn

### Put wf on luke and launch it (on a luke node or on dahu cluster)
rsync -auvP ~/projects/practicle_sessions/polya_signature/ luke:~/projects/practicle_sessions/polya_signature/
snakemake -s ~/projects/practicle_sessions/polya_signature/wf.py --cores 8 -pn
snakemake -s ~/projects/practicle_sessions/polya_signature/wf.py --cores 49 --cluster "oarsub --project epimed -l nodes=1/core={threads},walltime=6:00:00 " -pn

### Get results
mkdir -p ~/projects/datashare/GSE45332/raw/
rsync -auvP luke:~/projects/datashare/GSE45332/raw/*.html ~/projects/datashare/GSE45332/raw/.
rsync -auvP luke:~/projects/datashare/GSE45332/raw/multiqc_notrim* ~/projects/datashare/GSE45332/raw/
rsync -auvP luke:~/projects/datashare/GSE45332/*.txt ~/projects/datashare/GSE45332/.

