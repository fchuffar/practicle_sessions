"""
QC, alignment and count of RNA-seq data 
Celine Mandier, Florent Chuffart (EpiMed/IAB/INSERM)
17 july 2019
"""

localrules: target

rule target:
    threads: 1
    message: "-- Rule target completed. --"
    input: 
      os.path.expanduser("~/projects/datashare/GSE45332/raw/SRR787296_1_fastqc.zip"),
      os.path.expanduser("~/projects/datashare/GSE45332/raw/SRR787297_1_fastqc.zip"),
      os.path.expanduser("~/projects/datashare/GSE45332/raw/SRR787298_1_fastqc.zip"),
      os.path.expanduser("~/projects/datashare/GSE45332/raw/SRR787299_1_fastqc.zip"),
      os.path.expanduser("~/projects/datashare/GSE45332/raw/SRR787300_1_fastqc.zip"),
      os.path.expanduser("~/projects/datashare/GSE45332/raw/SRR787301_1_fastqc.zip"),
      os.path.expanduser("~/projects/datashare/GSE45332/raw/SRR787302_1_fastqc.zip"),
      os.path.expanduser("~/projects/datashare/GSE45332/raw/SRR787303_1_fastqc.zip"),
      os.path.expanduser("~/projects/datashare/GSE45332/raw/SRR787304_1_fastqc.zip"),
      os.path.expanduser("~/projects/datashare/GSE45332/raw/SRR787305_1_fastqc.zip"),
      
      os.path.expanduser("~/projects/datashare/GSE45332/GSM1102624_notrim_star_Homo_sapiens_hg19_geneswchrm_strandedyes_classiccounts.txt"),
      os.path.expanduser("~/projects/datashare/GSE45332/GSM1102625_notrim_star_Homo_sapiens_hg19_geneswchrm_strandedyes_classiccounts.txt"),
      os.path.expanduser("~/projects/datashare/GSE45332/GSM1102626_notrim_star_Homo_sapiens_hg19_geneswchrm_strandedyes_classiccounts.txt"),
      os.path.expanduser("~/projects/datashare/GSE45332/GSM1102627_notrim_star_Homo_sapiens_hg19_geneswchrm_strandedyes_classiccounts.txt"),
      os.path.expanduser("~/projects/datashare/GSE45332/GSM1102628_notrim_star_Homo_sapiens_hg19_geneswchrm_strandedyes_classiccounts.txt"),
      os.path.expanduser("~/projects/datashare/GSE45332/GSM1102629_notrim_star_Homo_sapiens_hg19_geneswchrm_strandedyes_classiccounts.txt"),
      os.path.expanduser("~/projects/datashare/GSE45332/GSM1102630_notrim_star_Homo_sapiens_hg19_geneswchrm_strandedyes_classiccounts.txt"),
      os.path.expanduser("~/projects/datashare/GSE45332/GSM1102631_notrim_star_Homo_sapiens_hg19_geneswchrm_strandedyes_classiccounts.txt"),
      os.path.expanduser("~/projects/datashare/GSE45332/GSM1102632_notrim_star_Homo_sapiens_hg19_geneswchrm_strandedyes_classiccounts.txt"),
      os.path.expanduser("~/projects/datashare/GSE45332/GSM1102633_notrim_star_Homo_sapiens_hg19_geneswchrm_strandedyes_classiccounts.txt"),
      
    shell:"""
multiqc --force -o ~/projects/datashare/GSE45332/raw -n multiqc_notrim \
  ~/projects/datashare/GSE45332/*_notrim_star_Homo_sapiens_hg19_Log.final.out \
  ~/projects/datashare/GSE45332/raw/*_*_fastqc.zip 

echo workflow \"align_GSE45332\" completed at `date` 
          """
rule fastqc:
    input:  fastqgz="{prefix}.fastq.gz"
    output: zip="{prefix}_fastqc.zip",
            html="{prefix}_fastqc.html"
    threads: 1
    shell:"""
    export PATH="/summer/epistorage/miniconda3/bin:$PATH"
    /summer/epistorage/miniconda3/bin/fastqc {input.fastqgz}
    """

rule index_genome:
    input:
      genome_fasta=os.path.expanduser("~/projects/datashare/genomes/{species}/UCSC/{index}/Sequence/WholeGenomeFasta/genome.fa"), 
      gtf=os.path.expanduser("~/projects/datashare/genomes/{species}/UCSC/{index}/Annotation/Genes/geneswchrm.gtf"),
    output: directory(os.path.expanduser("~/projects/datashare/genomes/{species}/UCSC/{index}/Sequence/StarIndex"))
    #priority: 0
    threads: 8
    shell:"""
mkdir -p {output}
/summer/epistorage/miniconda3/bin/STAR \
  --runThreadN `echo "$(({threads} * 2))"` \
  --runMode genomeGenerate \
  --genomeDir {output} \
  --genomeFastaFiles  {input.genome_fasta} \
  --sjdbGTFfile {input.gtf} \
  --sjdbOverhang 100
    """

rule align_trimed:
    input:
      fastqc_info="{prefix}/{sample}_{trim}_fqgz.info",
      star_index=os.path.expanduser("~/projects/datashare/genomes/{species}/UCSC/{index}/Sequence/StarIndex"),
      gtf=os.path.expanduser("~/projects/datashare/genomes/{species}/UCSC/{index}/Annotation/Genes/geneswchrm.gtf"),
    output:  "{prefix}/{sample}_{trim}_star_{species}_{index}_Aligned.sortedByCoord.out.bam"
    threads: 8
    shell:"""
cd {wildcards.prefix}
/summer/epistorage/miniconda3/bin/STAR \
  --runThreadN `echo "$(({threads} * 2))"` \
  --genomeDir  {input.star_index} \
  --sjdbGTFfile {input.gtf} \
  --readFilesCommand gunzip -c \
  --readFilesIn `cat {input.fastqc_info}` \
  --outFileNamePrefix {wildcards.prefix}/{wildcards.sample}_{wildcards.trim}_star_{wildcards.species}_{wildcards.index}_ \
  --outReadsUnmapped Fastx \
  --outSAMtype BAM SortedByCoordinate
/summer/epistorage/miniconda3/bin/samtools index {output}
    """
              
rule count_classic:
    input:
      bam_file="{prefix}/{sample}_{trim}_star_{species}_{index}_Aligned.sortedByCoord.out.bam",
      gtf_file= os.path.expanduser("~/projects/datashare/genomes/{species}/UCSC/{index}/Annotation/Genes/{gtf_prefix}.gtf")
    output: "{prefix}/{sample}_{trim}_star_{species}_{index}_{gtf_prefix}_stranded{stranded}_classiccounts.txt"
    priority: 50
    threads: 1
    shell:"""
/summer/epistorage/miniconda3/bin/htseq-count -t exon -f bam -r pos --stranded={wildcards.stranded} -m intersection-strict --nonunique none \
  {input.bam_file} \
  {input.gtf_file} \
  > {output}
    """

