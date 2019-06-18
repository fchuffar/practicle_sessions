import os 

localrules: target

rule target:
    threads: 1
    message: "-- Rule target completed. --"
    input: 
       os.path.expanduser("~/projects/datashare/GSE101086/raw/SRR5815193_1_fastqc.zip"),
       os.path.expanduser("~/projects/datashare/GSE101086/raw/SRR5815194_1_fastqc.zip"),
       os.path.expanduser("~/projects/datashare/GSE101086/raw/SRR5815195_1_fastqc.zip"),
       os.path.expanduser("~/projects/datashare/GSE101086/raw/SRR5815199_1_fastqc.zip"),
       os.path.expanduser("~/projects/datashare/GSE101086/raw/SRR5815200_1_fastqc.zip"),
       os.path.expanduser("~/projects/datashare/GSE101086/raw/SRR5815201_1_fastqc.zip"),

       os.path.expanduser("~/projects/datashare/GSE101086/GSM2699060_notrim_star_Schizosaccharomyces_pombe_ASM294v2_Aligned.sortedByCoord.out.bam"),
       os.path.expanduser("~/projects/datashare/GSE101086/GSM2699061_notrim_star_Schizosaccharomyces_pombe_ASM294v2_Aligned.sortedByCoord.out.bam"),
       os.path.expanduser("~/projects/datashare/GSE101086/GSM2699062_notrim_star_Schizosaccharomyces_pombe_ASM294v2_Aligned.sortedByCoord.out.bam"),
       os.path.expanduser("~/projects/datashare/GSE101086/GSM2699066_notrim_star_Schizosaccharomyces_pombe_ASM294v2_Aligned.sortedByCoord.out.bam"),
       os.path.expanduser("~/projects/datashare/GSE101086/GSM2699067_notrim_star_Schizosaccharomyces_pombe_ASM294v2_Aligned.sortedByCoord.out.bam"),
       os.path.expanduser("~/projects/datashare/GSE101086/GSM2699068_notrim_star_Schizosaccharomyces_pombe_ASM294v2_Aligned.sortedByCoord.out.bam"),
      
       os.path.expanduser("~/projects/datashare/GSE101086/GSM2699060_notrim_star_Schizosaccharomyces_pombe_ASM294v2_genes_strandedreverse_classiccounts.txt"),
       os.path.expanduser("~/projects/datashare/GSE101086/GSM2699061_notrim_star_Schizosaccharomyces_pombe_ASM294v2_genes_strandedreverse_classiccounts.txt"),
       os.path.expanduser("~/projects/datashare/GSE101086/GSM2699062_notrim_star_Schizosaccharomyces_pombe_ASM294v2_genes_strandedreverse_classiccounts.txt"),
       os.path.expanduser("~/projects/datashare/GSE101086/GSM2699066_notrim_star_Schizosaccharomyces_pombe_ASM294v2_genes_strandedreverse_classiccounts.txt"),
       os.path.expanduser("~/projects/datashare/GSE101086/GSM2699067_notrim_star_Schizosaccharomyces_pombe_ASM294v2_genes_strandedreverse_classiccounts.txt"),
       os.path.expanduser("~/projects/datashare/GSE101086/GSM2699068_notrim_star_Schizosaccharomyces_pombe_ASM294v2_genes_strandedreverse_classiccounts.txt"),

    shell:"""
multiqc --force -o  ~/projects/datashare/GSE101086/raw -n multiqc_notrim \
  ~/projects/datashare/GSE101086/raw/*_*_fastqc.zip \
  ~/projects/datashare/GSE101086/*_notrim_star_Schizosaccharomyces_pombe_ASM294v2_Log.final.out \
  
echo workflow \"align_pombe\" completed at `date` 
          """
rule fastqc:
    input:  fastqgz="{prefix}.fastq.gz"
    output: zip="{prefix}_fastqc.zip",
            html="{prefix}_fastqc.html"
    threads: 1
    shell:"fastqc {input.fastqgz}"
    
rule index_genome:
    input:
      genome_fasta= os.path.expanduser("~/projects/datashare/genomes/{species}/UCSC/{index}/Sequence/WholeGenomeFasta/genome.fa"), 
      gtf= os.path.expanduser("~/projects/datashare/genomes/{species}/UCSC/{index}/Annotation/Genes/genes.gtf"),
    output: directory(os.path.expanduser("~/projects/datashare/genomes/{species}/UCSC/{index}/Sequence/StarIndex"))
    threads: 8
    shell:    """
mkdir -p {output}
STAR \
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
      star_index= os.path.expanduser("~/projects/datashare/genomes/{species}/Ensembl/{index}/Sequence/StarIndex"),
      gtf= os.path.expanduser("~/projects/datashare/genomes/{species}/Ensembl/{index}/Annotation/Genes/genes.gtf"),
    output:  "{prefix}/{sample}_{trim}_star_{species}_{index}_Aligned.sortedByCoord.out.bam"
    threads: 8
    shell:"""
cd {wildcards.prefix}
STAR \
  --runThreadN `echo "$(({threads} * 2))"` \
  --genomeDir  {input.star_index} \
  --sjdbGTFfile {input.gtf} \
  --readFilesCommand gunzip -c \
  --readFilesIn `cat {input.fastqc_info}` \
  --outReadsUnmapped Fastx \
  --outFileNamePrefix {wildcards.prefix}/{wildcards.sample}_{wildcards.trim}_star_{wildcards.species}_{wildcards.index}_ \
  --outSAMtype BAM SortedByCoordinate
samtools index {output}
    """
              
rule count_classic_stranded:
    input:
      bam_file="{prefix}/{sample}_{trim}_star_{species}_{index}_Aligned.sortedByCoord.out.bam",
      gtf_file=  os.path.expanduser("~/projects/datashare/genomes/{species}/Ensembl/{index}/Annotation/Genes/{gtf_prefix}.gtf")
    output: "{prefix}/{sample}_{trim}_star_{species}_{index}_{gtf_prefix}_stranded{stranded}_classiccounts.txt"
    priority: 50
    threads: 1
    shell:"""
htseq-count -t exon -f bam -r pos --stranded={wildcards.stranded} -m intersection-strict --nonunique none \
  {input.bam_file} \
  {input.gtf_file} \
  > {output}
    """



rule bigwig_coverage:
    input:
      bam_file="{prefix}/{sample}_{trim}_star_{species}_{index}_Aligned.sortedByCoord.out.bam",
    output: "{prefix}/{sample}_{trim}_star_{species}_{index}_Aligned.sortedByCoord.out.bw"
    threads: 4
    shell:"""
/summer/epistorage/miniconda3/bin/bamCoverage \
  -b {input.bam_file} \
  --numberOfProcessors `echo "$(({threads} * 2))"` \
  --binSize 10 \
  --minMappingQuality 30 \
  --normalizeUsingRPKM \
  -o {output}

    """