def get_files(src_dir, src_suffix, dest_dir, dest_suffix):
  files = [f for f in os.listdir(src_dir) if re.match("^.*"+src_suffix+"$", f)]
  files = [x.replace(src_suffix, dest_suffix) for x in files ]
  return [os.path.join(dest_dir, f) for f in files]


localrules: target

rule target:
    threads: 1
    message: "-- Rule target completed. --"
    input: 
      "/bettik/fchuffar/datashare/GSE101086/GSM2699060_notrim_star_Schizosaccharomyces_pombe_ASM294v2_genes_strandedreverse_classiccounts.txt",
      "/bettik/fchuffar/datashare/GSE101086/GSM2699061_notrim_star_Schizosaccharomyces_pombe_ASM294v2_genes_strandedreverse_classiccounts.txt",
      "/bettik/fchuffar/datashare/GSE101086/GSM2699062_notrim_star_Schizosaccharomyces_pombe_ASM294v2_genes_strandedreverse_classiccounts.txt",
      "/bettik/fchuffar/datashare/GSE101086/GSM2699066_notrim_star_Schizosaccharomyces_pombe_ASM294v2_genes_strandedreverse_classiccounts.txt",
      "/bettik/fchuffar/datashare/GSE101086/GSM2699067_notrim_star_Schizosaccharomyces_pombe_ASM294v2_genes_strandedreverse_classiccounts.txt",
      "/bettik/fchuffar/datashare/GSE101086/GSM2699068_notrim_star_Schizosaccharomyces_pombe_ASM294v2_genes_strandedreverse_classiccounts.txt",

    shell:"""
multiqc --force -o ~/projects/practicle_sessions/data/GSE101086/ -n multiqc_notrim \
  /bettik/fchuffar/datashare/GSE101086/*_notrim_star_Schizosaccharomyces_pombe_ASM294v2_Log.final.out \
  /bettik/fchuffar/datashare/GSE101086/raw/*_*_fastqc.zip \

echo workflow \"align_heatshock\" completed at `date` 
          """
rule fastqc:
    input:  fastqgz="{prefix}.fastq.gz"
    output: zip="{prefix}_fastqc.zip",
            html="{prefix}_fastqc.html"
    threads: 1
    shell:"/summer/epistorage/miniconda3/bin/fastqc {input.fastqgz}"

rule align_trimed:
    input:
      # fqgz_file="{prefix}/{sample}_{trim}.fastq.gz",
      # fastqc_file="{prefix}/{sample}_{trim}_fastqc.zip",
      fastqc_info="{prefix}/{sample}_{trim}_fqgz.info",
      star_index="/home/fchuffar/projects/datashare/genomes/{species}/Ensembl/{index}/Sequence/StarIndex",
      gtf="/home/fchuffar/projects/datashare/genomes/{species}/Ensembl/{index}/Annotation/Genes/genes.gtf",
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
  --outSAMtype BAM SortedByCoordinate
/summer/epistorage/miniconda3/bin/samtools index {output}
    """
              
rule count_classic_stranded:
    input:
      bam_file="{prefix}/{sample}_{trim}_star_{species}_{index}_Aligned.sortedByCoord.out.bam",
      gtf_file= "/home/fchuffar/projects/datashare/genomes/{species}/Ensembl/{index}/Annotation/Genes/{gtf_prefix}.gtf"
    output: "{prefix}/{sample}_{trim}_star_{species}_{index}_{gtf_prefix}_stranded{stranded}_classiccounts.txt"
    priority: 50
    threads: 1
    shell:"""
/summer/epistorage/miniconda3/bin/htseq-count -t exon -f bam -r pos --stranded={wildcards.stranded} -m intersection-strict --nonunique none \
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