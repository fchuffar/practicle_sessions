"""
## Workflow for testing oarsub and snakemake
## This file must be put in your home.
## This workflow must be launched:
##   - as follow,
##   - on luke,
##   - from your home.


## Step 1. jobs run on luke
# cleaning files...
rm test_oarsub_snakemake_file_*
# run snakemake without oarsub...
snakemake -s demo_snakemake_oar.py -p
# get file content... (where jobs run?)
cat test_oarsub_snakemake_file_*


## Step 2. jobs run on computing elements (luke20, luke30...) 2 by 2 (option --cores 2, the maximum of ressources to use.)
# cleaning files...
rm test_oarsub_snakemake_file_*
# run snakemake with oarsub...
snakemake -s demo_snakemake_oar.py --cores 2 --cluster "oarsub --project epimed -l /core={threads},walltime=01:00:00 " -p
# get file content... (where jobs run?)
cat test_oarsub_snakemake_file_*


## Step 3. jobs run on computing elements (luke20, luke30...) alltogether (option --cores 100, limit is bigger than needs...)
# cleaning files...
rm test_oarsub_snakemake_file_*
# run snakemake with oarsub...
snakemake -s demo_snakemake_oar.py --cores 100 --cluster "oarsub --project epimed -l /core={threads},walltime=01:00:00 " -p
# get file content... (where jobs run?)
cat test_oarsub_snakemake_file_*
"""

localrules: target

rule target:
    threads: 1
    message: "-- Rule target completed. --"
    input:
        "test_oarsub_snakemake_file_01.txt",
        "test_oarsub_snakemake_file_02.txt",
        "test_oarsub_snakemake_file_03.txt",
        "test_oarsub_snakemake_file_04.txt",
        "test_oarsub_snakemake_file_05.txt",
        "test_oarsub_snakemake_file_06.txt",
        "test_oarsub_snakemake_file_07.txt",
        "test_oarsub_snakemake_file_08.txt",
        "test_oarsub_snakemake_file_09.txt",
        "test_oarsub_snakemake_file_10.txt"
    shell:"""
echo "workflow \"fch_demo_oar\" completed at `date`" 
          """

rule create_file:
    output:"test_oarsub_snakemake_file_{id}.txt"
    threads: 1
    shell:"echo $HOSTNAME > {output}"