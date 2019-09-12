## data description
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110637
#
# SRX3700627 SRR6727625 GSM3004546 TALL_JS_2 (polyA+ RNA)
# SRX3700628 SRR6727626 GSM3004547 TALL_JS_3 (polyA+ RNA)
# SRX3700629 SRR6727627 GSM3004548 TALL_JS_4 (polyA+ RNA)
# SRX3700630 SRR6727628 GSM3004549 TALL_JS_5 (polyA+ RNA)
# SRX3700635 SRR6727633 GSM3004554 TALL_JS_10 (polyA+ RNA)
# SRX3700636 SRR6727634 GSM3004555 TALL_JS_11 (polyA+ RNA)
# SRX3700638 SRR6727636 GSM3004557 TALL_JS_13 (polyA+ RNA)
# SRX3700639 SRR6727637 GSM3004558 TALL_JS_14 (polyA+ RNA)
# SRX3700643 SRR6727641 GSM3004562 TALL_JS_18 (polyA+ RNA)
# SRX3700648 SRR6727646 GSM3004567 TALL_JS_23 (polyA+ RNA)
# SRX3700651 SRR6727649 GSM3004570 TALL_JS_26 (polyA+ RNA)
# SRX3700656 SRR6727654 GSM3004575 TALL_JS_31 (polyA+ RNA)
# SRX3700659 SRR6727657 GSM3004578 TALL_JS_35 (polyA+ RNA)
# SRX3700660 SRR6727658 GSM3004579 TALL_JS_37 (polyA+ RNA)
# SRX3700661 SRR6727659 GSM3004580 TALL_JS_39 (polyA+ RNA)
# SRX3700663 SRR6727661 GSM3004582 TALL_JS_42 (polyA+ RNA)
# SRX3700665 SRR6727663 GSM3004584 TALL_JS_44 (polyA+ RNA)
# SRX3700672 SRR6727670 GSM3004591 TALL_JS_51 (polyA+ RNA)
# SRX3700676 SRR6727674 GSM3004595 TALL_JS_55 (polyA+ RNA)
# SRX3700677 SRR6727675 GSM3004596 TALL_JS_56 (polyA+ RNA)
# SRX3700736 SRR6727734 GSM3004620 TALL_JS_2 (total RNA)
# SRX3700737 SRR6727735 GSM3004621 TALL_JS_3 (total RNA)
# SRX3700738 SRR6727736 GSM3004622 TALL_JS_4 (total RNA)
# SRX3700739 SRR6727737 GSM3004623 TALL_JS_5 (total RNA)
# SRX3700740 SRR6727738 GSM3004624 TALL_JS_10 (total RNA)
# SRX3700741 SRR6727739 GSM3004625 TALL_JS_11 (total RNA)
# SRX3700742 SRR6727740 GSM3004626 TALL_JS_13 (total RNA)
# SRX3700743 SRR6727741 GSM3004627 TALL_JS_14 (total RNA)
# SRX3700744 SRR6727742 GSM3004628 TALL_JS_18 (total RNA)
# SRX3700745 SRR6727743 GSM3004629 TALL_JS_23 (total RNA)
# SRX3700746 SRR6727744 GSM3004630 TALL_JS_26 (total RNA)
# SRX3700747 SRR6727745 GSM3004631 TALL_JS_31 (total RNA)
# SRX3700748 SRR6727746 GSM3004632 TALL_JS_35 (total RNA)
# SRX3700749 SRR6727747 GSM3004633 TALL_JS_37 (total RNA)
# SRX3700750 SRR6727748 GSM3004634 TALL_JS_39 (total RNA)
# SRX3700751 SRR6727749 GSM3004635 TALL_JS_42 (total RNA)
# SRX3700752 SRR6727750 GSM3004636 TALL_JS_44 (total RNA)
# SRX3700753 SRR6727751 GSM3004637 TALL_JS_51 (total RNA)
# SRX3700754 SRR6727752 GSM3004638 TALL_JS_55 (total RNA)
# SRX3700755 SRR6727753 GSM3004639 TALL_JS_56 (total RNA)


## download fastq files
mkdir -p ~/projects/datashare/GSE110637/raw/ 
cd ~/projects/datashare/GSE110637/raw/ 

for srr in SRR6727625 SRR6727626 SRR6727627 SRR6727628 SRR6727633 SRR6727634 SRR6727636 SRR6727637 SRR6727641 SRR6727646 SRR6727649 SRR6727654 SRR6727657 SRR6727658 SRR6727659 SRR6727661 SRR6727663 SRR6727670 SRR6727674 SRR6727675 SRR6727734 SRR6727735 SRR6727736 SRR6727737 SRR6727738 SRR6727739 SRR6727740 SRR6727741 SRR6727742 SRR6727743 SRR6727744 SRR6727745 SRR6727746 SRR6727747 SRR6727748 SRR6727749 SRR6727750 SRR6727751 SRR6727752 SRR6727753
do
 echo ${srr}
 parallel-fastq-dump --threads 12 --tmpdir ./ --gzip --split-files --outdir ./ --sra-id ${srr}
done


# SR or PE?
ls -lha ~/projects/datashare/GSE110637/raw

#metadata
cd ~/projects/datashare/GSE110637/
echo raw/SRR6727625_1.fastq.gz raw/SRR6727625_2.fastq.gz > GSM3004546_notrim_fqgz.info
echo raw/SRR6727626_1.fastq.gz raw/SRR6727626_2.fastq.gz > GSM3004547_notrim_fqgz.info
echo raw/SRR6727627_1.fastq.gz raw/SRR6727627_2.fastq.gz > GSM3004548_notrim_fqgz.info
echo raw/SRR6727628_1.fastq.gz raw/SRR6727628_2.fastq.gz > GSM3004549_notrim_fqgz.info
echo raw/SRR6727633_1.fastq.gz raw/SRR6727633_2.fastq.gz > GSM3004554_notrim_fqgz.info
echo raw/SRR6727634_1.fastq.gz raw/SRR6727634_2.fastq.gz > GSM3004555_notrim_fqgz.info
echo raw/SRR6727636_1.fastq.gz raw/SRR6727636_2.fastq.gz > GSM3004557_notrim_fqgz.info
echo raw/SRR6727637_1.fastq.gz raw/SRR6727637_2.fastq.gz > GSM3004558_notrim_fqgz.info
echo raw/SRR6727641_1.fastq.gz raw/SRR6727641_2.fastq.gz > GSM3004562_notrim_fqgz.info
echo raw/SRR6727646_1.fastq.gz raw/SRR6727646_2.fastq.gz > GSM3004567_notrim_fqgz.info
echo raw/SRR6727649_1.fastq.gz raw/SRR6727649_2.fastq.gz > GSM3004570_notrim_fqgz.info
echo raw/SRR6727654_1.fastq.gz raw/SRR6727654_2.fastq.gz > GSM3004575_notrim_fqgz.info
echo raw/SRR6727657_1.fastq.gz raw/SRR6727657_2.fastq.gz > GSM3004578_notrim_fqgz.info
echo raw/SRR6727658_1.fastq.gz raw/SRR6727658_2.fastq.gz > GSM3004579_notrim_fqgz.info
echo raw/SRR6727659_1.fastq.gz raw/SRR6727659_2.fastq.gz > GSM3004580_notrim_fqgz.info
echo raw/SRR6727661_1.fastq.gz raw/SRR6727661_2.fastq.gz > GSM3004582_notrim_fqgz.info
echo raw/SRR6727663_1.fastq.gz raw/SRR6727663_2.fastq.gz > GSM3004584_notrim_fqgz.info
echo raw/SRR6727670_1.fastq.gz raw/SRR6727670_2.fastq.gz > GSM3004591_notrim_fqgz.info
echo raw/SRR6727674_1.fastq.gz raw/SRR6727674_2.fastq.gz > GSM3004595_notrim_fqgz.info
echo raw/SRR6727675_1.fastq.gz raw/SRR6727675_2.fastq.gz > GSM3004596_notrim_fqgz.info
echo raw/SRR6727734_1.fastq.gz raw/SRR6727734_2.fastq.gz > GSM3004620_notrim_fqgz.info
echo raw/SRR6727735_1.fastq.gz raw/SRR6727735_2.fastq.gz > GSM3004621_notrim_fqgz.info
echo raw/SRR6727736_1.fastq.gz raw/SRR6727736_2.fastq.gz > GSM3004622_notrim_fqgz.info
echo raw/SRR6727737_1.fastq.gz raw/SRR6727737_2.fastq.gz > GSM3004623_notrim_fqgz.info
echo raw/SRR6727738_1.fastq.gz raw/SRR6727738_2.fastq.gz > GSM3004624_notrim_fqgz.info
echo raw/SRR6727739_1.fastq.gz raw/SRR6727739_2.fastq.gz > GSM3004625_notrim_fqgz.info
echo raw/SRR6727740_1.fastq.gz raw/SRR6727740_2.fastq.gz > GSM3004626_notrim_fqgz.info
echo raw/SRR6727741_1.fastq.gz raw/SRR6727741_2.fastq.gz > GSM3004627_notrim_fqgz.info
echo raw/SRR6727742_1.fastq.gz raw/SRR6727742_2.fastq.gz > GSM3004628_notrim_fqgz.info
echo raw/SRR6727743_1.fastq.gz raw/SRR6727743_2.fastq.gz > GSM3004629_notrim_fqgz.info
echo raw/SRR6727744_1.fastq.gz raw/SRR6727744_2.fastq.gz > GSM3004630_notrim_fqgz.info
echo raw/SRR6727745_1.fastq.gz raw/SRR6727745_2.fastq.gz > GSM3004631_notrim_fqgz.info
echo raw/SRR6727746_1.fastq.gz raw/SRR6727746_2.fastq.gz > GSM3004632_notrim_fqgz.info
echo raw/SRR6727747_1.fastq.gz raw/SRR6727747_2.fastq.gz > GSM3004633_notrim_fqgz.info
echo raw/SRR6727748_1.fastq.gz raw/SRR6727748_2.fastq.gz > GSM3004634_notrim_fqgz.info
echo raw/SRR6727749_1.fastq.gz raw/SRR6727749_2.fastq.gz > GSM3004635_notrim_fqgz.info
echo raw/SRR6727750_1.fastq.gz raw/SRR6727750_2.fastq.gz > GSM3004636_notrim_fqgz.info
echo raw/SRR6727751_1.fastq.gz raw/SRR6727751_2.fastq.gz > GSM3004637_notrim_fqgz.info
echo raw/SRR6727752_1.fastq.gz raw/SRR6727752_2.fastq.gz > GSM3004638_notrim_fqgz.info
echo raw/SRR6727753_1.fastq.gz raw/SRR6727753_2.fastq.gz > GSM3004639_notrim_fqgz.info


cat *.info

# qc align count
# put wf on luke and luachn
cd ~/projects/heatshock/results/GSE110637/
snakemake -s ~/projects/heatshock/results/GSE110637/wf.py --cores 24 -pn
snakemake -s ~/projects/heatshock/results/GSE110637/wf.py --cores 49 --cluster "oarsub --project epimed -l nodes=1/core={threads},walltime=6:00:00 " -pn






