# Dirofilaria immitis populaiton genomics: genome improvement

### stephen doyle

- want to see if we can improve on the recently published assembly 
- will try 
    - recover genes missing from the assembly by comparison of BUSCOs between original Di WBP17 assembly and new assembly
    - scaffold assembly against Brugia malayi, which is properly chromosomal




## Recover missing genes
- there was a difference in the BUSCO scores from the original dirofilaria assembly and the updated assembly
- 

```bash
cd /nfs/users/nfs_s/sd21/lustre_link/dirofilaria_immitis/BUSCO/MISSING_GENES



ln -s ../BUSCO_di_ICBAS_JMDir_1.0_genome_nematoda_odb10/run_nematoda_odb10/missing_busco_list.tsv ICBAS_JMDir_1.0.missing.list

ln -s ../BUSCO_di_wbp17_genome_nematoda_odb10/run_nematoda_odb10/missing_busco_list.tsv di_wbp17.missing.list

ln -s ../BUSCO_di_wbp17_genome_nematoda_odb10/run_nematoda_odb10/

ln -s ../BUSCO_di_wbp17_genome_nematoda_odb10/run_nematoda_odb10/full_table.tsv di_wbp17.full_table.txt

ln -s ../../REFERENCE_GENOMES/dirofilaria_immitis.PRJEB1797.WBPS17.genomic.fa
samtools faidx dirofilaria_immitis.PRJEB1797.WBPS17.genomic.fa

# get list of busco IDs present in di_wbp17 but missing in ragtag
diff -y <(sort ICBAS_JMDir_1.0.missing.list) <(sort di_wbp17.missing.list) | grep "<" | cut -f1 | sort | uniq > di_wbp17.uniq.busco_ids.list

# 30 di_wbp17.uniq.busco_ids.list


# extract contig data per gene, and when there is more than one hit, take the longest contig ID
>tmp2
while read -r BUSCO_ID; do 
    grep "$BUSCO_ID" di_wbp17.full_table.txt | cut -f3 > tmp 
    grep -f tmp dirofilaria_immitis.PRJEB1797.WBPS17.genomic.fa.fai | sort -k2,2nr | head -n1 | cut -f1 >> tmp2
    done < di_wbp17.uniq.busco_ids.list

# using the longest contig ID, extract the contig fasta    
> di_wbp17.uniq.busco_ids.fasta
cat tmp2 | sort | uniq | while read CONTIG_ID; do 
    samtools faidx dirofilaria_immitis.PRJEB1797.WBPS17.genomic.fa ${CONTIG_ID} >> di_wbp17.uniq.busco_ids.fasta; 
    done

grep -c ">" di_wbp17.uniq.busco_ids.fasta
#> 26
#>> 26 contigs/scaffolds contain 30 BUSCO genes missing from the current assembly
```


## Merge the old genome and the recovered contigs
```bash

ln -s ../../REFERENCE_GENOMES/GCA_024305405.1_ICBAS_JMDir_1.0_genomic.fa

cat GCA_024305405.1_ICBAS_JMDir_1.0_genomic.fa di_wbp17.uniq.busco_ids.fasta > dimmitis.merged_genome.fa

assembly-stats dimmitis.merged_genome.fa
# stats for dimmitis.merged_genome.fa
# sum = 93204778, n = 136, ave = 685329.25, largest = 12535412
# N50 = 4253153, n = 7
# N60 = 3601687, n = 10
# N70 = 1986382, n = 13
# N80 = 1563077, n = 19
# N90 = 739389, n = 28
# N100 = 433, n = 136
# N_count = 94816
# Gaps = 468

```
- bit bigger now


## Run BUSCO to see the effect of merging
```bash
# load conda env
conda activate busco_5.4.3

bsub.py --queue long --threads 20 60 busco_di_merged_nematoda_odb10 \
    "busco --in dimmitis.merged_genome.fa --out BUSCO_di_merged_genome_nematoda_odb10 --mode genome --lineage_dataset /nfs/users/nfs_s/sd21/lustre_link/databases/busco/nematoda_odb10 --cpu 20 -f -r"


bsub.py --queue long --threads 20 60 busco_di_merged_eukaryota_odb10 \
    "busco --in dimmitis.merged_genome.fa --out BUSCO_di_merged_genome_eukaryota_odb10 --mode genome --lineage_dataset /nfs/users/nfs_s/sd21/lustre_link/databases/busco/eukaryota_odb10 --cpu 20 -f -r"

```

- nematoda: C:95.0%[S:87.8%,D:7.2%],F:1.8%,M:3.2%,n:3131
- eukaryota: C:95.3%[S:90.2%,D:5.1%],F:2.0%,M:2.7%,n:255

- improvement in BUSCO recovery, with the downside there is increase duplication




## Closing gaps using pacbio reads
```bash
bsub.py --threads 20 30 tgsgapcloser_pb "/nfs/users/nfs_s/sd21/lustre_link/software/GENOME_IMPROVEMENT/TGS-GapCloser/tgsgapcloser --scaff dimmitis.merged_genome.fa --reads ../../GENOME_IMPROVEMENT/PACBIO_DATA/SRR10533235_subreads.fastq --output pb --thread 20 --ne --tgstype pb"


```


## Scaffolding against Bmalayi

```bash 
ln -s ../../REFERENCE_GENOMES/brugia_malayi.PRJNA10729.WBPS17.genomic.fa

bsub.py 1 ragtag "ragtag.py scaffold brugia_malayi.PRJNA10729.WBPS17.genomic.fa pb.scaff_seqs"

assembly-stats ragtag.scaffold.fasta
#stats for ragtag.scaffold.fasta
#sum = 93191081, n = 73, ave = 1276590.15, largest = 27859139
#N50 = 15438497, n = 3
#N60 = 15438497, n = 3
#N70 = 15220724, n = 4
#N80 = 15220724, n = 4
#N90 = 15199470, n = 5
#N100 = 433, n = 73
#N_count = 93785
#Gaps = 338
```

## Run BUSCO to see the effect of merging
```bash
# load conda env
conda activate busco_5.4.3

bsub.py --queue long --threads 20 60 busco_di_ragtag_nematoda_odb10 \
    "busco --in ragtag.scaffold.fasta --out BUSCO_di_ragtag_genome_nematoda_odb10 --mode genome --lineage_dataset /nfs/users/nfs_s/sd21/lustre_link/databases/busco/nematoda_odb10 --cpu 20 -f -r"


bsub.py --queue long --threads 20 60 busco_di_ragtag_eukaryota_odb10 \
    "busco --in ragtag.scaffold.fasta  --out BUSCO_di_ragtag_genome_eukaryota_odb10 --mode genome --lineage_dataset /nfs/users/nfs_s/sd21/lustre_link/databases/busco/eukaryota_odb10 --cpu 20 -f -r"

```

C:95.0%[S:88.2%,D:6.8%],F:1.8%,M:3.2%,n:3131



## Purge duplicates 
- https://github.com/dfguan/purge_dups#pplg
```bash

minimap2 -xmap-pb ragtag.scaffold.fasta /nfs/users/nfs_s/sd21/lustre_link/dirofilaria_immitis/GENOME_IMPROVEMENT/PACBIO_DATA/SRR10533235_subreads.fastq | gzip -c - > paf.gz

/nfs/users/nfs_s/sd21/lustre_link/software/PACBIO/purge_dups/bin/pbcstat paf.gz 
#(produces PB.base.cov and PB.stat files)
/nfs/users/nfs_s/sd21/lustre_link/software/PACBIO/purge_dups/bin/calcuts PB.stat > cutoffs 2>calcults.log


/nfs/users/nfs_s/sd21/lustre_link/software/PACBIO/purge_dups/bin/split_fa ragtag.scaffold.fasta > ragtag.scaffold.split

minimap2 -xasm5 -DP ragtag.scaffold.split ragtag.scaffold.split | gzip -c - > ragtag.scaffold.split.self.paf.gz

/nfs/users/nfs_s/sd21/lustre_link/software/PACBIO/purge_dups/bin/purge_dups -2 -T cutoffs -c PB.base.cov ragtag.scaffold.split.self.paf.gz > dups.bed 2> purge_dups.log

/nfs/users/nfs_s/sd21/lustre_link/software/PACBIO/purge_dups/bin/get_seqs -e dups.bed ragtag.scaffold.fasta 

assembly-stats purged.fa
#stats for purged.fa
#sum = 88533216, n = 29, ave = 3052869.52, largest = 26831765
#N50 = 15418668, n = 3
#N60 = 15418668, n = 3
#N70 = 15199470, n = 4
#N80 = 15199470, n = 4
#N90 = 15131962, n = 5
#N100 = 522, n = 29
#N_count = 33899
#Gaps = 147


```

## Run BUSCO to see the effect of purging
```bash
# load conda env
conda activate busco_5.4.3

bsub.py --queue long --threads 20 60 busco_di_purge_nematoda_odb10 \
    "busco --in purged.fa --out BUSCO_di_purge_genome_nematoda_odb10 --mode genome --lineage_dataset /nfs/users/nfs_s/sd21/lustre_link/databases/busco/nematoda_odb10 --cpu 20 -f -r"


bsub.py --queue long --threads 20 60 busco_di_purge_eukaryota_odb10 \
    "busco --in purged.fa  --out BUSCO_di_purge_genome_eukaryota_odb10 --mode genome --lineage_dataset /nfs/users/nfs_s/sd21/lustre_link/databases/busco/eukaryota_odb10 --cpu 20 -f -r"

```

- nematoda: C:94.4%[S:92.5%,D:1.9%],F:1.6%,M:4.0%,n:3131
- eukaryota: C:95.3%[S:94.5%,D:0.8%],F:2.0%,M:2.7%,n:255

- slightly better overall
    - also higher duplications, but these have been reduced



- where "sequence_name_conversion.txt" contains:
```bash
Bm_013_RagTag_1	dirofilaria_immitis_scaffold_0024
Bm_v4_Chr1_scaffold_001_RagTag_1	dirofilaria_immitis_chr1
Bm_v4_Chr2_contig_001_RagTag_1	dirofilaria_immitis_chr2
Bm_v4_Chr3_scaffold_001_RagTag_1	dirofilaria_immitis_chr3
Bm_v4_Chr4_scaffold_001_RagTag_1	dirofilaria_immitis_chr4
Bm_v4_ChrX_scaffold_001_RagTag_1	dirofilaria_immitis_chrX
JAKNDB010000003.1_1	dirofilaria_immitis_scaffold_0001
JAKNDB010000005.1_1	dirofilaria_immitis_scaffold_0002
JAKNDB010000010.1_1	dirofilaria_immitis_scaffold_0003
JAKNDB010000011.1_1	dirofilaria_immitis_scaffold_0004
JAKNDB010000012.1_1	dirofilaria_immitis_scaffold_0005
JAKNDB010000016.1_1	dirofilaria_immitis_scaffold_0006
JAKNDB010000017.1_1	dirofilaria_immitis_scaffold_0007
JAKNDB010000021.1_1	dirofilaria_immitis_scaffold_0008
JAKNDB010000022.1_1	dirofilaria_immitis_scaffold_0009
JAKNDB010000023.1_1	dirofilaria_immitis_scaffold_0010
JAKNDB010000029.1_1	dirofilaria_immitis_scaffold_0011
JAKNDB010000030.1_1	dirofilaria_immitis_scaffold_0012
JAKNDB010000034.1_1	dirofilaria_immitis_scaffold_0013
JAKNDB010000035.1_1	dirofilaria_immitis_scaffold_0014
JAKNDB010000040.1_1	dirofilaria_immitis_scaffold_0015
JAKNDB010000048.1_1	dirofilaria_immitis_scaffold_0016
JAKNDB010000049.1_1	dirofilaria_immitis_scaffold_0017
JAKNDB010000057.1_1	dirofilaria_immitis_scaffold_0018
JAKNDB010000058.1_1	dirofilaria_immitis_scaffold_0019
JAKNDB010000062.1_1	dirofilaria_immitis_scaffold_0020
JAKNDB010000071.1_1	dirofilaria_immitis_scaffold_0021
JAKNDB010000091.1_1	dirofilaria_immitis_scaffold_0022
nDi.2.2.scaf00004_1	dirofilaria_immitis_scaffold_0023
```

## rename scaffolds
```bash
cp purged.fa purged.renamed.fa
while read OLD NEW; do 
    sed -i "s/${OLD}/${NEW}/g" purged.renamed.fa; 
    done <  sequence_name_conversion.txt