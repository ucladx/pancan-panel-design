## Shortlist targets for a hybridization capture-based cancer DNA sequencing panel

### Overview of required targets:
* Cancer genes: Coding exons of all major transcripts of ~1k genes
* Microsatellites: Frequently instable loci in tumor samples with MSI per PCR testing
* Non-coding loci: UTRs, introns, promoters, etc. with cancer-associated variants or breakpoints
* Germline SNPs: Commonly heterozygous ~16k loci to measure allelic imbalance and copy-number

### Cancer genes

Gather, rename, compare, review, and shortlist genes from 16 sources:
* fmi_cdx_324 - [FoundationOne CDx panel](https://info.foundationmedicine.com/hubfs/FMI%20Labels/FoundationOne_CDx_Label_Technical_Info.pdf)
* fmi_heme_407 - [FoundationOne Heme panel](https://www.foundationmedicine.be/content/dam/rfm/sample-reports/f1heme/Technical%20Specs%20F1Heme%202018_06_digital.pdf)
* tcga_pancan_299 - [TCGA PancanAtlas cancer drivers](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6029450/bin/NIHMS948705-supplement-8.xlsx)
* oncokids_170 - [CHLA OncoKids panel](https://ars.els-cdn.com/content/image/1-s2.0-S1525157818301028-mmc2.docx)
* goal_core_497 - [GOAL Consortium's core genes](https://github.com/GoalConsortium/goal_misc/blob/4d1190b/GOAL_Core_497_genes_sorted.bed)
* resistance_77 - From reviews of resistance to [immunotherapy](https://pubmed.ncbi.nlm.nih.gov/28187290) and [PARPi](https://pubmed.ncbi.nlm.nih.gov/32364757)
* tempus_xt_648 - [Tempus xT panel](https://www.tempus.com/wp-content/uploads/2019/11/xTGene-List_110419.pdf)
* pgdx_elio_505 - [PGDx elio tissue complete panel](https://www.personalgenome.com/assets/resources/PGDx_ETC_-_Gene_Tables_IVD-2_final.pdf)
* msk_impact_505 - [MSK-IMPACT panel](https://www.oncokb.org/cancerGenes)
* trusight_523 - [Illumina TruSight Oncology 500 panel](https://www.illumina.com/content/dam/illumina/gcs/assembled-assets/marketing-literature/trusight-oncology-500-data-sheet-1170-2018-010/trusight-oncology-500-and-ht-data-sheet-1170-2018-010.pdf)
* profile_447 - [Profile/OncoPanel version 3 panel](https://researchcores.partners.org/data/wiki_pages/97/POPv3_TEST_INFORMATION.pdf)
* mayo_neuro_148 - [Mayo Clinic's neuro-oncology panel](https://www.mayocliniclabs.com/it-mmfiles/Targeted_DNA_Gene_Regions_Interrogated_by_Neuro-Oncology_Panel.pdf)
* chop_solid_309 - [CHOP solid tumor panel](https://www.testmenu.com/chop/Tests/785964)
* intogen_568 - [Cancer drivers from intogen.org](https://www.intogen.org/download?file=IntOGen-Drivers-20200201.zip)
* chop_heme_211 - [CHOP heme malignancy panel](https://www.testmenu.com/chop/Tests/786447)
* mayo_heme_42 - [Mayo Clinic's OncoHeme panel](https://www.mayocliniclabs.com/it-mmfiles/Targeted_Genes_Interrogated_by_OncoHeme_Next-Generation_Sequencing.pdf)

Shortlist 1080 genes for a new panel as follows, and label them as "Y" under column "Selected" in [`data/cancer_genes_review.txt`](data/cancer_genes_review.txt):
1. Include all genes from: fmi_cdx_324, fmi_heme_407, tcga_pancan_299, oncokids_170, goal_core_497, and resistance_77.
2. Include genes seen in at least 2 of the 16 lists.
3. Include genes seen in at least 1 of the 16 lists, with >26 mutations per Mbp.
4. Include genes seen in at least 1 of the 16 lists, with >3.15 mean gain (CN>=2).
5. Include genes seen in at least 1 of the 16 lists, with >1.5 and <1.82 mean loss (CN<=2). Mean loss <1.5 was on sex chromosomes.
6. Exclude genes PDE4DIP and FAT3 to reduce DNA-seq costs. They are each ~20kbp and of unknown relevance to cancer.
7. Add [OncoKB genes](https://www.oncokb.org/cancerGenes) based on their descriptions, mut/CN burdens, and literature review (`ARID3A, ARID4B, ATP6AP1, ATP6V1B2, ATXN7, CRBN, CYP19A1, DEK, DKK4, EGR1, ERF, EZH1, EZHIP, FLI1, GAB1, GAB2, KBTBD4, KLF3, KNSTRN, LCK, LRP5, LRP6, LTB, MEF2D, MIDEAS, MIR142, NADK, PGBD5, PPP4R2, PRKD1, PTP4A1, PTPN1, RAC2, ROBO1, RRAS, SAMHD1, SERPINB3, SERPINB4, SESN2, SESN3, SETD1B, SETDB2, SMYD3, SP140, SPRTN, STAG1, STAT1, STAT2, STK19, TCL1B, TET3, WIF1, WWTR1`).
8. Add more cancer genes requested by UCLA colleagues.
    - `INO80` - DNA repair, regulates abundance and positioning of nucleosomes, mutated in 4% of DLBCL, algorithmically predicted to be a cancer driver per intogen.org
    - `LUC7L2` - Splicing factor, low expression in 14% of MDS patients, del(7q) is common and truncating mutations have been reported
    - `MBD4` - Germline loss predisposes to Uveal Melanoma and Leukemia, targeted by DFCI's OncoPanel (aka Profile)
    - `SRP72` - Germline loss predisposes to Familial Aplasia and Myelodysplasia, targeted by Mayo's OncoHeme panel
    - `FRK` - Somatic mutations in 6% of hepatocellular adenomas, and frequent rearrangements

Mutations per Mbp calculated using TCGA+TARGET MuTect2 MAFs from NCI GDC, and gene sizes from Gencode v35. Mean gain/loss calculated using TCGA+TARGET Gistic2 gene-level absolute CN from NCI GDC.

Extract the gene names and their Ensembl ENSG IDs:
```bash
cut -f1,2,9 data/cancer_genes_review.txt | grep -w Y$ | cut -f1,2 | sort > data/exon_targets_gene_list.txt
```

Create a BED file for these genes' coding regions with 2bp flanks, using their Gencode basic transcripts except level 3 (not verified nor curated):
```bash
gzip -dc /hot/tracks/gencode.v38.basic.annotation.gff3.gz | grep -w "$(cut -f2 data/exon_targets_gene_list.txt)" | perl -a -F'\t' -ne '%t=map{split("=")} split(";",$F[8]); if(($t{gene_type} eq "protein_coding" and $F[2] eq "CDS" and $t{level} ne "3" and $t{ID}!~m/PAR/) or ($t{gene_type}=~/lncRNA|miRNA|pseudogene/ and $F[2] eq "exon")){$F[3]-=3; $F[4]+=2; print join("\t",@F[0,3,4],$t{gene_name}.":".$F[2],@F[5,6])."\n"}' | sort -s -k1,1V -k2,2n -k3,3n | bedtools merge -i - -c 4 -o distinct > targets/exon_targets_grch38.bed
```

**Source:** ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.basic.annotation.gff3.gz

Try different combinations of transcripts/exons/annotations in steps above and measure total target space:
```bash
awk -F"\t" '{sum+=$3-$2} END {print sum}' targets/exon_targets_grch38.bed
```

2798232 bps in total; 2805929 if we include level 3 transcripts; 2986850 if we used 8bp flanks; 6117858 if we included all UTRs.

### Microsatellites

* GOAL consortium targets 28 microsatellites including the 7 sites targeted by the Promega PCR kit
* There will be hundreds more frequently mutable microsatellites in other regions we are targeting

Fetch loci of upstream/downstream probes around microsatellites targeted by GOAL consortium:
```bash
curl -sL https://github.com/GoalConsortium/goal_misc/raw/e9966b5/GOAL_GRCh38%2Bviral/Consortium_Probes_All_Final.probes_GRCh38%2Bviral.bed | cut -f-2 -d\| | grep -w MSI | sed -E 's/MSI\|//; s/$/:Microsatellite/' > targets/goal_msi_targets_grch38.bed
```

Missing downstream probe for MONO-27, NR-24, D18S58 and upstream probe for D10S197, D17S250. But these are short enough to be captured by 1 probe each.

### Non-coding loci

* GOAL consortium targets fusion breakpoints for `ALK, BRAF, CD74, EGFR, ETV6, FGFR1, FGFR2, FGFR3, MET, NTRK1, NTRK2, PAX8, RAF1, RET, ROS1, RSPO3, TFE3, TFEB`

Fetch loci of GOAL consortium probes targeting fusion breakpoints and merge probes <=60bp (half a probe) apart:
```bash
curl -sL https://github.com/GoalConsortium/goal_misc/raw/e9966b5/GOAL_GRCh38%2Bviral/Consortium_Probes_All_Final.probes_GRCh38%2Bviral.bed | cut -f1 -d\| | grep _Fusion | sed -E 's/_Fusion/:FusionSite/' | bedtools merge -i - -d 60 -c 4 -o distinct > targets/goal_fusion_targets_grch38.bed
```

For each gene in the panel, find associated [GeneHancer clusters](https://genome.ucsc.edu/cgi-bin/hgTrackUi?g=geneHancer#TRACK_HTML) with score>325, and then fetch loci of overlapping [ENCODE cCREs](https://genome.ucsc.edu/cgi-bin/hgTrackUi?g=encodeCcreCombined#TRACK_HTML) with score>450:
```bash
curl -sL 'https://api.genome.ucsc.edu/getData/track?genome=hg38;track=geneHancerInteractions' | jq -r '.geneHancerInteractions[] | [.geneHancerChrom,.geneHancerStart,.geneHancerEnd,.name,.score,.geneStrand] | @tsv' | perl -a -F'\t' -ne 'BEGIN{%gs=map{chomp; ($_,1)}`cut -f1 data/exon_targets_gene_list.txt`} $F[1]--; ($g)=split("/",$F[3]); print join("\t",@F) if($F[4]>325 && $gs{$g})' | sort -s -k1,1V -k2,2n -k3,3n > data/genehancer_regions_grch38.bed
curl -sL 'https://api.genome.ucsc.edu/getData/track?genome=hg38;track=encodeCcreCombined' | jq -r '.encodeCcreCombined[] | [.chrom,.chromStart,.chromEnd,.ucscLabel,.name,.score,.strand] | @tsv' | perl -a -F'\t' -ne '$F[1]--; print join("\t",@F[0..3,5,6])' | sort -s -k1,1V -k2,2n -k3,3n | bedtools intersect -f 1 -wo -a - -b data/genehancer_regions_grch38.bed | perl -a -F'\t' -ne '($g)=split("/",$F[9]); print join("\t",@F[0..2],"$g:$F[3]",@F[4,5])."\n" if($F[4]>450)' > data/encode_ccre_grch38.bed
```

The ENCODE cCREs cover ~800Kbp which is too large. Reduce this to ~5Kbp by targeting only cCREs associated with APC, FOXA1, PMS2, PTEN, and TERT:
```bash
perl -a -F'\t' -pe '($g,$t)=split(":",$F[3]); $_="" unless($g=~m/^(APC|FOXA1|PMS2|PTEN|TERT)$/)' data/encode_ccre_grch38.bed > targets/non_coding_targets_grch38.bed
```

For each gene in the panel, target the ends of 5'UTRs of MANE transcripts where mutations could cause loss of function:
```bash
gzip -dc /hot/tracks/gencode.v38.basic.annotation.gff3.gz | grep -w "$(cut -f2 data/exon_targets_gene_list.txt)" | perl -a -F'\t' -ne '%t=map{split("=")} split(";",$F[8]); if($t{gene_type} eq "protein_coding" and $F[2] eq "five_prime_UTR" and $t{tag}=~/MANE/ and $t{ID}!~m/PAR/){print join("\t",$F[0],$F[3]-2,($F[3]+118<$F[4]?$F[3]+118:$F[4]),$t{gene_name}.":5pUTR")."\n".join("\t",$F[0],($F[4]-118>$F[3]?$F[4]-118:$F[3]),$F[4]+2,$t{gene_name}.":5pUTR")."\n"}' | sort -s -k1,1V -k2,2n -k3,3n | bedtools merge -i - -c 4 -o distinct >> targets/non_coding_targets_grch38.bed
sort -su -k1,1V -k2,2n -k3,3n targets/non_coding_targets_grch38.bed -o targets/non_coding_targets_grch38.bed
```

Fetch loci of Pathogenic and Likely Pathogenic (P/LP) mutations with decent evidence from ClinVar related to cancer:
```bash
curl -sL 'https://api.genome.ucsc.edu/getData/track?genome=hg38;track=clinvarMain' | jq -r '.clinvarMain[] | [.chrom,.chromStart,.chromEnd,._variantId,._clinSignCode,.reviewStatus,.phenotypeList] | @tsv' | perl -a -F'\t' -ne 'print if($F[4]=~m/^P|LP$/ && $F[5]=~m/practice guideline|expert panel|multiple submitters/i) && $F[6]=~m/cancer|lynch|neoplas|tumor|adenoma|carcinoma|li-fraumeni|polyposis|hippel-lindau/i' | cut -f1-4 | sort -su -k1,1V -k2,2n -k3,3n > data/clinvar_plp_muts_grch38.bed
```

Target the subset of P/LP ClinVar mutations that do not overlap targeted exons:
```bash
bedtools subtract -a data/clinvar_plp_muts_grch38.bed -b targets/exon_targets_grch38.bed | sed 's/$/:ClinVar/' | sort -su -k1,1V -k2,2n -k3,3n >> targets/non_coding_targets_grch38.bed
sort -s -k1,1V -k2,2n -k3,3n targets/non_coding_targets_grch38.bed -o targets/non_coding_targets_grch38.bed
```

Manually added the following into `targets/non_coding_targets_grch38.bed`:
* Breakpoints of MSH2 inversion (PMID: 18335504, 12203789, 24114314)
* Breakpoints of PMS2 retrotransposon insertion and intronic regions homologous to PMS2CL pseudogene (PMID: 29792936)
* Breakpoints of 40Kbp duplication between GREM1 and SCG5 (PMID: 22561515, 26493165)
* Pathogenic germline variants near cancer susceptibility genes (ClinVar, [April 2020](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/archive/variant_summary_2020-04.txt.gz))

Downloaded TSV-format catalog data from ebi.ac.uk/gwas for these genome-wide association studies listing GRCh38 genomic loci per trait-associated SNP:
* 265 SNP loci for Prostate cancer risk prediction (PMID: 33398198)
* 461 SNP loci for Pan-cancer risk prediction (PMID: 32887889)

Manually converted into `targets/gwas_targets_grch38.bed` with 711 unique loci. Data provenance not possible, and parsing out loci was messy and needed regexes. So, cannot document steps for reproducibility.

### Germline SNPs

* Ensure frequent heterozygosity across each human subpopulation, to ensure uniform sensitivity of allelic imbalance (AI)
* Usable to detect loss of heterozygosity (LOH), copy-number aberration (CNA), homologous recombination deficiency (HRD)
* SNPs per gene for gene-level AI/CNA, evenly distributed across genome for arm-level CNA, at telomeres to detect telomeric AI
* Usable for QCs steps like fingerprinting, detecting tumor-normal swaps, cross contamination, and identifying contaminant samples

gnomAD v3 lists GRCh38 variants from 71702 genomes (and no exomes). For each variant, they list AC, AN, and nhomalt per subpopulation. Frequency of heterozygosity per subpopulation can be calculated as: `(AC - 2*nhomalt) / (AN/2)`

Use a script to find ~2m candidate SNPs with a frequency of heterozygotes >25% in all 9 gnomAD v3 subpopulations:
```bash
python3 scripts/select_gnomad_snps.py --gnomad-vcf /hot/ref/gnomad/gnomad.genomes.r3.0.sites.vcf.bgz --max-snps 2100000 --output-bed data/snp_candidates_grch38.bed
sort -su -k1,1V -k2,2n -k3,3n data/snp_candidates_grch38.bed -o data/snp_candidates_grch38.bed
```

**Source:** https://storage.googleapis.com/gcp-public-data--gnomad/release/3.0/vcf/genomes/gnomad.genomes.r3.0.sites.vcf.bgz

We captured and sequenced ~42k of these SNPs in the v1 design using 16 FFPE samples and analyzed them. Per testing, major contributors to off-bait capture are Alu repeats.

Exclude SNPs <90bp from an Alu repeat, or within genomic loci flagged by vendor as contributors to off-bait capture:
```bash
curl -sL https://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/rmsk.txt.gz | gzip -dc | perl -a -F'\t' -ne 'print join("\t",@F[5..7],"$F[12]:$F[10]",$F[2],$F[9])."\n"' > data/repeat_masker_grch38.bed
grep Alu: data/repeat_masker_grch38.bed | grep -P "^chr(\d+|X|Y)\t" | bedtools slop -b 90 -g /hot/ref/GRCh38_Verily_v1.genome.fa.fai -i | cut -f1-3 - data/vendor_flagged_probes.bed data/vendor_flagged_targets.bed | sort -su -k1,1V -k2,2n -k3,3n | bedtools subtract -a data/snp_candidates_grch38.bed -b - > data/snp_candidates_good_grch38.bed
```

Using per-target Picard HsMetrics across 16 samples (MQ>=1), choose SNPs with >0.5 and <1.5 normalized depth in >=12 samples:
```bash
grep -hv ^chrom /hot/bams/mdl_cancer_v1.2/*.mq1.hsmetrics_by_target | cut -f1-3,8 | awk -F'\t' '{if($4>0.5 && $4<1.5){c[$1"\t"$2-1"\t"$3]++}} END{for(v in c){if(c[v]>=12){print v}}}' | sort -su -k1,1V -k2,2n -k3,3n | bedtools window -w 60 -a - -b data/snp_candidates_good_grch38.bed | cut -f4- | sort -su -k1,1V -k2,2n -k3,3n > data/snp_candidates_good_grch38.bed_
mv -f data/snp_candidates_good_grch38.bed_ data/snp_candidates_good_grch38.bed
```

Target 2474 SNPs within exon targets or <240bp near them, which gives 664 genes at least 1 SNP likely to help detect LOH (don't choose >1 SNP <200bp apart):
```bash
bedtools window -w 240 -a targets/exon_targets_grch38.bed -b data/snp_candidates_good_grch38.bed | cut -f5-8 | sort -su -k1,1V -k2,2n -k3,3n | bedtools spacing -i - | awk -F'\t' '{if($5>=200) print}' | cut -f1-4 | sed 's/$/:SNP_LOH/' > targets/snp_targets_grch38.bed
```

Add 13526 more SNPs (16k total, sufficient for HRD) that are most distant from their nearest SNPs:
```bash
grep SNP_LOH$ targets/snp_targets_grch38.bed | bedtools slop -b 200 -g /hot/ref/GRCh38_Verily_v1.genome.fa.fai -i | bedtools subtract -a data/snp_candidates_good_grch38.bed -b - | bedtools spacing -i - | sort -k7,7rn | head -n13526 | cut -f1-4 | sed 's/$/:SNP_CNV/' >> targets/snp_targets_grch38.bed
sort -s -k1,1V -k2,2n -k3,3n targets/snp_targets_grch38.bed -o targets/snp_targets_grch38.bed
```

### Probe design

Combined all targets into a single merged BED file:
```bash
cat targets/*_targets_grch38.bed | sort -s -k1,1V -k2,2n -k3,3n > data/ucla_mdl_targets_grch38.bed
```

Estimated how many probes will be needed for 1x tiling (targets <=120bp get one 120bp probe each, others get total bps ÷ 120):
```bash
bedtools merge -i data/ucla_mdl_targets_grch38.bed -c 4 -o distinct | awk -F"\t" '{len=$3-$2; sum+=(len<120?120:len)} END {print sum/120}'
```

At this point, we sent our targets to the custom panel vendor in BED format. They ran bioinformatics tools to design the tiling and content of 120bp probes across our targets. They also flagged tricky targets that are most likely to cause off-bait capture (usually homology with common genomic repeats), and sent us `data/vendor_flagged_targets.bed`. Over 8k of these were intergenic SNPs that we removed from the design, since the remaining 42k were sufficient for HRD/LOH detection. We reviewed the remaining non-SNP targets based on importance.

Review non-SNP targets with partial or no coverage by vendor's probe design:
```bash
echo -e "Region\tLabels\tLength\tSkipped_Length\tFraction_Skipped\tReason_to_Keep" > data/ucla_mdl_tricky_targets_grch38.txt
bedtools intersect -wo -a data/ucla_mdl_targets_grch38.bed -b data/vendor_flagged_targets.bed | perl -ane '$l=$F[2]-$F[1]; $s=$F[6]-$F[5]; print join("\t","$F[0]:$F[1]-$F[2]",$F[3],$l,$s,$s/$l,"")."\n" unless($F[3]=~m/^(rs\d+|\.)$/)' >> data/ucla_mdl_tricky_targets_grch38.txt
```

Shortlisted 32 important targets in `data/ucla_mdl_tricky_targets_grch38.txt` we asked vendor to capture anyway, at the cost of some off-bait reads. Final target/bait loci of manufactured probes from vendor are stored at `ucla_mdl_cancer_ngs_v1_baits.grch38.bed` and `ucla_mdl_cancer_ngs_v1_targets.grch38.bed`. These are useful for calculating hybrid selection metrics. Note that bait loci are merged and precise tiling/genomic loci of each 120bp bait is not indicated. Most vendors consider this proprietary, but will share this with their clients under an NDA.

### Testing

Captured 2 pools of 8 FFPE samples each and sequenced on a HiSeq 2500 Rapid at 2x100bp. Ran GATK's best-practice for secondary analysis of TN-pairs using [Sarek v2.7.1](https://nf-co.re/sarek/2.7.1). The 16 BAMs after Picard MarkDuplicates and GATK BQSR are stored at `/hot/bams/mdl_cancer_v1/*.bam`. Per Picard HsMetrics, overall capture efficiency (`FOLD_ENRICHMENT`) was poor due to abundance of off-bait reads, mostly from intergenic SNPs. On vendor's recommendation, we tried it again with post-hyb wash temperature increased from 68deg to 70deg. `FOLD_ENRICHMENT` roughly doubled, though could be better. BAMs for these are stored at `/hot/bams/mdl_cancer_v1.2/*.bam`. These were sequenced on a NovaSeq 6000 SP at 2x150bp and shared with vendor for further analysis. They ran scripts that find the closest matching probe per off-bait read, and then ranked all probes by percent-contribution to the total number of off-bait reads in `data/ucla_mdl_cancer_ngs_v1_ranked_probes.bed`.

Create Picard-friendly interval lists for baits/targets, and generate HsMetrics per target requiring MQ>=1:
```bash
picard BedToIntervalList --INPUT ucla_mdl_cancer_ngs_v1_baits.grch38.bed --OUTPUT ucla_mdl_cancer_ngs_v1_baits.grch38.ilist --SEQUENCE_DICTIONARY /hot/bams/mdl_cancer_v1.2/NGSPanel_FFPE001.md.bam
picard BedToIntervalList --INPUT ucla_mdl_cancer_ngs_v1_targets.grch38.bed --OUTPUT ucla_mdl_cancer_ngs_v1_targets.grch38.ilist --SEQUENCE_DICTIONARY /hot/bams/mdl_cancer_v1.2/NGSPanel_FFPE001.md.bam
find /hot/bams/mdl_cancer_v1.2 -name "*.bam" | parallel -j16 picard CollectHsMetrics --INPUT {} --OUTPUT {.}.mq1.hsmetrics --PER_TARGET_COVERAGE {.}.mq1.hsmetrics_by_target --TARGET_INTERVALS ucla_mdl_cancer_ngs_v1_targets.grch38.ilist --BAIT_INTERVALS ucla_mdl_cancer_ngs_v1_baits.grch38.ilist --REFERENCE_SEQUENCE /hot/ref/GRCh38_Verily_v1.genome.fa --INCLUDE_INDELS --MINIMUM_MAPPING_QUALITY 1
```

Find the type of repeats that contribute the most to off-bait reads:
```bash
bedtools intersect -wao -f 0.25 -a data/ucla_mdl_cancer_ngs_v1_ranked_probes.bed -b data/repeat_masker_grch38.bed | awk -F"\t" 'OFS="\t" {sum[$8]+=$4} END {for (i in sum) print sum[i],i}' | sort -k1,1rg | less
```
77% of off-bait reads are from targets overlapping `Alu` repeats of which 45% are `AluS` repeats, 23% are `AluY`, and 7% are `AluJ`. Another 12% of off-bait reads are from simple 2-mer repeats like `(TG)n`, `(AC)n`, `(CA)n`, and `(GT)n`.

Create a BED file of the 4479 probes that contributed >0.001% of off-bait reads:
```bash
awk '{if($4>0.00001){print}}' data/ucla_mdl_cancer_ngs_v1_ranked_probes.bed > data/vendor_flagged_probes.bed
```
