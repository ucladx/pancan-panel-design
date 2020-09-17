## Shortlist targets for a hybridization capture-based cancer DNA sequencing panel

### Overview of required targets:
* Cancer genes: Coding exons of all major isoforms of ~1k genes
* Germline SNPs: Commonly heterozygous ~50k loci to measure allelic imbalance and copy-number
* Microsatellites: Frequently instable loci in tumor samples with MSI per PCR testing
* Non-coding loci: UTRs, introns, promoters, etc. with known variants or breakpoints
* Viral sequences: With known relevance to cancer diagnosis/prognosis/treatment

### Cancer genes

Gather, rename, compare, review, and shortlist genes from 16 sources:
* fmi_cdx_324 - [FoundationOne CDx panel](https://info.foundationmedicine.com/hubfs/FMI%20Labels/FoundationOne_CDx_Label_Technical_Info.pdf)
* fmi_heme_407 - [FoundationOne Heme panel](https://assets.ctfassets.net/w98cd481qyp0/42r1cTE8VR4137CaHrsaen/baf91080cb3d78a52ada10c6358fa130/FoundationOne_Heme_Technical_Specifications.pdf)
* tcga_pancan_299 - [TCGA PancanAtlas cancer drivers](https://www.cell.com/cms/10.1016/j.cell.2018.02.060/attachment/cf6b14b1-6af1-46c3-a2c2-91008c78e87f/mmc1.xlsx)
* oncokids_170 - [CHLA OncoKids panel](https://ars.els-cdn.com/content/image/1-s2.0-S1525157818301028-mmc2.docx)
* goal_core_497 - [GOAL Consortium's core genes](https://github.com/GoalConsortium/goal_misc/blob/4d1190b/GOAL_Core_497_genes_sorted.bed)
* resistance_77 - From reviews of resistance to [immunotherapy](https://pubmed.ncbi.nlm.nih.gov/28187290) and [PARPi](https://pubmed.ncbi.nlm.nih.gov/32364757)
* tempus_xt_648 - [Tempus xT panel](https://www.tempus.com/wp-content/uploads/2020/05/xTGene-List_110419.pdf)
* pgdx_elio_505 - [PGDx elio tissue complete panel](https://www.personalgenome.com/assets/resources/PGDx_ETC_-_Gene_Tables_IVD-2_final.pdf)
* msk_impact_505 - [MSK-IMPACT panel](https://www.oncokb.org/cancerGenes)
* trusight_523 - [Illumina TruSight Oncology 500 panel](https://www.illumina.com/content/dam/illumina-marketing/documents/products/datasheets/trusight-oncology-500-and-ht-data-sheet-1170-2018-010.pdf)
* profile_447 - [Profile/OncoPanel version 3 panel](https://researchcores.partners.org/data/wiki_pages/97/POPv3_TEST_INFORMATION.pdf)
* mayo_neuro_148 - [Mayo Clinic's neuro-oncology panel](https://www.mayocliniclabs.com/it-mmfiles/Targeted_DNA_Gene_Regions_Interrogated_by_Neuro-Oncology_Panel.pdf)
* chop_solid_309 - [CHOP solid tumor panel](https://www.testmenu.com/chop/Tests/785964)
* intogen_568 - [Cancer drivers from intogen.org](https://www.intogen.org/download?file=IntOGen-Drivers-20200201.zip)
* chop_heme_211 - [CHOP heme malignancy panel](https://www.testmenu.com/chop/Tests/786447)
* mayo_heme_42 - [Mayo Clinic's OncoHeme panel](https://www.mayocliniclabs.com/it-mmfiles/Targeted_Genes_Interrogated_by_OncoHeme_Next-Generation_Sequencing.pdf)

Shortlisted 1075 genes for a new panel as follows, and labeled them as "Y" under column "Selected" in [`data/cancer_genes_review.txt`](data/cancer_genes_review.txt):
1. Include all genes from: fmi_cdx_324, fmi_heme_407, tcga_pancan_299, oncokids_170, goal_core_497, and resistance_77.
2. Include genes seen in at least 2 of the 16 lists.
3. Include genes seen in at least 1 of the 16 lists, with >26 mutations per Mbp.
4. Include genes seen in at least 1 of the 16 lists, with >3.15 mean gain (CN>=2).
5. Include genes seen in at least 1 of the 16 lists, with >1.5 and <1.82 mean loss (CN<=2). Mean loss <1.5 was on sex chromosomes.
6. Exclude genes PDE4DIP and FAT3 to reduce DNA-seq costs. They are each ~20kbp and of unknown relevance to cancer.
7. Add more cancer genes using [OncoKB](https://www.oncokb.org/cancerGenes) summaries, mut/CN burdens, and literature review (`ARID3A, ARID4B, ATP6AP1, ATP6V1B2, ATXN7, CRBN, CYP19A1, DEK, DKK4, EGR1, ERF, EZH1, EZHIP, FLI1, GAB1, GAB2, KBTBD4, KLF3, KNSTRN, LCK, LRP5, LRP6, LTB, MEF2D, MIDEAS, MIR142, NADK, PGBD5, PPP4R2, PRKD1, PTP4A1, PTPN1, RAC2, ROBO1, RRAS, SAMHD1, SERPINB3, SERPINB4, SESN2, SESN3, SETD1B, SETDB2, SMYD3, SP140, SPRTN, STAG1, STAT1, STAT2, STK19, TCL1B, TET3, WIF1, WWTR1`).

**Note:** Mutations per Mbp calculated using TCGA+TARGET MuTect2 MAFs from NCI GDC, and gene sizes from Gencode v35, and mean gain/loss calculated using TCGA+TARGET Gistic2 gene-level absolute CN from NCI GDC.

Extract the gene names and their Ensembl ENSG IDs:
```bash
cut -f1,2,9 data/cancer_genes_review.txt | grep -w Y$ | cut -f1,2 | sort > data/exon_targets_gene_list.txt
```

Create a BED file for these genes' exons with 2bp flanks, using their Gencode basic isoforms except level 3 (not verified nor curated):
```bash
gzip -dc /mdl/gencode/gencode.v35.basic.annotation.gff3.gz | grep -w "$(cut -f2 data/exon_targets_gene_list.txt)" | perl -a -F'\t' -ne '%t=map{split("=")} split(";",$F[8]); if(($t{gene_type} eq "protein_coding" and $F[2] eq "CDS" and $t{level} ne "3" and $t{ID}!~m/PAR/) or ($t{gene_type}=~/lncRNA|miRNA|pseudogene/ and $F[2] eq "exon")){$F[3]-=3; $F[4]+=2; print join("\t",@F[0,3,4],$t{gene_name},@F[5,6])."\n"}' | sort -s -k1,1V -k2,2n -k3,3n | bedtools merge -i - -c 4 -o distinct > data/exon_targets_grch38.bed
```

**Source:** ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/gencode.v35.basic.annotation.gff3.gz

Try different combinations of isoforms/exons/annotations in steps above and measure total target space:
```bash
awk -F"\t" '{sum+=$3-$2} END {print sum}' data/exon_targets_grch38.bed
```

**Note:** 2780910 bps in total; 2967992 if we used 8bp flanks; 5985251 if we included all UTRs.

### Germline SNPs

* Ensure frequent heterozygosity across each human subpopulation, to ensure uniform sensitivity of allelic imbalance (AI)
* Usable to detect loss of heterozygosity (LOH), copy-number aberration (CNA), homologous recombination deficiency (HRD)
* SNPs per gene for gene-level AI/CNA, evenly distributed across genome for arm-level CNA, at telomeres to detect telomeric AI
* Usable for QCs steps like fingerprinting, detecting tumor-normal swaps, cross contamination, and identifying contaminant samples

gnomAD v3 lists GRCh38 variants from 71702 genomes (and no exomes). For each variant, they list AC, AN, and nhomalt per subpopulation. Frequency of heterozygosity per subpopulation can be calculated as: `(AC - 2*nhomalt) / (AN/2)`

Wrote a script and used it to shortlist 50k SNPs with a frequency of heterozygotes >10% in all 9 gnomAD v3 subpopulations:
```bash
python3 bin/select_gnomad_snps.py --gnomad-vcf /mdl/gnomad/gnomad.genomes.r3.0.sites.vcf.bgz --max-snps 50000 --output-bed data/snp_targets_grch38.bed
```

**Source:** https://storage.googleapis.com/gnomad-public/release/3.0/vcf/genomes/gnomad.genomes.r3.0.sites.vcf.bgz

### Microsatellites

* GOAL consortium targets ACTC, BAT-25, BAT-26, BAT-34C4, BAT-40, D10S197, D11S925, D13S175, D17S250, D18S35, D18S474, D18S55, D18S58, D18S64, D18S69, D2S123, D3S1478, D3S1766, D5S107, D5S346, D6S260, HSPH1_T17, MONO-27, NR-21, NR-24, NR-27, PENTA-C, PENTA-D
* There should be sufficiently more microsatellites in coding regions to detect MSI with more sensitivity

Fetch loci of upstream/downstream probes around microsatellites targeted by GOAL consortium:
```bash
curl -sL https://github.com/GoalConsortium/goal_misc/raw/e9966b5/GOAL_GRCh38%2Bviral/Consortium_Probes_All_Final.probes_GRCh38%2Bviral.bed | cut -f-2 -d\| | grep -w MSI | sed -E 's/MSI\|//' > data/goal_msi_targets_grch38.bed
```

**Note:** Missing downstream probe for MONO-27, NR-24, D18S58 and upstream probe for D10S197, D17S250 - but these are short enough to be captured by 1 probe each

### Non-coding loci

* GOAL consortium targets fusion breakpoints for ALK, BRAF, CD74, EGFR, ETV6, FGFR1, FGFR2, FGFR3, MET, NTRK1, NTRK2, PAX8, RAF1, RET, ROS1, RSPO3, TFE3, TFEB

Fetch loci of GOAL consortium probes targeting fusion breakpoints and merge probes <=60bp (half a probe) apart:
```bash
curl -sL https://github.com/GoalConsortium/goal_misc/raw/e9966b5/GOAL_GRCh38%2Bviral/Consortium_Probes_All_Final.probes_GRCh38%2Bviral.bed | cut -f1 -d\| | grep _Fusion | bedtools merge -i - -d 60 -c 4 -o distinct > data/goal_fusion_targets_grch38.bed
```

* Promoter regions of TERT and APC
* Introns and UTRs of cancer genes that contain pathogenic variants per ClinVar and literature review
* Breakpoints of MSH2 inversion (PMID: 18335504, 12203789, 24114314)
* Breakpoints of PMS2 insertion, including those homologous to PMS2CL pseudogene (PMID: 29792936)
* Breakpoints of 40Kbp duplication between upstream promoter of GREM1 and 3' end of SCG5 (PMID: 22561515, 26493165)

### Viral sequences

* GOAL consortium targets BKV, EBV1, EBV2, HHV8, HPV11, HPV16, HPV18, HPV31, HPV33, HPV35, HPV52, HPV56, HPV58, HPV59, HPV6b, HPyV7, HTLV1, HTLV2, JCV, MCPyV

Fetch viral sequences targeted by GOAL consortium, bgzip, and index:
```bash
curl -sL https://github.com/GoalConsortium/goal_misc/raw/e9966b5/GOAL_GRCh38%2Bviral/goal_viral_genomes.fasta -o data/goal_viral_targets.fa
bgzip -@4 data/goal_viral_targets.fa
samtools faidx data/goal_viral_targets.fa.gz
```

* Detect presence/absence of viral DNA - integration site detection or quantification not necessary
* We must dilute virus probes to avoid under-capture of human targets in tumors with virus content
