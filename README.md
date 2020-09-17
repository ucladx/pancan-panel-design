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

### Germline SNPs

* Ensure frequent heterozygosity across each human subpopulation, to ensure uniform sensitivity of allelic imbalance (AI)
* Usable to detect loss of heterozygosity (LOH), copy-number aberration (CNA), homologous recombination deficiency (HRD)
* SNPs per gene for gene-level AI/CNA, evenly distributed across genome for arm-level CNA, at telomeres to detect telomeric AI
* Usable for QCs steps like fingerprinting, detecting tumor-normal swaps, cross contamination, and identifying contaminant samples

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
