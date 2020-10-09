# PGC stage 1 analysis workflow

# Table of Contents
* [Stage 1 analysis overview](#stage-1-analysis-overview)
* [Resources](#resources)
* [Stage 1 workflow requirements](#stage-1-workflow-requirements)
* [Case-control workflow](#case-control-workflow)
* [Family-based workflow](#family-based-workflow)
* [Trio / family data considerations](#Trio-family-data-considerations)


## Stage 1 analysis overview

Main objectives
 * Match phenotype data to genotype data
 * Provide transparent documentation of QC steps
 * Remove poor quality SNPs
 * Remove poor quality samples
 * Cases and controls are matched for broad/continental ancestry
 * If multiple case/control clusters, dataset is split into separate continental ancestry subgroups
 * If population based sample, related individuals are flagged and filtered appropriately
 * If family-based sample, family relationships confirmed and unconfirmed families filtered appropriately

Common challenges
 * Dataset is too small to run meaningful QC
 * File format is inconsistent with PLINK expectations
 * Unusual ID characters stalling analysis
 * Determining useful PC cutoff filters 
 * Mixed ancestry in sample - do you filter to single continental ancestry or split into ancestry subgroups?

Deliverables for each dataset
 * QCed PLINK files (.bed / .bim / .fam files) uploaded to specified QC directory
 * Completed QC report (google doc templates below) containing:
   * Sample counts before and after QC
   * Description of analyses done
   * Screenshots of Ricopili QC table, Post-QC Manhattan/QQ plot, and relevant PC scatterplots
   * Samples removed due to cryptic relatedness, pedigree error, and population stratification
 * [Case-control QC document template](https://docs.google.com/document/d/1zWY_EX7yWH2ZQi2-ob5rppGAlrd3jH9yXVDhI3wZPNM/edit?usp=sharing)
 * [Family-based QC document template](https://docs.google.com/document/d/1N_WoZyUJcJQuEFLqlp60nYbO23tc3Rc58yabUasD7G4/edit?usp=sharing)

Imputation
 * Currently, imputation is run at the discretion of the Stage 2 analyst
 * Stage 2 analysts can instruct Stage 1 analysts to run imputation, so long as they provide the location of the preferred reference panel and commands to run imputation on the server


## Resources

[Active list of PGC reps/analysts](https://docs.google.com/spreadsheets/d/1CpHDhW0RYVTxBvo9GzrsS38HU_62cR3GZ_3KknBumR4/edit#gid=311537287)
 * Google spreadsheet of PGC DRC reps, stage 1, and stage 2 analysts
 * PGC DRC Google Folder: https://drive.google.com/drive/folders/1xd67r2Qj6GUeIinDCiRHfT1D--aWgmrc

[Stage 1 analyst assignments](https://docs.google.com/spreadsheets/d/1rfl8o7Wm8qAOsoTzy7E9uCA1s2EY-aOu7EXZovFpems/edit?usp=sharing)
 * Running list of assigned datasets
 * Confirmation of data uploads
 * NOTE: Reach out to Don Hucks if your assigned dataset record needs to be updated

[Uploading data to LISA server](https://docs.google.com/document/d/1JYptCBxtcv9cLXhWCzydf_rqd_5xECfmtjdlxO12PEY/edit?usp=sharing)
 * Instructions for PGC collaborators to upload genetic data to LISA server

[PGC DRC directory structure](https://docs.google.com/document/d/1SSBpc9dGnUAGm1zeFyKf_QZ7MgVkDMLNetoALuUq-1o/edit?usp=sharing)
 * Contains list of stage 1 analyst emails and logins
 * Outlines DRC directory access and structure
 * Current directory make-up:
```
/home
└── pgcdac
    └── DWFV2CJb8Piv_0116_pgc_data
        └── pgcdrc
            └── mdd
                ├── incoming_datasets
                │   └── exam1
                │   └── exam2
                ├── working
                │   └── wave2
                │       └── exam1
                └── stage2_ready
                    └── wave2
                        └── exam1
```

[Ricopili website](https://sites.google.com/a/broadinstitute.org/ricopili/home)
 * Primary website for downloading and using Ricopili
 * Read through the entire site before running Ricopili
 * Google docs of additional analyses are available in the HowTo section (left sidebar)

[ext folder](https://github.com/psychiatric-genomics-consortium/pgc_stage1_analysis/ext)
 * External files
 * Ricopili_BioArxiv_2019.pdf
 * Sample Picopili pipeline workflow
 * Sample Ricopili pipeline workflow

[src folder](https://github.com/psychiatric-genomics-consortium/pgc_stage1_analysis/src)
 * Source scripts
 * Portable R scripts to process Ricopili output
 * See README in folder for more details


## Stage 1 workflow requirements

Prior to running through the workflow, the stage 1 analyst should have information on:
 * **Dataset ID**
   * TIP: This should be either provided by the data intake form or data provider
 * **Genotype platform**
   * Often this will be available in the Genotyping platform (GP) delivery bucket or in the [Data intake inquiry form](https://docs.google.com/spreadsheets/d/1StQeoldBPkV0uxpHoS7wWxR20ir6cIABKhbY985lc9s/edit?usp=sharing)
 * **Continental ancestry** (or mixed ancestry)
   * Common ancestry abbreviations: eur / afr / asn / amr
 * **Number of cases**
 * **Number of controls**
 * **Number of females**
 * **Number of males**
 * **Ascertainment** (unrelated case-control, trios, extended pedigrees, or mixed)

You can derive counts from the phenotype files received, but do not proceed until you have all this information.


# Case-control workflow:

**_Create a QC report document in Google docs_**
 * Make a copy of this [Case-control QC document template](https://docs.google.com/document/d/1zWY_EX7yWH2ZQi2-ob5rppGAlrd3jH9yXVDhI3wZPNM/edit?usp=sharing)
 * Rename to particular dataset starting with 5-letter identifier (e.g. cloz1_scz_PGC_S1_QC_report)
 * Provide sample breakdown of incoming dataset
 * Write ANY notes that make this dataset unique (e.g. "This dataset has both unrelated and family-based samples" or "This is the 4th wave of data from the XXX cohort")

**_QC prior to Ricopili_**
 * Removing samples for phenotypic reasons (unknown phenotype, splitting by phenotypes, etc..)
 * If there are multiple phenotype designations, note how designations were collapsed to case/control labels
 * Updating FID/IID due to typos
 * **NOTE: This set of PLINK files is your Ricopili-ready dataset, and iterations of preimp-QC/PCA will occur on this dataset**

**_Run Ricopili pre-imputation QC on dataset_**
 * Create qc subdirectory (e.g. qc1) and link current PLINK files
 * Run preimp-qc module
 * Keep track of the QC run in QC_CYCLE of .names file (most datasets require at least 2 runs)
 * Attach screenshot of Flags, General Info, Manhattan-Plot post-QC from qc[run].pdf to QC report

**_Run Ricopili PCA on dataset without reference data_**
 * Create pca subdirectory (e.g. pca1/) and soft-link to PLINK QC files
 * Run pcaer command with --prefercase (so controls are preferentially removed in relatedness check)
 * **In the first PCA run**, Find out which samples are removed with relatedness checking (1st/2nd col. in /pcaer_[data]/[data].mepr.famex file)
 * Map these onto the genetic relatedness file to see how related they are (/pcaer_[data]/[data].mepr.genome file)
 * Report the number of samples removed and list of individual IDs (if not too long) due to cryptic relatedness in QC report 
 * **TIP: Read the .famex and .genome file into [related_samples_filter.R](https://github.com/psychiatric-genomics-consortium/pgc_stage1_analysis/blob/main/src/related_samples_filter.R) to get quick output of these results**
 * Manually check first 10 PCs for multiple clusters or unusual clusters, particularly case/control biased clusters
 * Attach screenshots of relevant PCA plots (at minimum PCA1/PCA2 graph) to QC report

**_Run Ricopili PCA with reference data_**
 * Create pca subdirectory (e.g. refpca1/) and soft-link to PLINK QC files AND reference files
 * Run pcaer command with --prefercase and include dataset and reference .bim file
 * **NOTE:** Ignore the .famex file in this case. While different individuals maybe filtered in this .famex file due to relatedness, this should have affect the final list of samples filtered
 * Manually check first 10 PCs for multiple clusters or unusual clusters, particularly case/control biased clusters
 * Attach screenshots of relevant PCA plots (at minimum PCA1/PCA2 graph) to QC report

**_Determine whether dataset should be split into multiple ancestries and/or PCA outliers filtered out_**
 * Set cutoffs to filter out PCA outliers in PC axis. Some may require both low/high cutoffs to remove outliers, some require only one.
 * List number of samples removed at each cutoff in QC document
 * If splitting into mutliple ancestries, label samples with imputed continental ancestry
 * **TIP: Input cutoffs into a small file and run [PC_filter_casecon.R](https://github.com/psychiatric-genomics-consortium/pgc_stage1_analysis/blob/main/src/PC_filter_casecon.R) to get list of samples filtered and improved graphs. Keep re-running to refine cutoffs quickly too!**

**_Filtering samples and re-running QC and PCA_**
 * Track or retain FID/IID of cryptic related samples and PCA outliers
 * Remove/split samples from original PLINK file and re-run preimp-QC/PCA
 * Continue as necessary until case/control or families no longer have PCA outliers
 * **NOTE: Through each iteration of QC/PCA, keep an ID list of samples removed due to cryptic relatedness/PCA, and remove them from the original PLINK files in the next iteration. In this respect, samples removed for genotyping rate, sex check, and Fhet are not removed in advance, but removed each time preimp-QC is run.**

**_Notify project manager / stage 2 analyst / data group leader_**
 * Fill out remainder of QC checklist
 * Fill in a highlight the Post-QC sample breakdown in the QC document
 * Create .pdf of both QC document and S1 checklist and keep a copy in your github subdirectory
 * Notify project manager / stage 2 analyst / PI about QCed dataset
 * Add to your own analysis log that the QCed dataset is complete


# Family-based workflow:

**_Create a QC report document in Google docs_**
 * Make a copy of this [Family-based QC document template](https://docs.google.com/document/d/1N_WoZyUJcJQuEFLqlp60nYbO23tc3Rc58yabUasD7G4/edit?usp=sharing)
 * Rename to particular dataset starting with 5-letter identifier (e.g. cloz1_scz_PGC_S1_QC_report)
 * Provide sample breakdown of incoming dataset
 * Write ANY notes that make this dataset unique (e.g. "This dataset has both unrelated and family-based samples" or "This is the 4th wave of data from the XXX cohort")

**_QC prior to Ricopili_**
 * Removing samples for phenotypic reasons (unknown phenotype, incomplete family, splitting by phenotypes, etc..)
 * If there are multiple phenotype designations, note how designations were collapsed to case/control labels
 * Updating FID/IID due to typos

**_PRE-PEDIGREE CONFIRMATION: Run Ricopili pre-imputation QC on dataset_**
 * Create qc subdirectory (e.g. qc0) and link current PLINK files
 * Keep track of the QC run in QC_CYCLE of .names file
 * Run preimp-qc module with `--trio` flag
 * Attach screenshot of Flags, General Info, Family size breakdown, and Manhattan-Plot post-QC from qc[run].pdf to QC report

**_PRE-PEDIGREE CONFIRMATION: Run Ricopili PCA with reference data_**
 * Create pca subdirectory (e.g. refpca0/) and soft-link to PLINK QC files AND reference files
 * Run pcaer command (NOTE: do not use --prefer-case as it may list excessive individuals in .famex file)
 * Run additional PCA projection in PLINK to project PCs for related samples based on unrelated set
 * **TIP: Running [family_based_pca.R](https://github.com/psychiatric-genomics-consortium/pgc_stage1_analysis/blob/main/src/family_based_pca.R) with `--refdata TRUE` will run PCA projection for related samples and create a PDF with 2-dimensional PCA scatterplots** 
 * Manually check first 6 PCs for distinct or unusual clusters
 * Attach screenshots of relevant PCA plots (at minimum PCA1/PCA2 graph) to QC report

**_PRE-PEDIGREE CONFIRMATION: Determine whether dataset should be split into multiple ancestries and/or PCA outliers filtered out_**
 * Set cutoffs to filter out PCA outliers in PC axis. Some may require both low/high cutoffs to remove outliers, some require only one.
 * List number of samples removed at each cutoff in QC document
 * If splitting into mutliple ancestries, label samples with imputed continental ancestry
 * **TIP: Input cutoffs into a small file and run [PC_filter_related.R](https://github.com/psychiatric-genomics-consortium/pgc_stage1_analysis/blob/main/src/PC_filter_related.R) to get list of samples filtered and improved graphs. Keep re-running to refine cutoffs quickly too!**

**_Confirming pedigrees using genotype data_**
 * Confirm parents/proband/sibling relationships using Identity-By-Descent (IBD) estimates from genotype data
 * Perform strict QC, LD-pruning, and running the [--genome](https://www.cog-genomics.org/plink/1.9/ibd) command in PLINK
 * Flag/filter out samples that:
   * Are unrelated despite having the same Family ID (nonFID)
   * Are unable to confirm parental status (wrong parents)
   * Likely parents that have not been marked as such (unmarked parents)
   * Show cryptic relatedness (crossFID)
 * The filtering process should retain as many individuals when flagging pedigree errors (e.g. only remove one sample from a pair with cryptic relatedness)
 * **TIP: Running [check_ped.sh](https://github.com/Nealelab/pgc_stage1_analysis/blob/master/nbaya/scripts/bash/check_ped.sh) for family-based data to comfirm relatedness will do this all for you**
 * Filter out PCA outliers and related samples on the pre-QC PLINK files with updated FIDs created by Ricopili in the qc0/ subdirectory
 * **NOTE: You'll know you're working with the proper PLINK files when the .bed and .bim file are soft-linked, but not the .fam file (as FIDs were altered)**
 * **This set of samples will mark the beginning of preimp-QC/PCA iterations for this dataset**


**_Run Ricopili pre-imputation QC on dataset_**
 * Create qc subdirectory (e.g. qc1) and link to PLINK files with PCA outliers/unconfirmed relatedness samples removed
 * Keep track of the QC run in QC_CYCLE of .names file (most datasets require at least 2 runs)
 * Run preimp-qc module with `--trio` flag
 * Attach screenshot of Flags, General Info, Family size breakdown, and Manhattan-Plot post-QC from qc[run].pdf to QC report

**_Run Ricopili PCA on dataset without reference data_**
 * Create pca subdirectory (e.g. pca1/) and soft-link to PLINK QC files
 * Run pcaer command (NOTE: do not use --prefer-case as it may list excessive individuals in .famex file)
 * Run additional PCA projection in PLINK to project PCs for related samples based on unrelated set
 * **TIP: Running [family_based_pca.R](https://github.com/psychiatric-genomics-consortium/pgc_stage1_analysis/blob/main/src/family_based_pca.R) will run PCA projection for related samples and create a PDF with 2-dimensional PCA scatterplots** 
 * Manually check first 10 PCs for multiple clusters or unusual clusters, particularly case/control biased clusters
 * Attach screenshots of relevant PCA plots (at minimum PCA1/PCA2 graph) to QC report

**_Run Ricopili PCA with reference data_**
 * Create pca subdirectory (e.g. refpca1/) and soft-link to PLINK QC files AND reference files
 * Run pcaer command (NOTE: do not use --prefer-case as it may list excessive individuals in .famex file)
 * Run additional PCA projection in PLINK to project PCs for related samples based on unrelated set
 * **TIP: Running [family_based_pca.R](https://github.com/psychiatric-genomics-consortium/pgc_stage1_analysis/blob/main/src/family_based_pca.R) with `--refdata TRUE` will run PCA projection for related samples and create a PDF with 2-dimensional PCA scatterplots** 
 * Manually check first 10 PCs for multiple clusters or unusual clusters, particularly case/control biased clusters
 * Attach screenshots of relevant PCA plots (at minimum PCA1/PCA2 graph) to QC report

**_Determine PCA outliers to filtered out_**
 * Set cutoffs to filter out PCA outliers in PC axis. Some may require both low/high cutoffs to remove outliers, some require only one.
 * List number of samples removed at each cutoff in QC document
 * **TIP: Input cutoffs into a small file and run [PC_filter_related.R](https://github.com/psychiatric-genomics-consortium/pgc_stage1_analysis/blob/main/src/PC_filter_related.R) to get list of samples filtered and improved graphs. Keep re-running to refine cutoffs quickly too!**

**_Filtering samples and re-running QC and PCA_**
 * Track or retain FID/IID of cryptic related samples and PCA outliers
 * Remove/split samples from original PLINK file and re-run preimp-QC/PCA
 * Continue as necessary until case/control or families no longer have PCA outliers
 * **NOTE: Through each iteration of QC/PCA, keep an ID list of samples removed due to cryptic relatedness/PCA, and remove them from the original PLINK files in the next iteration. In this respect, samples removed for genotyping rate, sex check, and Fhet are not removed in advance, but removed each time preimp-QC is run.**

**_Notify project manager / stage 2 analyst / data group leader_**
 * Fill out remainder of QC checklist
 * Fill in a highlight the Post-QC sample breakdown in the QC document
 * Create .pdf of both QC document and S1 checklist and keep a copy in your github subdirectory
 * Notify project manager / stage 2 analyst / PI about QCed dataset
 * Add to your own analysis log that the QCed dataset is complete



## Trio / family data considerations

It is crucial that you know whether the dataset is unrelated case/conntrol or family-based, as QC will differ dramtically depending on study design.

Among Family-based studies, there are two main strategies for GWAS meta-analysis:
 1. Convert trios into case/pseudo-control
 2. Maintain family structure and use [Picopili](https://github.com/Nealelab/picopili)
NOTE: Both strategies are determined by the Stage 2 analyst, and trios will not be converted to case/pseudo-controls until after the imputation step 

**GOAL**: Confirm pedigree format, 

Confirming pedigree format
 * Parental ID columns are correctly filled in (MID / PID in .fam file)
 * Phenotypic sex matches parental status (proband can only have one biological mother and one biological father)
 * Table of family structures confirmed / created
   * Count of family sizes (duos / trios / quads / extended families along with unrelated case/control)
   * This should be reported in QC document prior to QC/PCA

Confirming relationships using identity-by-descent estimation
 * Ricopili QC will filter out offspring that show excessive mendelian errors with a parent
 * Running the [--genome](https://www.cog-genomics.org/plink/1.9/ibd) command on an LD-pruned SNP set 

Running Ricopili with trio data
 * Run Ricopili pre-imputation QC with `--trio` flag
   * Adds mendelian error checking to the QC process
   * Creates table of family sizes
   * Runs initial case/pseudo-control GWAS excluding complete trios (often unreliable if predominantly family-based data)   
 * Run Ricopili PCA with parents
   * Running PCA with `--trio` flag will restrict to offspring
   * Using parents will often increase the power to generate ancestry-informative PCAs
   * Can impute the offspring's PCs as the midpoint of the parents





