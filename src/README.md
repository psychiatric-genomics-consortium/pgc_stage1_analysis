## Howrigan - Stage 1 QC reusable code

**script: family_based_pca.R**
 * Runs PCA on related samples by projecting PCs from an unrelated set to related pairs
 * Uses files generated from the Ricopili pcaer module
 * Generates scatterplot plot of PCs
 * Run `Rscript family_based_pca.R --help` to view all options

Example script with no reference data:
```
Rscript family_based_pca.R \
--data_id [5-letter ID used] \
--pcaer_dir /psych/ripke/stage1_pgc/howrigan/path/to/pcaer_dataID \
--pc_count 10 \
--graph TRUE
```

Example script with reference data and specifying related samples to remove:     
```
Rscript family_based_pca.R \
--data_id data1.1KG \
--refdata TRUE \
--pcaer_dir /psych/ripke/stage1_pgc/howrigan/path/to/pcaer_data1.1KG \
--famex /psych/ripke/stage1_pgc/howrigan/path/to/data1.founders \
--pc_count 20 \
--graph TRUE
```


**script: related_samples_filter.R**
 * Read in .genome and .famex file to get related samples filtered in Ricopili pcaer module
 * Writes breakdown of samples flagged/filtered in standard output and to file `[dataID].related_filtered.txt`
 * Writes out filtered FID/IID to output file `[dataID].related_filtered.ids`

Example command using pcaer module output
```
Rscript related_samples_filter.R \
pcaer_[dataID]/[dataID].mepr.genome \
pcaer_[dataID]/[dataID].famex \
[dataID]
```


**script: PC_filter_casecon.R**
 * Designed for .mds_cov file from Ricopili
 * Quickly apply PC cutoffs and flag PC outliers after running Ricopili pcaer module
 * Creates graphs with filtered samples highlighted
 * Uses a separate file with PC/low/high columns
 * NA values will be ignored
 * Writes out cutoffs and samples flagged to standard output
 * See top of script for more detailed instructions

Example command using pcaer module output
```
Rscript PC_filter_casecon.R \
--data_id abcd1 \
--cutoff_file cutoff.txt \
--pcaer_dir pcaer_abcd1 \
--refdata FALSE \
--graph TRUE
```


**script: PC_filter_related.R**
 * Use after running family_based_pca.R and determining samples to fitler
 * Quickly apply PC cutoffs and flag PC outliers
 * Creates graphs with filtered samples highlighted
 * Uses a separate file with PC/low/high columns
 * NA values will be ignored
 * Writes out cutoffs and samples flagged to standard output
 * See top of script for more detailed instructions

Example command using pcaer module output
```
Rscript PC_filter_related.R \
--data_id abcd1 \
--cutoff_file cutoff.txt \
--pcaer_dir pcaer_abcd1 \
--pca_file abcd1.related_pca.eigenvec \
--refdata FALSE \
--graph TRUE
```


**file: cutoff.txt**
 * Sample cutoff.txt file to copy and use for PC_filter scripts
 * Fill in NAs with cutoff values where desired

Or just copy/paste directly from here:
```
PC   low   high
C1   NA   NA
C2   NA   NA
C3   NA   NA
C4   NA   NA 
C5   NA   NA 
C6   NA   NA  
C7   NA   NA  
C8   NA   NA
C9   NA   NA
C10   NA   NA
```

