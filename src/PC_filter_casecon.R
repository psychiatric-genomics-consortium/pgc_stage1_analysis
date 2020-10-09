## USAGE:
# Rscript PC_filter_casecon.R --data_id abcd1 --cutoff_file cutoff.txt

## Full example (with default parameters):
# Rscript PC_filter_casecon.R \
# --data_id abcd1 \
# --cutoff_file cutoff.txt \
# --pcaer_dir pcaer_abcd1 \
# --refdata FALSE \
# --graph TRUE

# Example cutoff.txt file:
# PC   low    high
# C1   0.01   0.1
# C2   -0.2   0.5
# C5   NA     0.4

## NOTES: 
## Script will output filtered IDs and zipped pdf in whatever directory the script was run
## Put in full file path for --pcaer_dir if running from a separate directory than the normal pcaer output directory
## Cutoff file must adhere to this format (header as PC/low/high, PC column as C1, C4, etc...)
## Can tolerate missing PCs and PCs that have NA in both low and high row 
## To create graphs without making any cutoffs, do not include the --cutoff_file flag
## To just get filtered samples from cutoffs, use --graph FALSE
## If --ref_data TRUE, looks for specific 4 population 1KG samples used in PGC stage1 pipeline


#!/bin/Rscript

## Load getopt
if(require("getopt")){
    print("getopt is loaded correctly")
} else {
    print("trying to install getopt")
    install.packages("getopt")
    if(require(getopt)){
        print("getopt installed and loaded")
    } else {
        stop("could not install getopt")
    }
}

#we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
  'verbose', 'v', 2, "integer",
  'help'   , 'h', 0, "logical",
  'data_id'  , 'i', 1, "character",
  'pcaer_dir'   , 'd', 1, "character",
  'refdata'         , 'r', 1, "logical",
  'cutoff_file'     , 'c', 1, "character",  
  'graph'     , 'g', 2, "logical"  
), byrow=TRUE, ncol=4)
opt = getopt(spec)

# if help was asked for print a friendly message 
# and exit with a non-zero error code
if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=TRUE))
  print(spec)
  q(status=1)
}

#set some reasonable defaults for the options that are needed,
#but were not specified.
cutoffs <- TRUE
if ( is.null(opt$verbose ) ) { opt$verbose = TRUE }
if ( is.null(opt$data_id    ) ) { opt$data_id    = 'dataset'     }
if ( is.null(opt$pcaer_dir      ) ) { opt$pcaer_dir      = paste0('pcaer_',opt$data_id)     }
if ( is.null(opt$refdata ) ) { opt$refdata = FALSE }
if ( is.null(opt$cutoff_file ) ) { cutoffs <- FALSE }
if ( is.null(opt$cutoff_file ) ) { opt$cutoff_file = ''}
if ( is.null(opt$graph   ) ) { opt$graph   = TRUE    }

#print verbose output unless FALSE
if ( opt$verbose ) { 
    write("verbose output...",stdout()) 
    print(opt)
    cat(paste0('pcaer directory/data_id is = ',opt$pcaer_dir,"/",opt$data_id),'\n')
    cat(paste0('cutoffs=',cutoffs,' and file is ',opt$cutoff_file),'\n')    
}



## read in and apply cutoffs
if (cutoffs == TRUE) { 

  ## read in file:
  cf <- read.table(opt$cutoff_file,h=T,stringsAsFactors=F)

  if (nrow(cf) == 0 & opt$verbose) { 
    cat('unable to read in cutoff file, please check input \n')
    q(status=1) 
  } ## END if

  if (colnames(cf)[1] != 'PC' & opt$verbose) { 
    cat('cutoff file column header should be 1)PC 2)low 3)high \n')
    q(status=1) 
  } ## END if

  if (ncol(cf) != 3 & opt$verbose) { 
    cat('cutoff file should have three columns: 1)PC 2)low 3)high \n')
    q(status=1) 
  } ## END if

  ## read in unrelated .fam file
  fam <- read.table(paste0(opt$pcaer_dir,'/',opt$data_id,'.menv.fam'),h=F,stringsAsFactors=F)

  if (nrow(fam) == 0 & opt$verbose) { 
    cat('unable to read in fam file, please check input \n')
    q(status=1) 
  } ## END if

  ## read in PC file
  pc <- read.table(paste0(opt$pcaer_dir,'/',opt$data_id,'.menv.mds'),h=T,stringsAsFactors=F)

  if (nrow(pc) == 0 & opt$verbose) { 
    cat('unable to read in pc file, please check input \n')
    q(status=1) 
  } ## END if

  ## get parameters of PC file
  pc_names <- names(pc)[4:ncol(pc)]

  ## remove reference pops
  chk <- strsplit(pc$FID,'_'); id <- sapply(chk,'[',2)
  pc <- subset(pc,id!='pop')
  ## get Ricopili FID cas/con
  tmp <- strsplit(pc$FID,'_',fixed=T)
  pc$AFF <- sapply(tmp,'[',1)

  ## get list of PCs to apply cutoffs
  cf2 <- subset(cf,!is.na(cf$low) | !is.na(cf$high))
  pc_list <- cf2$PC

  ## LooP through each PC that has cutoffs
  for (i in 1:length(pc_list)) {
  
  cutoff_pc <- cf2$PC[i] 
  cutoff_low <- cf2$low[i]
  cutoff_high <- cf2$high[i] 
  pc_in <- eval(parse(text=paste0('pc$',cutoff_pc)))
  ## getting IDs outside the cutoffs 
  low_indx <- pc_in < cutoff_low
  if (is.na(cutoff_low)) { low_indx <- rep(FALSE,length(pc_in)) } 
  high_indx <- pc_in > cutoff_high
  if (is.na(cutoff_high)) { high_indx <- rep(FALSE,length(pc_in)) } 
  indx <- low_indx==T | high_indx==T

  if (i==1) { removed_ids <- pc$IID[indx] }
  if (i >1) { removed_ids <- c(removed_ids,pc$IID[indx]) }

  ## Print out results
  cat('\nPC = ',cutoff_pc,'\n')
  cat('low-end cut = ',cutoff_low,'\n')
  cat('samples flagged = ',sum(low_indx),'\n')
  cat('high-end cut = ',cutoff_high,'\n')
  cat('samples flagged = ',sum(high_indx),'\n')
 
  } ## END of i LooP


## ---- Print out final results 

cat('\n',length(removed_ids),'total samples flagged\n')
 
## table breakdown 
cat('\nTimes flagged:\n')
print(table(table(removed_ids)))

## Samples removed
removed_ids_2 <- unique(removed_ids)
cat('\n',length(removed_ids_2),'unique samples flagged as PC outliers\n')

pc2 <- subset(pc,pc$IID %in% removed_ids_2)

## case/con breakdown 
if (sum(pc2$AFF=='cas') > 0 & sum(pc2$AFF=='con') > 0) {
cat('\nFlagged case/control comparison:\n')
cat(sum(pc2$AFF=='cas'),'cases flagged out of',sum(pc$AFF=='cas'),'(',round(sum(pc2$AFF=='cas')/sum(pc$AFF=='cas'),2),')\n')
cat(sum(pc2$AFF=='con'),'controls flagged out of',sum(pc$AFF=='con'),'(',round(sum(pc2$AFF=='con')/sum(pc$AFF=='con'),2),')\n')
cat('prop test p-value:',prop.test(c(sum(pc2$AFF=='cas'),sum(pc2$AFF=='con')),c(sum(pc$AFF=='cas'),sum(pc$AFF=='con')))$p.value,'\n')
} ## END if

## write out samples to FILTER
write.table(pc2[,c('FID','IID')],paste0(opt$data_id,'.flagged.ids'),col=T,row=F,quo=F,sep='\t')

cat('\nPC-flagged samples written out to:',paste0(opt$data_id,'.flagged.ids'),'\n\n\n')


} ## END of cutoffs file evaluation





## ---- GRAPH of pca (no reference data)


if (opt$graph==TRUE & opt$refdata==FALSE) {

cat('\nprinting graph of PCs\n')

## read in unrelated .fam file
fam <- read.table(paste0(opt$pcaer_dir,'/',opt$data_id,'.menv.fam'),h=F,stringsAsFactors=F)

if (nrow(fam) == 0 & opt$verbose) { 
  cat('unable to read in fam file, please check input \n')
  q(status=1) 
} ## END if

## read in PC file
pc <- read.table(paste0(opt$pcaer_dir,'/',opt$data_id,'.menv.mds'),h=T,stringsAsFactors=F)

if (nrow(pc) == 0 & opt$verbose) { 
  cat('unable to read in pc file, please check input \n')
  q(status=1) 
} ## END if

if (cutoffs == TRUE) { 
  ## read in file:
  cf <- read.table(opt$cutoff_file,h=T,stringsAsFactors=F)
}


## get parameters of PC file
pc_names <- names(pc)[4:ncol(pc)]

## find out which PCs to create lines
if (cutoffs == TRUE) { cf_lines <- pc_names %in% cf$PC }

## merge with menv fam file
pca <- merge(pc,fam,by.x='IID',by.y='V2')

pc_xaxis <- pc_names[seq(1,length(pc_names),2)]
pc_yaxis <- pc_names[seq(2,length(pc_names),2)]

aff <- pca$V6
chk <- strsplit(pca$FID,'_'); id <- sapply(chk,'[',2)
aff[id=='pop'] <- 'pop'
aff[id!='pop' & aff==0] <- 'mis'
aff[id!='pop' & aff==-9] <- 'mis'
aff[aff==2] <- 'cas'
aff[aff==1] <- 'con'

legend_id <- substr(pc$FID,start=1,stop=17)
## update fam to cas/con
grp_id <- substr(legend_id,start=1,stop=3)
grp_id[grp_id=='fam' & aff=='cas'] <- 'cas'
grp_id[grp_id=='fam' & aff=='con'] <- 'con'
end_id <- substr(legend_id,start=4,stop=nchar(legend_id))
legend_id <- paste0(grp_id,end_id)

## create label / color table
ltable <- table(legend_id)
if (length(ltable) < 2) { lcols <- c('red') }
if (length(ltable)==2) { lcols <- c('red','blue') }
if (length(ltable) > 2) { add <- length(ltable)-2; lcols <- c('red','blue',rep('black',add)) }
ltable <- cbind.data.frame(ltable,lcols)
names(ltable) <- c('label','count','colors')

cols <- rep('red',nrow(pc))
cols[aff=='cas'] <- 'red'
cols[aff=='con'] <- 'blue'


## -------- PLOT
pdf(paste0(opt$data_id,'.caseControl_pca.2d.pdf'),height=8,width=10)

par(xpd=NA,mar=c(5, 4, 4, 12) + 0.1)

## Loop through PCs
for (i in 1:length(pc_xaxis)) {

  eval(parse(text=paste0('pcx <- pca$',pc_xaxis[i])))
  eval(parse(text=paste0('pcy <- pca$',pc_yaxis[i])))

  ## set up legend coordinates
  xleg <- range(pcx)[2] + 0.05*(range(pcx)[2] - range(pcx)[1])
  yleg <- range(pcy)[2] - 0.05*(range(pcy)[2] - range(pcy)[1])

  plot(pcx,pcy,xlab=pc_xaxis[i],ylab=pc_yaxis[i],col=cols,cex=0.5,cex.axis=0.7)

  ## grid lines
  abline(v=seq(round(range(pcx)[1],2),round(range(pcx)[2],2),0.01),col='lightgrey',xpd=F)
  abline(h=seq(round(range(pcy)[1],2),round(range(pcy)[2],2),0.01),col='lightgrey',xpd=F)

  ## add case/control again
  points(pcx[aff=='con'],pcy[aff=='con'],col=cols[aff=='con'],cex=0.5)
  points(pcx[aff=='cas'],pcy[aff=='cas'],col=cols[aff=='cas'],cex=0.5)

  ## copy original legend table
  ltable2 <- ltable

  ## Get x-axis cutoffs
  if (cutoffs == TRUE) {
    if ((pc_xaxis[i] %in% cf$PC) == TRUE) {
    pc_row <- which(cf$PC %in% pc_xaxis[i])
    pc_row_low <- cf$low[pc_row]; pc_row_high <- cf$high[pc_row]
    pc_filter <- 0
    ## DRAW LINE and highlight points
    if (!is.na(pc_row_low)) {
      abline(v=pc_row_low,col='darkred',lwd=1.5,xpd=FALSE)
      points(pcx[pcx < pc_row_low],pcy[pcx < pc_row_low],col='palegreen',pch=1,cex=1) 
      pc_filter <- pc_filter + sum(pcx < pc_row_low) }
    if (!is.na(pc_row_high)) { 
      abline(v=pc_row_high,col='darkred',lwd=1.5,xpd=FALSE)
      points(pcx[pcx > pc_row_high],pcy[pcx > pc_row_high],col='palegreen',pch=1,cex=1) 
      pc_filter <- pc_filter + sum(pcx > pc_row_high) }
    ## Add to legend
    xcount <- c(ltable2$count,pc_filter)
    xlabel <- c(as.character(ltable2$label),paste0(cf$PC[pc_row],' filtered = ',pc_filter))
    xcolors <- c(as.character(ltable2$colors),'palegreen')
    ltable2 <- cbind.data.frame(xlabel,xcount,xcolors)
    names(ltable2) <- c('label','count','colors')
  }} ## END if

  ## Get y-axis cutoffs
  if (cutoffs == TRUE) {
    if ((pc_yaxis[i] %in% cf$PC) == TRUE) {
    pc_row <- which(cf$PC %in% pc_yaxis[i])
    pc_row_low <- cf$low[pc_row]; pc_row_high <- cf$high[pc_row]
    pc_filter <- 0
    ## DRAW LINE and highlight points
    if (!is.na(pc_row_low)) { 
      abline(h=pc_row_low,col='darkred',lwd=1.5,xpd=FALSE)
      points(pcx[pcy < pc_row_low],pcy[pcy < pc_row_low],col='chartreuse3',pch=1,cex=1)
      pc_filter <- pc_filter + sum(pcy < pc_row_low) }
    if (!is.na(pc_row_high)) { 
      abline(h=pc_row_high,col='darkred',lwd=1.5,xpd=FALSE)
      points(pcx[pcy > pc_row_high],pcy[pcy > pc_row_high],col='chartreuse3',pch=1,cex=1)
      pc_filter <- pc_filter + sum(pcy > pc_row_high) }
    ## Add to legend
    ycount <- c(ltable2$count,pc_filter)
    ylabel <- c(as.character(ltable2$label),paste0(cf$PC[pc_row],' filtered = ',pc_filter))
    ycolors <- c(as.character(ltable2$colors),'chartreuse3')
    ltable2 <- cbind.data.frame(ylabel,ycount,ycolors)
    names(ltable2) <- c('label','count','colors')
  }} ## END if

  ## print legend
  legend(x=xleg,y=yleg,legend=ltable2$label,fill=as.character(ltable2$colors),bty='n',bg='white')

} ## END of i LooP

dev.off()

## gzip file
system(paste0('gzip ',opt$data_id,'.caseControl_pca.2d.pdf'))

if (opt$verbose) { cat('\n',paste0('pdf of PCs: ',opt$data_id,'.caseControl_pca.2d.pdf.gz'),'\n') }


} ## END if





## ---- GRAPH of pca (WITH reference data)


if (opt$graph==TRUE & opt$refdata==TRUE) {

cat('\nprinting graph of PCs\n')

## read in unrelated .fam file
fam <- read.table(paste0(opt$pcaer_dir,'/',opt$data_id,'.menv.fam'),h=F,stringsAsFactors=F)

if (nrow(fam) == 0 & opt$verbose) { 
  cat('unable to read in fam file, please check input \n')
  q(status=1) 
} ## END if

## read in PC file
pc <- read.table(paste0(opt$pcaer_dir,'/',opt$data_id,'.menv.mds'),h=T,stringsAsFactors=F)

if (nrow(pc) == 0 & opt$verbose) { 
  cat('unable to read in pc file, please check input \n')
  q(status=1) 
} ## END if

if (cutoffs == TRUE) { 
  ## read in file:
  cf <- read.table(opt$cutoff_file,h=T,stringsAsFactors=F)
}


## get parameters of PC file
pc_names <- names(pc)[4:ncol(pc)]

## find out which PCs to create lines
if (cutoffs == TRUE) { cf_lines <- pc_names %in% cf$PC }

## merge with menv fam file
pca <- merge(pc,fam,by.x='IID',by.y='V2')

pc_xaxis <- pc_names[seq(1,length(pc_names),2)]
pc_yaxis <- pc_names[seq(2,length(pc_names),2)]

aff <- pca$V6
chk <- strsplit(pca$FID,'_'); id <- sapply(chk,'[',2)
aff[id=='pop'] <- 'pop'
aff[id!='pop' & aff==0] <- 'mis'
aff[id!='pop' & aff==-9] <- 'mis'
aff[aff==2] <- 'cas'
aff[aff==1] <- 'con'

legend_id <- substr(pca$FID,start=1,stop=17)
legend_id[legend_id=='mis_pop_afri_afr_'] <- '1kg_pop_afri'
legend_id[legend_id=='mis_pop_amer_amr_'] <- '1kg_pop_amer'
legend_id[legend_id=='mis_pop_euro_eur_'] <- '1kg_pop_euro'
legend_id[legend_id=='mis_pop_asia_asn_'] <- '1kg_pop_asia'
## update fam to cas/con
grp_id <- substr(legend_id,start=1,stop=3)
grp_id[grp_id=='fam' & aff=='cas'] <- 'cas'
grp_id[grp_id=='fam' & aff=='con'] <- 'con'
end_id <- substr(legend_id,start=4,stop=nchar(legend_id))
legend_id <- paste0(grp_id,end_id)

## create label / color table
ltable <- table(legend_id)
if (length(ltable)==5) { lcols <- c('purple','grey','pink','orange','red') }
if (length(ltable)==6) { lcols <- c('purple','grey','pink','orange','red','blue') }
if (length(ltable) > 6) { add <- length(ltable)-6; lcols <- c('purple','grey','pink','orange','red','blue',rep('black',add)) }
ltable <- cbind.data.frame(ltable,lcols)
names(ltable) <- c('label','count','colors')

cols <- rep('red',nrow(pc))
cols[aff=='cas'] <- 'red'
cols[aff=='con'] <- 'blue'
cols[legend_id=='1kg_pop_afri'] <- 'purple'
cols[legend_id=='1kg_pop_amer'] <- 'grey'
cols[legend_id=='1kg_pop_euro'] <- 'orange'
cols[legend_id=='1kg_pop_asia'] <- 'pink'


## -------- PLOT
pdf(paste0(opt$data_id,'.caseControl_pca.2d.pdf'),height=8,width=10)

par(xpd=NA,mar=c(5, 4, 4, 12) + 0.1)

## Loop through PCs
for (i in 1:length(pc_xaxis)) {

  eval(parse(text=paste0('pcx <- pca$',pc_xaxis[i])))
  eval(parse(text=paste0('pcy <- pca$',pc_yaxis[i])))

  ## set up legend coordinates
  xleg <- range(pcx)[2] + 0.05*(range(pcx)[2] - range(pcx)[1])
  yleg <- range(pcy)[2] - 0.05*(range(pcy)[2] - range(pcy)[1])

  plot(pcx,pcy,xlab=pc_xaxis[i],ylab=pc_yaxis[i],col=cols,cex=0.5,cex.axis=0.7)

  ## grid lines
  abline(v=seq(round(range(pcx)[1],2),round(range(pcx)[2],2),0.01),col='lightgrey',xpd=F)
  abline(h=seq(round(range(pcy)[1],2),round(range(pcy)[2],2),0.01),col='lightgrey',xpd=F)

  ## add case/control again
  points(pcx[aff=='mis'],pcy[aff=='mis'],col=cols[aff=='mis'],cex=0.5)
  points(pcx[aff=='con'],pcy[aff=='con'],col=cols[aff=='con'],cex=0.5)
  points(pcx[aff=='cas'],pcy[aff=='cas'],col=cols[aff=='cas'],cex=0.5)

  ## copy original legend table
  ltable2 <- ltable

  ## Get x-axis cutoffs
  if (cutoffs == TRUE) {
    if ((pc_xaxis[i] %in% cf$PC) == TRUE) {
    pc_row <- which(cf$PC %in% pc_xaxis[i])
    pc_row_low <- cf$low[pc_row]; pc_row_high <- cf$high[pc_row]
    pc_filter <- 0
    ## DRAW LINE and highlight points
    if (!is.na(pc_row_low)) {
      abline(v=pc_row_low,col='darkred',lwd=1.5,xpd=FALSE)
      points(pcx[pcx < pc_row_low & aff!='pop'],pcy[pcx < pc_row_low & aff!='pop'],col='palegreen',pch=1,cex=1) 
      pc_filter <- pc_filter + sum(pcx[aff!='pop'] < pc_row_low) }
    if (!is.na(pc_row_high)) { 
      abline(v=pc_row_high,col='darkred',lwd=1.5,xpd=FALSE)
      points(pcx[pcx > pc_row_high & aff!='pop'],pcy[pcx > pc_row_high & aff!='pop'],col='palegreen',pch=1,cex=1) 
      pc_filter <- pc_filter + sum(pcx[aff!='pop'] > pc_row_high) }
    ## Add to legend
    xcount <- c(ltable2$count,pc_filter)
    xlabel <- c(as.character(ltable2$label),paste0(cf$PC[pc_row],' filtered = ',pc_filter))
    xcolors <- c(as.character(ltable2$colors),'palegreen')
    ltable2 <- cbind.data.frame(xlabel,xcount,xcolors)
    names(ltable2) <- c('label','count','colors')
  }} ## END if

  ## Get y-axis cutoffs
  if (cutoffs == TRUE) { 
    if ((pc_yaxis[i] %in% cf$PC) == TRUE) {
    pc_row <- which(cf$PC %in% pc_yaxis[i])
    pc_row_low <- cf$low[pc_row]; pc_row_high <- cf$high[pc_row]
    pc_filter <- 0
    ## DRAW LINE and highlight points
    if (!is.na(pc_row_low)) { 
      abline(h=pc_row_low,col='darkred',lwd=1.5,xpd=FALSE)
      points(pcx[pcy < pc_row_low & aff!='pop'],pcy[pcy < pc_row_low & aff!='pop'],col='chartreuse3',pch=1,cex=1)
      pc_filter <- pc_filter + sum(pcy[aff!='pop'] < pc_row_low) }
    if (!is.na(pc_row_high)) { 
      abline(h=pc_row_high,col='darkred',lwd=1.5,xpd=FALSE)
      points(pcx[pcy > pc_row_high & aff!='pop'],pcy[pcy > pc_row_high & aff!='pop'],col='chartreuse3',pch=1,cex=1)
      pc_filter <- pc_filter + sum(pcy[aff!='pop'] > pc_row_high) }
    ## Add to legend
    ycount <- c(ltable2$count,pc_filter)
    ylabel <- c(as.character(ltable2$label),paste0(cf$PC[pc_row],' filtered = ',pc_filter))
    ycolors <- c(as.character(ltable2$colors),'chartreuse3')
    ltable2 <- cbind.data.frame(ylabel,ycount,ycolors)
    names(ltable2) <- c('label','count','colors')
  }} ## END if

  ## print legend
  legend(x=xleg,y=yleg,legend=ltable2$label,fill=as.character(ltable2$colors),bty='n',bg='white')

} ## END of i LooP

dev.off()

## gzip file
system(paste0('gzip ',opt$data_id,'.caseControl_pca.2d.pdf'))

if (opt$verbose) { cat('\n',paste0('pdf of PCs: ',opt$data_id,'.caseControl_pca.2d.pdf.gz'),'\n') }


} ## END if

