## === Rscript: related_samples_filter.R

## USAGE:
# Rscript related_samples_filter.R \
# [ricopili .genome file] \
# [ricopili .famex file] \
# [output filepath]


## EXAMPLE:
# Rscript related_samples_filter.R \
# pcaer_dataID/dataID.mepr.genome \
# pcaer_dataID/dataID.mepr.famex \
# dataID


args <- commandArgs(TRUE)
gn_file <- args[1]
fx_file <- args[2]
output_filepath <- args[3]

if (file.info(fx_file)$size == 0) {
  cat('\n no IDs in .famex file to remove \n')
  q(status=1)
}

gn <- read.table(gn_file,h=F,stringsAsFactors=F)
gn <- subset(gn,gn$V1!='FID1')
names(gn) <- c("FID1","IID1","FID2","IID2","RT","EZ","Z0","Z1","Z2","PI_HAT","PHE","DST","PPC","RATIO")

fx <- read.table(fx_file,stringsAsFactors=F)

## restructured file

# - FID1 / IID1 / AFF1
# - FID2 / IID2 / AFF2
# - PI_HAT / IID_removed / removed_aff

## get original FID
tmp <- strsplit(gn$FID1,'*',fixed=T)
gn$FID3 <- sapply(tmp,'[',2)
tmp <- strsplit(gn$FID2,'*',fixed=T)
gn$FID4 <- sapply(tmp,'[',2)

## get Ricopili FID cas/con
tmp <- strsplit(gn$FID1,'_',fixed=T)
gn$AFF1 <- sapply(tmp,'[',1)
tmp <- strsplit(gn$FID2,'_',fixed=T)
gn$AFF2 <- sapply(tmp,'[',1)

## filter .genome to Ricopili default
gn2 <- subset(gn,gn$PI_HAT > 0.2)


gn2$removed_IID <- NA
gn2$removed_AFF <- NA
gn2$retained_IID <- NA
gn2$retained_AFF <- NA

## index on IID1
gn2$removed_IID[gn2$IID1 %in% fx$V2] <- gn2$IID1[gn2$IID1 %in% fx$V2]
gn2$removed_AFF[gn2$IID1 %in% fx$V2] <- gn2$AFF1[gn2$IID1 %in% fx$V2]
gn2$retained_IID[!gn2$IID1 %in% fx$V2] <- gn2$IID1[!gn2$IID1 %in% fx$V2]
gn2$retained_AFF[!gn2$IID1 %in% fx$V2] <- gn2$AFF1[!gn2$IID1 %in% fx$V2]
## index on IID2
gn2$removed_IID[gn2$IID2 %in% fx$V2] <- gn2$IID2[gn2$IID2 %in% fx$V2]
gn2$removed_AFF[gn2$IID2 %in% fx$V2] <- gn2$AFF2[gn2$IID2 %in% fx$V2]
gn2$retained_IID[!gn2$IID2 %in% fx$V2] <- gn2$IID2[!gn2$IID2 %in% fx$V2]
gn2$retained_AFF[!gn2$IID2 %in% fx$V2] <- gn2$AFF2[!gn2$IID2 %in% fx$V2]


gn3 <- gn2[,c('FID3','IID1','AFF1','FID4','IID2','AFF2','PI_HAT','removed_IID','removed_AFF')]
colnames(gn3)[1] <- 'FID1'
colnames(gn3)[4] <- 'FID2'

## write file:
write.table(gn3,paste0(output_filepath,'.related_filtered.txt'),col=T,row=F,quo=F,sep='\t')


## Removed FID/IID
gn2$removed_FID[gn2$IID1 %in% fx$V2] <- gn2$FID1[gn2$IID1 %in% fx$V2]
gn2$removed_FID[gn2$IID2 %in% fx$V2] <- gn2$FID2[gn2$IID2 %in% fx$V2]
gn4 <- gn2[,c('removed_FID','removed_IID')]

## write file:
write.table(gn4,paste0(output_filepath,'.related_filtered.ids'),col=F,row=F,quo=F,sep='\t')


## Print out results
cat('Related_samples_filter output:\nSend errors to: howrigan@broadinstitute.org\n\n')

cat(paste(nrow(gn),'pairs in .genome file (PI_HAT > .1)\n'))
cat(paste(nrow(fx),'samples removed in .famex file (PI_HAT > .2)\n'))

ids <- c(gn$IID1,gn$IID2)
cat(paste('table of IIDs in .genome file:\n'))
print(table(table(ids)))

dups <- as.vector(ids[duplicated(ids)])
dup_removed <- dups %in% fx$V2
xx <- cbind.data.frame(dups,dup_removed); names(xx) <- c('IID','removed')
cat('\nduplicated IDs:\n')
xx

## cases/controls removed 
cat('\nREMOVED case/control count:\n')
table(gn2$removed_AFF[!duplicated(gn2$removed_IID)])

## cases/controls retained 
cat('\nRETAINED case/control count:\n')
table(gn2$retained_AFF[!duplicated(gn2$retained_IID)])

## file output
cat(paste0('\n\n\noutput_file:',output_filepath,'.related_filtered.txt\n'))
system(paste0('column -t ',output_filepath,'.related_filtered.txt'))


## === END of Rscript: related_samples_filter.R

