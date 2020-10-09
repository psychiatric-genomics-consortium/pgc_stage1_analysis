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
  'famex'   , 'f', 1, "character",
  'refdata'         , 'r', 1, "logical",
  'pc_count'     , 'p', 1, "integer",
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
if ( is.null(opt$verbose ) ) { opt$verbose = TRUE }
if ( is.null(opt$data_id    ) ) { opt$data_id    = 'dataset'     }
if ( is.null(opt$pcaer_dir      ) ) { opt$pcaer_dir      = '.'     }
if ( is.null(opt$pc_count   ) ) { opt$pc_count   = 20    }
if ( is.null(opt$graph   ) ) { opt$graph   = FALSE    }
if ( is.null(opt$refdata ) ) { opt$refdata = FALSE }

#print verbose output unless FALSE
if ( opt$verbose ) { 
    write("verbose output...",stdout()) 
    print(opt)
    cat(paste0(opt$pcaer_dir,"/",opt$data_id),'\n')
}


#read in .fam file
fam <- read.table(paste0(opt$pcaer_dir,'/',opt$data_id,'.mepr.fam'),stringsAsFactors=F)

if (nrow(fam) == 0 & opt$verbose) { 
    cat('unable to read in .mepr.fam file, please check data_id and pcaer_dir input\n')
    q(status=1) 
} ## END if
if (nrow(fam) > 0 & opt$verbose) { 
    cat('read in fam file\n') 
    cat(paste0(nrow(fam),' samples\n'))
    cat(paste0(sum(fam$V6==1),' controls\n'))
    cat(paste0(sum(fam$V6==2),' cases\n\n'))
} ## END if


#read in .famex file
if (is.null(opt$famex)) { 

famex <- read.table(paste0(opt$pcaer_dir,'/',opt$data_id,'.mepr.famex'),stringsAsFactors=F)

if (nrow(famex) == 0 & opt$verbose) { 
    cat('unable to read in default .mepr.famex file, please check data_id and pcaer_dir input\n')
    q(status=1) 
} ## END if
if (nrow(famex) > 0 & opt$verbose) { 
    cat('read in DEFAULT .famex file\n') 
    cat(paste0(nrow(famex),' related samples to project for PCA\n\n'))
} ## END if
} ## END if

if (!is.null(opt$famex)) { 

famex <- read.table(paste0(opt$famex),stringsAsFactors=F)

if (nrow(famex) == 0 & opt$verbose) { 
    cat('unable to read in specified .famex file, please check data_id and pcaer_dir input\n')
    q(status=1) 
} ## END if
if (nrow(famex) > 0 & opt$verbose) { 
    cat('read in SPECIFIED .famex file\n') 
    cat(paste0(nrow(famex),' related samples to project for PCA\n\n'))
} ## END if

} ## END if


fam$related <- 'unrelated_set'
fam$related[fam$V2 %in% famex$V2] <- 'related_set'

fam2 <- fam[,c('V1','V2','related')]

write.table(fam2,'family_based_pca.clst',col=F,row=F,quo=F,sep='\t')


## --- Run PCA on cluster

system(paste0("/home/unix/sripke/plink_src/plink_1.9_newest/plink --bfile ",opt$pcaer_dir,"/",opt$data_id,".mepr --pca ",opt$pc_count," 'header' --within family_based_pca.clst --pca-cluster-names unrelated_set --out ",opt$data_id,".related_pca"))




## ---- GRAPH of pca (no reference data)


if (opt$graph==TRUE & opt$refdata==FALSE) {

## read in data
pc <- read.table(paste0(opt$data_id,'.related_pca.eigenvec'),h=T,stringsAsFactors=F)
pc_names <- names(pc)[3:ncol(pc)]
pc_xaxis <- pc_names[seq(1,length(pc_names),2)]
pc_yaxis <- pc_names[seq(2,length(pc_names),2)]

## merge with mepr fam file
pc <- merge(pc,fam,by.x='IID',by.y='V2')

# tmp <- strsplit(pc$FID,'_')
# aff <- sapply(tmp,'[',1)
aff <- pc$V6
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
pdf(paste0(opt$data_id,'.related_pca.2d.pdf'),height=8,width=10)

par(xpd=NA,mar=c(5, 4, 4, 12) + 0.1)

## Loop through PCs
for (i in 1:length(pc_xaxis)) {

  eval(parse(text=paste0('pcx <- pc$',pc_xaxis[i])))
  eval(parse(text=paste0('pcy <- pc$',pc_yaxis[i])))

  ## set up legend coordinates
  xleg <- range(pcx)[2] + 0.05*(range(pcx)[2] - range(pcx)[1])
  yleg <- range(pcy)[2] - 0.05*(range(pcy)[2] - range(pcy)[1])

  plot(pcx,pcy,xlab=pc_xaxis[i],ylab=pc_yaxis[i],col=cols,cex=0.5,cex.axis=0.7)

  ## grid lines
  abline(v=seq(round(range(pcx)[1],2),round(range(pcx)[2],2),0.01),col='lightgrey',xpd=F)
  abline(h=seq(round(range(pcy)[1],2),round(range(pcy)[2],2),0.01),col='lightgrey',xpd=F)

  ## add case/control again
  points(pcx[aff=='cas'],pcy[aff=='cas'],col=cols[aff=='cas'],cex=0.5)
  points(pcx[aff=='con'],pcy[aff=='con'],col=cols[aff=='con'],cex=0.5)

  legend(x=xleg,y=yleg,legend=ltable$label,fill=as.character(ltable$colors),bty='n',bg='white')

} ## END of i LooP

dev.off()

## gzip file
system(paste0('gzip ',opt$data_id,'.related_pca.2d.pdf'))

if (opt$verbose) { cat('\n',paste0('pdf of PCs: ',opt$data_id,'.related_pca.2d.pdf.gz'),'\n') }


} ## END if





## ---- GRAPH of pca (with 1KG reference data)


if (opt$graph==TRUE & opt$refdata==TRUE) {

## read in data
pc <- read.table(paste0(opt$data_id,'.related_pca.eigenvec'),h=T,stringsAsFactors=F)
pc_names <- names(pc)[3:ncol(pc)]
pc_xaxis <- pc_names[seq(1,length(pc_names),2)]
pc_yaxis <- pc_names[seq(2,length(pc_names),2)]

## merge with mepr fam file
pc <- merge(pc,fam,by.x='IID',by.y='V2')
aff <- pc$V6
aff[aff==2] <- 'cas'
aff[aff==1] <- 'con'


# tmp <- strsplit(pc$FID,'_')
# aff <- sapply(tmp,'[',1)

legend_id <- substr(pc$FID,start=1,stop=17)
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
pdf(paste0(opt$data_id,'.related_pca.2d.pdf'),height=8,width=10)

par(xpd=NA,mar=c(5, 4, 4, 12) + 0.1)

## Loop through PCs
for (i in 1:length(pc_xaxis)) {

  eval(parse(text=paste0('pcx <- pc$',pc_xaxis[i])))
  eval(parse(text=paste0('pcy <- pc$',pc_yaxis[i])))

  ## set up legend coordinates
  xleg <- range(pcx)[2] + 0.05*(range(pcx)[2] - range(pcx)[1])
  yleg <- range(pcy)[2] - 0.05*(range(pcy)[2] - range(pcy)[1])

  plot(pcx,pcy,xlab=pc_xaxis[i],ylab=pc_yaxis[i],col=cols,cex=0.5,cex.axis=0.7)

  ## grid lines
  abline(v=seq(round(range(pcx)[1],2),round(range(pcx)[2],2),0.01),col='lightgrey',xpd=F)
  abline(h=seq(round(range(pcy)[1],2),round(range(pcy)[2],2),0.01),col='lightgrey',xpd=F)

  ## add case/control again
  points(pcx[aff=='cas'],pcy[aff=='cas'],col=cols[aff=='cas'],cex=0.5)
  points(pcx[aff=='con'],pcy[aff=='con'],col=cols[aff=='con'],cex=0.5)

  legend(x=xleg,y=yleg,legend=ltable$label,fill=as.character(ltable$colors),bty='n',bg='white')

} ## END of i LooP

dev.off()

## gzip file
system(paste0('gzip ',opt$data_id,'.related_pca.2d.pdf'))

if (opt$verbose) { cat('\n',paste0('pdf of PCs: ',opt$data_id,'.related_pca.2d.pdf.gz'),'\n') }

} ## END if



#signal success and exit.
if (opt$verbose) { cat('\nSUCCESS!\n\n') }
q(status=0)



## ================= END of family_based_pca.R


