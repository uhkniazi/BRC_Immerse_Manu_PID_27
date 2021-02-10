# File: 10_deAnalysis.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: Modelling and selecting DE genes
# Date: 10/2/2021

## load the data
source('header.R')

## load the data
library(RMySQL)

db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)

# get the query
g_did
q = paste0('select MetaFile.* from MetaFile
           where (MetaFile.idData = 50) AND (MetaFile.comment like "%count%")')
dfSample = dbGetQuery(db, q)
dfSample
n = paste0(dfSample$location, dfSample$name)
load(n[1])

## load the metadata i.e. covariates
q = paste0('select Sample.* from Sample where Sample.idData = 50')
dfSample = dbGetQuery(db, q)
dim(dfSample)
dfSample = na.omit(dfSample)
dim(dfSample)
head(dfSample)
# close connection after getting data
dbDisconnect(db)

## make count matrix
names(lCounts)
mCounts = do.call(cbind, lCounts)
colnames(mCounts) = names(lCounts)

# sanity check
identical(dfSample$id, as.integer(colnames(mCounts)))

mData = mCounts
dim(mData)

i = grep('ercc', rownames(mData), ignore.case = T)
mErcc = mData[i, ]
mData = mData[-i,]

## test 2 different filters
i = rowMeans(mErcc)
table(i < 3)

## similar one used in RUVSeq vignette
i2 = apply(mErcc, 1, function(x) length(x[x>5]) >= 2)
table(i2)

table(names(i[i < 3]) %in% names(i2[!i2]))

mErcc = mErcc[!(i < 3), ]
dim(mErcc)

# drop the rows where average across rows is less than 3
i = rowMeans(mData)
table( i < 3)
# FALSE  TRUE 
# 15665 11530
mData = mData[!(i< 3),]
dim(mData)
# [1] 15665    56

ivProb = apply(mData, 1, function(inData) {
  inData[is.na(inData) | !is.finite(inData)] = 0
  inData = as.logical(inData)
  lData = list('success'=sum(inData), fail=sum(!inData))
  return(mean(rbeta(1000, lData$success + 0.5, lData$fail + 0.5)))
})

hist(ivProb)

#### check ercc quality, this was done in EDA script

library(DESeq2)
mData = rbind(mData, mErcc)
i = grep('ercc', rownames(mData), ignore.case = T)
sf = estimateSizeFactorsForMatrix(mData, controlGenes = i)

######### do a comparison with deseq2
str(dfSample)
dfDesign = data.frame(Treatment = factor(dfSample$group1), 
                      #Patient=fPid,
                      row.names=colnames(mData))

i = grep('ercc', rownames(mData), ignore.case = T)
oDseq = DESeqDataSetFromMatrix(mData, dfDesign, design = ~ Treatment)
oDseq = estimateSizeFactors(oDseq, controlGenes=i)
plot(sf, estimateSizeFactors(oDseq, controlGenes=i))
oDseq = DESeq(oDseq)

plotDispEsts(oDseq)
levels(dfDesign$Treatment)
oRes = results(oDseq, contrast = c('Treatment', 'SepT0', 'CT0'))
plotMA(oRes)
dfResults = as.data.frame(oRes)
i = grep('ercc', rownames(dfResults), ignore.case=T)
dfResults = dfResults[-i,]
dfResults$adj.P.Val = p.adjust(dfResults$pvalue, method='BH')

### plot the results
dfResults$logFC = log(2^dfResults$log2FoldChange)
dfResults$P.Value = dfResults$pvalue
library(org.Hs.eg.db)
## remove X from annotation names
dfResults$ind = as.character(rownames(dfResults))
df = AnnotationDbi::select(org.Hs.eg.db, keys = as.character(dfResults$ind), columns = 'SYMBOL', keytype = 'ENTREZID')
i = match(dfResults$ind, df$ENTREZID)
df = df[i,]
dfResults$SYMBOL = df$SYMBOL
identical(dfResults$ind, df$ENTREZID)
identical(dfResults$SYMBOL, df$SYMBOL)
## produce the plots 
f_plotVolcano(dfResults, 'Sepsis vs Cardiac')#, fc.lim=c(-2.5, 2.5))
f_plotVolcano(dfResults, 'Sepsis vs Cardiac', fc.lim=range(dfResults$logFC))

table(dfResults$adj.P.Val < 0.01)
quantile(round(abs(dfResults$logFC),3), 0:10/10)
table(dfResults$adj.P.Val < 0.01 & abs(dfResults$logFC) > 0.3)
## save the results 
write.csv(dfResults, file='dataID50/results/DEAnalysis_SpikeinSepT0vsCT0.xls')
