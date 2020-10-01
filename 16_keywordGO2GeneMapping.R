# File: 16_keywordGO2GeneMapping.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: using specific keywords to search database ids, reverse map the ids to genes
# Date: 1/10/2020

source('header.R')

#### To translate the data into a more meaningful biological context and to 
# characterize more thoroughly sets of functionally related genes
# organize the differentially expressed datasets into gene ontology groupings (figures)?
library(GOstats)
library(org.Hs.eg.db)

goTest = function(cvSeed, univ = keys(org.Hs.eg.db, 'ENTREZID')){
  ## set up universe background
  dfUniv = AnnotationDbi::select(org.Hs.eg.db, keys = univ, columns = c('GO'), keytype = 'ENTREZID')
  dfUniv = na.omit(dfUniv)
  univ = unique(dfUniv$ENTREZID)
  
  ## make hypergeometric test object for each type, CC, BP and MF
  params = new('GOHyperGParams', geneIds=unique(cvSeed),
               annotation='org.Hs.eg.db',
               universeGeneIds=univ,
               ontology='BP',
               pvalueCutoff= 0.01,
               conditional=FALSE,
               testDirection='over')
  
  oGOStat = tryCatch(hyperGTest(params), error=function(e) NULL)
  return(oGOStat)
}

## load the common binary matrix of DE genes created in earlier results
dfCommonGenes = read.csv('results/commonDEGenes.xls', header=T, row.names=1)
head(dfCommonGenes)

## as the analysis is done on gene symbols, remove duplicates
f = duplicated(dfCommonGenes$Symbol)
table(f)
dfCommonGenes = dfCommonGenes[!f,]

## replace the transcript ids with enterez ids
dfAnnotation = read.csv(file.choose(), stringsAsFactors = F, sep='\t', header=F)
table(rownames(dfCommonGenes) %in% dfAnnotation$V1)
dfAnnotation = dfAnnotation[dfAnnotation$V1 %in% rownames(dfCommonGenes), ]
dfCommonGenes = dfCommonGenes[rownames(dfCommonGenes) %in% dfAnnotation$V1, ]
table(rownames(dfCommonGenes) %in% dfAnnotation$V1)
## put the annotation table in same order
i = match(rownames(dfCommonGenes), dfAnnotation$V1)
dfAnnotation = dfAnnotation[i,]
# sanity check
identical(rownames(dfCommonGenes), dfAnnotation$V1)
rownames(dfCommonGenes) = as.character(dfAnnotation$V2)
identical(rownames(dfCommonGenes), as.character(dfAnnotation$V2))
head(dfCommonGenes)

## perform analysis on most common groups
iGroups = sort(table(dfCommonGenes$groups), decreasing = T)
iGroups = names(iGroups)

lGO.results = lapply(iGroups, function(group){
  return(goTest(rownames(dfCommonGenes)[dfCommonGenes$groups == group]))
})

names(lGO.results) = iGroups

oFile.go = file('results/GO_groups.csv', 'wt')
temp = sapply(as.character(iGroups), function(group){
  p1 = paste('Contrast Comparison group ', group)
  df = summary(lGO.results[[group]])
  p2 = paste(colnames(df), collapse = ',')
  writeLines(p1, oFile.go)
  writeLines(p2, oFile.go)
  sapply(1:20, function(x){
    p3 = gsub(',', replacement = '-', df[x,])
    p3 = paste(p3, collapse=',')
    writeLines(p3, oFile.go)
  })
})

close(oFile.go)

#######################################################
####### find particular types of genes from pathway keyword search
#dfGO = AnnotationDbi::select(org.Hs.eg.db, keys = 'HEPHL1', columns = c('GO'), keytype = 'SYMBOL')
dfGO = AnnotationDbi::select(org.Hs.eg.db, keys = rownames(dfCommonGenes), columns = c('GO'), keytype = 'ENTREZID')
dfGO = dfGO[dfGO$ONTOLOGY == 'BP', ]
dfGO = na.omit(dfGO)
dim(dfGO)

library(GO.db)
columns(GO.db)
dfGO = AnnotationDbi::select(GO.db, keys=as.character(unique(dfGO$GO)), columns=columns(GO.db), keytype='GOID')
dim(dfGO)
## keyword search
i = grep('cornification', dfGO$DEFINITION, ignore.case = T)
length(i)
temp = dfGO[i,]

# i = grep('blood cell', dfGO$DEFINITION, ignore.case = T)
# length(i)
# temp = rbind(temp, dfGO[i,])
dfGO.keyword = temp

### work back to original gene list
dfGO = AnnotationDbi::select(org.Hs.eg.db, keys = dfGO.keyword$GOID, keytype = c('GO'), columns = 'ENTREZID')
dim(dfGO)
head(dfGO)
table(dfGO$ENTREZID %in% rownames(dfCommonGenes))
cvGenes.keyword = dfGO$ENTREZID[(dfGO$ENTREZID %in% rownames(dfCommonGenes))]
cvGenes.keyword = unique(cvGenes.keyword)
length(cvGenes.keyword)
## go up to stan section to load the d.bk dataframe
d.bk = d[as.character(d$split.2) %in% cvGenes.keyword,]
d.bk = droplevels.data.frame(d.bk)
library(org.Hs.eg.db)
df = AnnotationDbi::select(org.Hs.eg.db, keys = as.character(d.bk$split.2), columns = 'SYMBOL', keytype = 'ENTREZID')
i = match(as.character(d.bk$split.2), df$ENTREZID)
df = df[i,]
d.bk$SYMBOL = df$SYMBOL
identical(as.character(d.bk$split.2), df$ENTREZID)
head(d.bk)
d.bk$coef = colMeans(mCoef[,d.bk$cols])
d.bk$differentiated = factor(d.bk$fBatch, levels=c('Non-lesional', 'Lesional'))

library(lattice)
xyplot(coef ~ differentiated | SYMBOL, data=d.bk, type=c('l', 'p'), scales=list(relation='free', x=list(cex=0.7), y=list(cex=0.7)), 
       ylab='Model Estimated log Average', main=list(label='profile of GO keyword: transcription factor', cex=0.8),
       xlab='Condition')

## cluster the data on trends of expression
m = split(d.bk$coef, f = d.bk$SYMBOL)
m = do.call(rbind, m)
hc = hclust(dist(t(scale(t(m)))))
plot(hc, main='clustered')
c = cutree(hc, k = 2)
table(c)

## expand the cluster variable after matching it with symbol
i = match(as.character(d.bk$SYMBOL), names(c))
d.bk$SYMBOL.cluster = factor(c[i]):factor(d.bk$SYMBOL)

xyplot(coef ~ differentiated | SYMBOL.cluster, data=d.bk, type=c('l', 'p'), scales=list(relation='free', x=list(cex=0.7), y=list(cex=0.7)), 
       ylab='Model Estimated log Average', main=list(label='profile of GO keyword: Cornification', cex=0.8),
       xlab='Condition')

# dim(m)
# 
# ## set names for columns
# cn = split(d.bk$fBatch, f = d.bk$SYMBOL)
# colnames(m) = as.character(cn[[1]])
# head(m)
# df = stack(data.frame(m))
# head(df)
# 
# df = data.frame(df, clustering=factor(c), symbols=rownames(m))
# head(df)
# #df$ind = factor(gsub('X', '', df$ind))
# l = factor(df$clustering:df$symbols)
# l = levels(l)
# l = gsub('\\d:', '', l)
# xyplot(values ~ ind | clustering:symbols, data=df, type=c('l', 'p'), scales=list(relation='free', x=list(cex=0.7), y=list(cex=0.7)), 
#        ylab='Model Estimated log Deflections from Intercept', main=list(label='Tooth Development Genes DE expressed at 3 time points in WT', cex=0.8),
#        strip=strip.custom(factor.levels=l))

#######################################################


################################################################################