# File: 10_deAnalysis.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: Modelling and selecting DE genes
# Date: 21/9/2020

## load the data
source('header.R')

## load the data
library(RMySQL)

db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)

# get the query
g_did
q = paste0('select MetaFile.* from MetaFile
           where (MetaFile.idData = 49) AND (MetaFile.comment like "%count%")')
dfSample = dbGetQuery(db, q)
dfSample
n = paste0(dfSample$location, dfSample$name)
load(n[2])

## load the metadata i.e. covariates
q = paste0('select Sample.* from Sample where Sample.idData = 49')
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
mData = mData[-i,]
# drop the rows where average across rows is less than 3
i = rowMeans(mData)
table( i < 3)
# FALSE  TRUE 
# 15391 11804 
mData = mData[!(i< 3),]
dim(mData)
# [1] 15391    60

ivProb = apply(mData, 1, function(inData) {
  inData[is.na(inData) | !is.finite(inData)] = 0
  inData = as.logical(inData)
  lData = list('success'=sum(inData), fail=sum(!inData))
  return(mean(rbeta(1000, lData$success + 0.5, lData$fail + 0.5)))
})

hist(ivProb)

library(DESeq2)
sf = estimateSizeFactorsForMatrix(mData)
mData.norm = sweep(mData, 2, sf, '/')

identical(colnames(mData.norm), as.character(dfSample$id))

## delete sample section after testing
mData.norm = round(mData.norm, 0)

# split the data into parts to run the model faster
# as too many transcripts
# iIndex = 1:nrow(mData.norm)
# iSplit = cut(iIndex, breaks = 5, include.lowest = T, labels = 1:5)
# set.seed(123)
# i = sample(1:nrow(mData.norm), 100, replace = F)
# dfData = data.frame(t(mData.norm[i,]))
#dfData = data.frame(t(mData.norm[iSplit == 6,]))

dfData = data.frame(t(mData.norm))
dim(dfData)
dfData = stack(dfData)

## create covariates for modelling
str(dfSample)
fPid = factor(gsub('T0|T1', '', dfSample$title))
nlevels(fPid)
fTime = factor(dfSample$group1)
nlevels(fTime)
dfData$fTreatment = fTime
dfData$fPid = fPid
dfData = droplevels.data.frame(dfData)

dfData$Coef.1 = factor(dfData$fTreatment:dfData$ind)
dfData$Coef.2 = factor(dfData$fPid:dfData$ind)
str(dfData)

# # setup the model
# library(lme4)
# fit.lme1 = glmer.nb(values ~ 1 + (1 | Coef.1) + (1 | Coef.2), data=dfData)
# summary(fit.lme1)
# fit.lme2 = glmer.nb(values ~ 1 + (1 | Coef.1), data=dfData)
# summary(fit.lme2)
# anova(fit.lme1, fit.lme2)
# ran = ranef(fit.lme1, condVar=F)
# 
# plot(log(fitted(fit.lme1)), resid(fit.lme1), pch=20, cex=0.7)
# lines(lowess(log(fitted(fit.lme1)), resid(fit.lme1)), col=2)

## setup the stan model
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

stanDso = rstan::stan_model(file='nbinomResp2RandomEffectsMultipleScales.stan')

## calculate hyperparameters for variance of coefficients
# l = gammaShRaFromModeSD(sd(log(dfData$values+0.5)), 2*sd(log(dfData$values+0.5)))
# # ## set initial values
# ran = ranef(fit.lme1)
# r1 = ran$Coef
# r2 = ran$Coef.adj1
# r3 = ran$Coef.adj2
# 
# initf = function(chain_id = 1) {
#   list(sigmaRan1 = 1, sigmaRan2=1)
# }

## subset the data to get the second level of nested parameters
## this is done to avoid loops in the stan script to map the scale parameters
## of each ind/gene to the respective set of coefficients for jitters
d = dfData[!duplicated(dfData$Coef.1), ]
d2 = dfData[!duplicated(dfData$Coef.2), ]

lStanData = list(Ntotal=nrow(dfData), 
                 Nclusters1=nlevels(dfData$Coef.1),
                 Nclusters2=nlevels(dfData$Coef.2),
                 NScaleBatches1 = nlevels(dfData$ind), # to add a separate scale term for each gene
                 NScaleBatches2 = nlevels(dfData$ind), # to add a separate scale term for each gene
                 NgroupMap1=as.numeric(dfData$Coef.1),
                 NgroupMap2=as.numeric(dfData$Coef.2),
                 NBatchMap1=as.numeric(d$ind), # this is where we use the second level mapping
                 NBatchMap2=as.numeric(d2$ind), # this is where we use the second level mapping
                 Nphi=nlevels(dfData$ind),
                 NphiMap=as.numeric(dfData$ind),
                 y=dfData$values, 
                 #gammaShape=l$shape, gammaRate=l$rate,
                 intercept = mean(log(dfData$values+0.5)), intercept_sd= sd(log(dfData$values+0.5))*3)

initf = function(chain_id = 1) {
  list(sigmaRan1 = rep(1, times=lStanData$NScaleBatches1),
       sigmaRan2= rep(0.1, times=lStanData$NScaleBatches2),
       rGroupsJitter1_scaled = rep(0, times=lStanData$Nclusters1),
       rGroupsJitter2_scaled = rep(0, times=lStanData$Nclusters2),
       phi_scaled=rep(15, times=lStanData$Nphi))
}


ptm = proc.time()

fit.stan = sampling(stanDso, data=lStanData, iter=800, chains=4,
                    pars=c('sigmaRan1',
                           #'sigmaRan2',
                           'phi',
                           #'mu',
                           'rGroupsJitter1',
                           #'rGroupsJitter2',
                           'betas'
                           #'phi_scaled'
                    ),
                    cores=4, #control=list(adapt_delta=0.99, max_treedepth = 11),
                    init=initf)
#save(fit.stan, file='results/fit.stan.salmon.nb_21Sep.rds')
ptm.end = proc.time()
print(fit.stan, c('sigmaRan1'), digits=3)
print(fit.stan, c('phi'), digits=3)
print(fit.stan, c('rGroupsJitter1'))
traceplot(fit.stan, c('sigmaRan1[1]'))
traceplot(fit.stan, c('sigmaRan1[2]'))
traceplot(fit.stan, c('rGroupsJitter1[1]', 'sigmaRan1[1]'))

## get the coefficient of interest - Modules in our case from the random coefficients section
mCoef = extract(fit.stan)$rGroupsJitter1
dim(mCoef)
# # ## get the intercept at population level
# iIntercept = as.numeric(extract(fit.stan)$betas)
# ## add the intercept to each random effect variable, to get the full coefficient
# mCoef = sweep(mCoef, 1, iIntercept, '+')

## function to calculate statistics for differences between coefficients
getDifference = function(ivData, ivBaseline){
  stopifnot(length(ivData) == length(ivBaseline))
  # get the difference vector
  d = ivData - ivBaseline
  # get the z value
  z = mean(d)/sd(d)
  # get 2 sided p-value
  p = pnorm(-abs(mean(d)/sd(d)))*2
  return(list(z=z, p=p))
}

## split the data into the comparisons required
d = data.frame(cols=1:ncol(mCoef), mods=levels(dfData$Coef.1))
# the split is done below on : symbol, but factor name has a : symbol due
# to creation of interaction earlier, do some acrobatics to sort that issue
## split this factor into sub factors
f = strsplit(as.character(d$mods), ':')
d = cbind(d, do.call(rbind, f))
head(d)
#d$`1` = d$`1`:d$`2`
#d = d[,-4]
colnames(d) = c(colnames(d)[1:2], c('fBatch', 'ind'))
str(d)
d$split = factor(d$ind)

levels(d$fBatch)
## repeat this for each comparison

## get a p-value for each comparison
l = tapply(d$cols, d$split, FUN = function(x, base='T0', deflection='T1') {
  c = x
  names(c) = as.character(d$fBatch[c])
  dif = getDifference(ivData = mCoef[,c[deflection]], ivBaseline = mCoef[,c[base]])
  r = data.frame(ind= as.character(d$ind[c[base]]), coef.base=mean(mCoef[,c[base]]), 
                 coef.deflection=mean(mCoef[,c[deflection]]), zscore=dif$z, pvalue=dif$p)
  r$difference = r$coef.deflection - r$coef.base
  #return(format(r, digi=3))
  return(r)
})

dfResults = do.call(rbind, l)
dfResults$adj.P.Val = p.adjust(dfResults$pvalue, method='BH')

### plot the results
dfResults$logFC = dfResults$difference
dfResults$P.Value = dfResults$pvalue
library(org.Mm.eg.db)
## remove X from annotation names
dfResults$ind = gsub('X', '', as.character(dfResults$ind))
df = AnnotationDbi::select(org.Mm.eg.db, keys = as.character(dfResults$ind), columns = 'SYMBOL', keytype = 'ENTREZID')
i = match(dfResults$ind, df$ENTREZID)
df = df[i,]
dfResults$SYMBOL = df$SYMBOL
identical(dfResults$ind, df$ENTREZID)
## produce the plots 
f_plotVolcano(dfResults, 'ind:WT vs ni:WT')#, fc.lim=c(-2.5, 2.5))
f_plotVolcano(dfResults, 'ind:WT vs ni:WT', fc.lim=range(dfResults$logFC))

m = tapply(dfData$values, dfData$ind, mean)
i = match(rownames(dfResults), names(m))
m = m[i]
identical(names(m), rownames(dfResults))
plotMeanFC(log(m), dfResults, 0.01, 'lei vs nl')
table(dfResults$adj.P.Val < 0.01)
quantile(round(abs(dfResults$logFC),3), 0:10/10)
table(dfResults$adj.P.Val < 0.01 & abs(dfResults$logFC) > 0.5)
## save the results 
write.csv(dfResults, file='results/DEAnalysis_ind:MELWTVSn.i:MELWT.xls')

######### do a comparison with deseq2
str(dfSample.2)
f = as.character(dfSample.2$fReplicates)
f = gsub('(.+)-\\d+$', replacement = '\\1', f)
dfDesign = data.frame(Treatment = factor(dfSample.2$group1):factor(dfSample.2$group2), 
                      Patient=factor(f),
                      row.names=colnames(mData))

oDseq = DESeqDataSetFromMatrix(mData, dfDesign, design = ~ Treatment)# + Patient)
oDseq = DESeq(oDseq)

plotDispEsts(oDseq)
levels(dfDesign$Treatment)
oRes = results(oDseq, contrast = c('Treatment', 'ind:MELWT', 'n.i:MELWT'))
plotMA(oRes)
temp = as.data.frame(oRes)
i = match((dfResults$ind), rownames(temp))
temp = temp[i,]
identical((dfResults$ind), rownames(temp))
plot(dfResults$logFC, log(2^temp$log2FoldChange), pch=20)
table(oRes$padj < 0.01)
write.csv(oRes, file='results/DESeq.xls')

r1 = dfResults
r2 = oRes

r1 = r1[order(abs(r1$logFC), decreasing = T),]
head(r1)
r2$logFC = log(2^r2$log2FoldChange)
r2 = r2[order(abs(r2$logFC), decreasing = T),]

r1.top = r1$ind[1:100]
r2.top = rownames(r2)[1:100]
table(r1.top %in% r2.top)