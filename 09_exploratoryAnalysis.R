# File: 09_exploratoryAnalysis.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: quality checks on the count matrix with covariates
# Date: 17/9/2020

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
load(n[1])

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
mCounts = oTxImport$counts
colnames(mCounts)

# sanity check
identical(dfSample$id, as.integer(colnames(mCounts)))

mData = mCounts
dim(mData)

## some EDA on raw data before merging replicates
library(downloader)
url = 'https://raw.githubusercontent.com/uhkniazi/CDiagnosticPlots/master/CDiagnosticPlots.R'
download(url, 'CDiagnosticPlots.R')

# load the required packages
source('CDiagnosticPlots.R')
# delete the file after source
unlink('CDiagnosticPlots.R')

oDiag.1 = CDiagnosticPlots(log(mData+0.5), 'Raw')
# the batch variable we wish to colour by, 
# this can be any grouping/clustering in the data capture process
str(dfSample)
dfSample$description[1]
#d = strsplit(d, ';')
#fBatch = factor(sapply(d, function(x) x[2]))
fPid = factor(gsub('T0|T1', '', dfSample$title))
nlevels(fPid)
fTime = factor(dfSample$group1)
nlevels(fTime)
# choose appropriate factor
fBatch = factor(fTime)
levels(fBatch)

boxplot.median.summary(oDiag.1, fBatch, legend.pos = 'topright', axis.label.cex = 0.5)
plot.mean.summary(oDiag.1, fBatch, axis.label.cex = 0.5)
plot.sigma.summary(oDiag.1, fBatch, axis.label.cex = 0.5)
plot.missing.summary(oDiag.1, fBatch, axis.label.cex = 0.5, cex.main=1)
plot.PCA(oDiag.1, fBatch, cex.main=1)
plot.dendogram(oDiag.1, fBatch, labels_cex = 0.8, cex.main=0.7)
## change parameters 
l = CDiagnosticPlotsGetParameters(oDiag.1)
l$PCA.jitter = F
l$HC.jitter = F
oDiag.1 = CDiagnosticPlotsSetParameters(oDiag.1, l)
plot.PCA(oDiag.1, fBatch, legend.pos = 'bottom')
plot.dendogram(oDiag.1, fBatch, labels_cex = 0.7)

## normalise the data, with and without ERCCs
load(file.choose())
mErcc = do.call(cbind, lCounts) 
i = grep('ercc', rownames(mErcc), ignore.case = T)
mErcc = mErcc[i,]
colnames(mErcc) = names(lCounts)
identical(colnames(mCounts), colnames(mErcc))
# drop the rows where average across rows is less than 3
i = rowMeans(mData)
table( i < 3)
# FALSE   TRUE 
# 67688 161066 
mData = mData[!(i< 3),]
dim(mData)
# [1] 67688    60

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

mData.ercc = rbind(mData, mErcc)
sf = estimateSizeFactorsForMatrix(mData.ercc, controlGenes = grep('ercc', rownames(mData.ercc), ignore.case = T))
mData.norm.2 = sweep(mData.ercc, 2, sf, '/')

identical(colnames(mData.norm), as.character(dfSample$id))
identical(colnames(mData.norm.2), as.character(dfSample$id))

## compare the normalised data 
oDiag.2 = CDiagnosticPlots(log(mData.norm+0.5), 'Normalised 1')
oDiag.3 = CDiagnosticPlots(log(mData.norm.2+0.5), 'Normalised ercc')

## change parameters 
l = CDiagnosticPlotsGetParameters(oDiag.2)
l$PCA.jitter = F
l$HC.jitter = F

oDiag.2 = CDiagnosticPlotsSetParameters(oDiag.2, l)
oDiag.3 = CDiagnosticPlotsSetParameters(oDiag.3, l)

# the batch variable we wish to colour by, 
fBatch = fTime
levels(fBatch)
## compare the 2 methods using various plots
par(mfrow=c(1,2))
boxplot.median.summary(oDiag.2, fBatch, legend.pos = 'topright', axis.label.cex = 0.5)
boxplot.median.summary(oDiag.3, fBatch, legend.pos = 'topright', axis.label.cex = 0.5)

plot.mean.summary(oDiag.2, fBatch, axis.label.cex = 1)
plot.mean.summary(oDiag.3, fBatch, axis.label.cex = 1)

plot.sigma.summary(oDiag.2, fBatch, axis.label.cex = 0.5)
plot.sigma.summary(oDiag.3, fBatch, axis.label.cex = 0.5)

plot.missing.summary(oDiag.2, fBatch, axis.label.cex = 0.5, cex.main=1)
plot.missing.summary(oDiag.3, fBatch, axis.label.cex = 0.5, cex.main=1)

plot.PCA(oDiag.2, fBatch, cex=1, csLabels = as.character(fPid), legend.pos = 'topleft')
plot.PCA(oDiag.3, fBatch, cex=1, csLabels = as.character(fPid), legend.pos = 'topright')
plot.dendogram(oDiag.2, fBatch, labels_cex = 0.7)
plot.dendogram(oDiag.3, fBatch, labels_cex = 1)

## import the ERCC mix data from sequencing lab
## data not available

######## modelling of PCA components to assign sources of variance to covariates in the design
par(mfrow=c(1,1))
plot(oDiag.2@lData$PCA$sdev)
plot.PCA(oDiag.2, fBatch)
mPC = oDiag.2@lData$PCA$x[,1:3]

## try a linear mixed effect model to account for varince
library(lme4)
dfData = data.frame(mPC)
dfData = stack(dfData)
str(dfData)
dfData$values = as.numeric(scale(dfData$values))
library(lattice)
densityplot(~ values, data=dfData)
densityplot(~ values | ind, data=dfData, scales=list(relation='free'))

str(dfSample)
dfData$fTreatment = fTime
dfData$fPid = fPid

densityplot(~ values | ind, groups=fTreatment, data=dfData, auto.key = list(columns=2), scales=list(relation='free'))
xyplot(values ~ fTreatment | ind, groups=fPid, data=dfData, type=c('p', 'l'), scales=list(relation='free'))
# format data for modelling
dfData$Coef.1 = factor(dfData$fTreatment:dfData$ind)
dfData$Coef.2 = factor(dfData$fPid:dfData$ind)
str(dfData)

fit.lme1 = lmer(values ~ 1  + (1 | Coef.1), data=dfData)
fit.lme2 = lmer(values ~ 1  + (1 | Coef.1) + (1 | Coef.2), data=dfData)

anova(fit.lme1, fit.lme2)

summary(fit.lme2)
plot((fitted(fit.lme2)), resid(fit.lme2), pch=20, cex=0.7)
lines(lowess((fitted(fit.lme2)), resid(fit.lme2)), col=2)
hist(dfData$values, prob=T)
lines(density(fitted(fit.lme2)))

## fit model with stan with various model sizes
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(rethinking)

stanDso = rstan::stan_model(file='tResponsePartialPooling.stan')

######## models of 3 sizes using stan
str(dfData)
m1 = model.matrix(values ~ Coef.1 - 1, data=dfData)
m2 = model.matrix(values ~ Coef.2 - 1, data=dfData)
m = cbind(m1, m2)

lStanData = list(Ntotal=nrow(dfData), Ncol=ncol(m), X=m,
                 NscaleBatches=2, NBatchMap=c(rep(1, times=nlevels(dfData$Coef.1)),
                                              rep(2, times=nlevels(dfData$Coef.2))),
                 y=dfData$values)

fit.stan.2 = sampling(stanDso, data=lStanData, iter=1000, chains=2, pars=c('betas', 'populationMean', 'sigmaPop', 'sigmaRan',
                                                                           'nu', 'mu', 'log_lik'),
                      cores=2, control=list(adapt_delta=0.99, max_treedepth = 12))
print(fit.stan.2, c('populationMean', 'sigmaPop', 'sigmaRan', 'nu', 'betas'), digits=3)

traceplot(fit.stan.2, 'populationMean')
traceplot(fit.stan.2, 'sigmaPop')
traceplot(fit.stan.2, 'sigmaRan')

## just using the one covariate
m = model.matrix(values ~ Coef.1 - 1, data=dfData)

lStanData = list(Ntotal=nrow(dfData), Ncol=ncol(m), X=m,
                 NscaleBatches=1, NBatchMap=c(rep(1, times=nlevels(dfData$Coef.1))),
                 y=dfData$values)

fit.stan.1 = sampling(stanDso, data=lStanData, iter=1000, chains=2, pars=c('betas', 'populationMean', 'sigmaPop', 'sigmaRan',
                                                                           'nu', 'mu', 'log_lik'),
                      cores=2)
print(fit.stan.1, c('populationMean', 'sigmaPop', 'sigmaRan', 'nu', 'betas'), digits=3)


## some model scores and comparisons
compare(fit.stan.2, fit.stan.1)
plot(compare(fit.stan.2, fit.stan.1))

# plot(LOOPk(fit.stan.2) ~ WAIC(fit.stan.2, pointwise = T))

## plot coefficients
# ct = coeftab(fit.stan.4, fit.stan.2, fit.stan.1)
# rn = rownames(ct@coefs)
# i = grep('betas', rn)
# plot(ct, pars=rn[i[1:10]])
# plot(ct, pars=rn[i])
############### new simulated data
###############
### generate some posterior predictive data
## generate random samples from alternative t-distribution parameterization
## see https://grollchristian.wordpress.com/2013/04/30/students-t-location-scale/
rt_ls <- function(n, df, mu, a) rt(n,df)*a + mu
## follow the algorithm in section 14.3 page 363 in Gelman 2013
simulateOne = function(mu, sigma, nu){
  yrep = rt_ls(length(mu), nu, mu,  sigma)
  return(yrep)
}

## sample n values, 1000 times
mDraws.sim = matrix(NA, nrow = nrow(dfData), ncol=300)
l = extract(fit.stan.2)
for (i in 1:300){
  p = sample(1:nrow(l$mu), 1)
  mDraws.sim[,i] = simulateOne(l$mu[p,], 
                               l$sigmaPop[p],
                               l$nu[p])
}

dim(mDraws.sim)
plot(density(dfData$values), main='posterior predictive density plots, model 1')
apply(mDraws.sim, 2, function(x) lines(density(x), lwd=0.5, col='lightgrey'))
lines(density(dfData$values))

## plot residuals
plot(dfData$values - colMeans(l$mu) ~ colMeans(l$mu))
lines(lowess(colMeans(l$mu), dfData$values - colMeans(l$mu)))
apply(l$mu[sample(1:nrow(l$mu), 100),], 1, function(x) {
  lines(lowess(x, dfData$values - x), lwd=0.5, col=2)
})

## plot the original PCA and replicated data
plot(dfData$values[dfData$ind == 'PC1'], dfData$values[dfData$ind == 'PC2'], 
     col=c(1,2)[as.numeric(dfData$fTreatment[dfData$ind == 'PC1'])], main='PCA Components - original and simulated',
     xlab='PC1', ylab='PC2')
points(rowMeans(mDraws.sim)[dfData$ind == 'PC1'], rowMeans(mDraws.sim)[dfData$ind == 'PC2'],
       col=c(1,2)[as.numeric(dfData$fTreatment[dfData$ind == 'PC1'])], pch='1')

plot(dfData$values[dfData$ind == 'PC1'], dfData$values[dfData$ind == 'PC2'], 
     col=c(1,2)[as.numeric(dfData$fTreatment[dfData$ind == 'PC1'])], main='PCA Components - original and model 3',
     xlab='PC1', ylab='PC2', xlim=c(-3, 3), ylim=c(-2, 2))

apply(mDraws.sim, 2, function(x) {
  points(x[dfData$ind == 'PC1'], x[dfData$ind == 'PC2'],
         col=c(1,2)[as.numeric(dfData$fTreatment[dfData$ind == 'PC1'])], pch=20)
})


# ############## differences in coefficients
# ## get the coefficient of interest - Modules in our case from the random coefficients section
# mCoef = extract(fit.stan.2)$betas
# dim(mCoef)
# ## get the intercept at population level
# iIntercept = as.numeric(extract(fit.stan.2)$populationMean)
# ## add the intercept to each coefficient, to get the full coefficient
# mCoef = sweep(mCoef, 1, iIntercept, '+')
# 
# ## split the data into the comparisons required
# d = data.frame(cols=1:ncol(mCoef), mods=c(levels(dfData$Coef.3), levels(dfData$Coef.4)))#, levels(dfData$Coef.4)))
# # the split is done below on : symbol
# ## split this factor into sub factors
# f = strsplit(as.character(d$mods), ':')
# d = cbind(d, do.call(rbind, f))
# head(d)
# 
# d[d$`2` == 'PC1',]
# ## main effects + interactions
# tapply(dfData$values, dfData$Coef.3, mean)
# iLei.New.PC1 = rowSums(mCoef[,c(1, 5)])
# iNl.New.PC1 = rowSums(mCoef[,c(3, 5)])
# iLei.Old.PC1 = rowSums(mCoef[,c(1, 7)])
# iNl.Old.PC1 = rowSums(mCoef[,c(3, 7)])
# 
# # iWT.PC1.av = rowMeans(cbind(iWT.B2.PC1, iWT.B1.PC1))
# # iKO.PC1.av = rowMeans(cbind(iKO.B2.PC1, iKO.B1.PC1))
# 
# ## main effects
# tapply(dfData$values, dfData$Coef.1, mean)
# summary(mCoef[,1])
# summary(mCoef[,3])
# 
# tapply(dfData$values, dfData$Coef.2, mean)
# mean(mCoef[,5])
# mean(mCoef[,7])
# ##########################################

m = cbind(extract(fit.stan.4)$sigmaRan, extract(fit.stan.4)$sigmaPop) 
dim(m)
m = log(m)
colnames(m) = c('Treatment', 'Genotype', 'TrGt', 'TechnicalRep', 'Residual')
pairs(m, pch=20, cex=0.5, col='grey')

df = stack(data.frame(m[,-5]))
histogram(~ values | ind, data=df, xlab='Log SD', scales=list(relation='free'))

## calculate bayesian p-value for this test statistic
getPValue = function(Trep, Tobs){
  left = sum(Trep <= Tobs)/length(Trep)
  right = sum(Trep >= Tobs)/length(Trep)
  return(min(left, right))
}
## define some test quantities to measure the lack of fit
## define a test quantity T(y, theta)
## variance
T1_var = function(Y) return(var(Y))

## min quantity
T1_min = function(Y){
  return(min(Y))
} 

## max quantity
T1_max = function(Y){
  return(max(Y))
} 

## mean quantity
T1_mean = function(Y){
  return(mean(Y))
} 

## mChecks
ivResp = dfData$values
mChecks = matrix(NA, nrow=4, ncol=1)
rownames(mChecks) = c('Variance', 'Max', 'Min', 'Mean')
colnames(mChecks) = c('model 1')

t1 = apply(mDraws.sim, 2, T1_var)
mChecks['Variance', 1] = getPValue(t1, var(ivResp))

## testing for outlier detection i.e. the minimum value show in the histograms earlier
t1 = apply(mDraws.sim, 2, T1_min)
t2 = T1_min(ivResp)
mChecks['Min',1] = getPValue(t1, t2)

## maximum value
t1 = apply(mDraws.sim, 2, T1_max)
t2 = T1_max(ivResp)
mChecks['Max', 1] = getPValue(t1, t2)

## mean value
t1 = apply(mDraws.sim, 2, T1_mean)
t2 = T1_mean(ivResp)
mChecks['Mean', 1] = getPValue(t1, t2)

mChecks