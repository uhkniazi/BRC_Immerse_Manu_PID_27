# File: header.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: global variables
# Date: 7/9/2020


## variables
g_pid = 27
g_did = 50
gcswd = getwd()
gcRemoteDir = "/run/user/1000/gvfs/sftp:host=login.rosalind.kcl.ac.uk,user=k1625253/users/k1625253/scratch/old-scratch_rosalind-legacy-import_2020-01-28/Data/ProjectsData/BRC_Immerse_Manu_PID_27"

p.old = par()

###### utility functions

f_LoadObject = function(r.obj.file)
{
  # temp environment to load object
  env <- new.env()
  # read object
  nm <- load(r.obj.file, env)[1]
  return(env[[nm]])
}


# utility function to calculate gamma distribution used as a prior for scale
gammaShRaFromModeSD = function( mode , sd ) {
  # function changed a little to return jeffery non-informative prior
  if ( mode <=0 || sd <=0 ) return( list( shape=0.5 , rate=0.0001 ) ) 
  rate = ( mode + sqrt( mode^2 + 4 * sd^2 ) ) / ( 2 * sd^2 )
  shape = 1 + mode * rate
  return( list( shape=shape , rate=rate ) )
}

f_plotVolcano = function(dfGenes, main, p.adj.cut = 0.1, fc.lim = c(-3, 3)){
  p.val = -1 * log10(dfGenes$P.Value)
  fc = dfGenes$logFC
  # cutoff for p.value y.axis
  y.cut = -1 * log10(0.01)
  col = rep('lightgrey', times=length(p.val))
  c = which(dfGenes$adj.P.Val < p.adj.cut)
  col[c] = 'red'
  plot(fc, p.val, pch=20, xlab='Fold Change', ylab='-log10 P.Value', col=col, main=main, xlim=fc.lim)
  abline(v = 0, col='grey', lty=2)
  abline(h = y.cut, col='red', lty=2)
  # second cutoff for adjusted p-values
  y.cut = quantile(p.val[c], probs=0.95)
  abline(h = y.cut, col='red')
  # identify these genes
  g = which(p.val > y.cut)
  lab = dfGenes[g, 'SYMBOL']
  text(dfGenes$logFC[g], y = p.val[g], labels = lab, pos=2, cex=0.6)
}

plotMeanFC = function(m, dat, p.cut, title){
  col = rep('grey', length.out=(nrow(dat)))
  col[dat$adj.P.Val < p.cut] = 'red'
  #mh = cut(m, breaks=quantile(m, 0:50/50), include.lowest = T)
  plot(m, dat$logFC, col=col, pch=20, main=title, xlab='log Mean', ylab='logFC', ylim=c(-2, 2), cex=0.6)
}

plot.ercc.proportions = function(ivTotal, ivErcc, fGroupings, ...){
  pErcc = ivErcc/ivTotal * 100
  c = rainbow(nlevels(fGroupings))
  barplot(pErcc[order(fGroupings)], col=c[as.numeric(fGroupings)[order(fGroupings)]],
          ...)
  abline(h = mean(pErcc), lty=2)
}

plot.ercc.MD = function(mLogCounts, fGroupings, ...){
  if (nlevels(fGroupings) != 2) stop('grouping factor should be 2 levels only')
  mResults = apply(mLogCounts, 1, function(x){
    tapply(x, fGroupings, mean)
  })
  mResults = t(mResults)
  ## difference
  d = mResults[,2] - mResults[,1]
  ## choose colour
  c = rep('grey', length.out=length(d))
  c[grep('ercc', names(d), ignore.case = T)] = 'blue'
  plot(rowMeans(mLogCounts), d, pch=20, col=c, xlab='mean log(count+1)', 
       ylab='log difference', ...)
  abline(h = 0, lty=2)
  i = grep('ercc', rownames(mLogCounts), ignore.case = T)
  lines(lowess(rowMeans(mLogCounts)[-i], d[-i]), col='black')
  lines(lowess(rowMeans(mLogCounts)[i], d[i]), col='blue')
}

plot.erccAbsent.MD = function(mLogCounts, fGroupings, ...){
  if (nlevels(fGroupings) != 2) stop('grouping factor should be 2 levels only')
  mResults = apply(mLogCounts, 1, function(x){
    tapply(x, fGroupings, mean)
  })
  mResults = t(mResults)
  ## difference
  d = mResults[,2] - mResults[,1]
  ## choose colour
  c = rep('grey', length.out=length(d))
  plot(rowMeans(mLogCounts), d, pch=20, col=c, xlab='mean log(count+1)', 
       ylab='log difference', ...)
  abline(h = 0, lty=2)
  lines(lowess(rowMeans(mLogCounts), d), col='black')
}