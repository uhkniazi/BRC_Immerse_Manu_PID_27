# File: 11_vennDiagram.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: using the results for each contrast, create venn diagrams and some plots
# Date: 10/2/2021

source('header.R')

lFiles = list.files('dataID50/results/', pattern='DEAnalysis*', full.names = T, ignore.case = T)

ldfData = lapply(lFiles, function(x) as.data.frame(read.csv(x, header=T, row.names=1, stringsAsFactors = F)))
names(ldfData) = lFiles
sapply(ldfData, nrow)

# put everything in one order by row names
rn = rownames(ldfData[[1]])
head(rn)

ldfData = lapply(ldfData, function(df){
  df = df[rn,]
})

sapply(ldfData, function(df) identical(rownames(df), rn))

cvTitle = gsub('dataID50/results//DEAnalysis_', '', names(ldfData))
cvTitle = gsub('.xls', '', cvTitle)

## redefine plot volcano
f_plotVolcano = function(dfGenes, main, p.adj.cut = 0.1, fc.lim = c(-3, 3)){
  p.val = -1 * log10(dfGenes$P.Value)
  fc = dfGenes$logFC
  # cutoff for p.value y.axis
  y.cut = -1 * log10(0.01)
  col = rep('lightgrey', times=length(p.val))
  c = which(dfGenes$adj.P.Val < p.adj.cut)
  col[c] = 'pink4'
  col[which(dfGenes$adj.P.Val < p.adj.cut & (dfGenes$logFC) < 0)] = 'skyblue'
  col[which(dfGenes$adj.P.Val < p.adj.cut & (dfGenes$logFC) > 1.5)] = 'red'
  col[which(dfGenes$adj.P.Val < p.adj.cut & (dfGenes$logFC) < -1.5)] = 'royalblue4'
  plot(fc, p.val, pch=20, xlab='log Fold Change', ylab='-log10 P.Value', col=col, main=main, xlim=fc.lim, cex=0.6)
  abline(v = 0, col='grey', lty=2)
  #abline(h = y.cut, col='red', lty=2)
  # second cutoff for adjusted p-values
  y.cut = quantile(p.val[c], probs=0.99)
  #abline(h = y.cut, col='red')
  # identify these genes
  g = which(p.val > y.cut)
  g = c(g, which(dfGenes$logFC < -1.5 & dfGenes$adj.P.Val < p.adj.cut))
  g = c(g, which(dfGenes$logFC > 2 & dfGenes$adj.P.Val < p.adj.cut))
  lab = dfGenes[g, 'SYMBOL']
  text(dfGenes$logFC[g], y = p.val[g], labels = lab, pos=2, cex=0.8)
}


for (i in seq_along(cvTitle)){
  df = ldfData[[1]]
  hist(df$logFC, xlab='Log Fold Change', ylab='', main=paste('SepT0 vs CarT0',sep=''))
  f_plotVolcano(df, main = 'SepT0 vs CarT0', 0.01, fc.lim=c(-3, 3))
}

legend('topleft', legend = c('FDR=1% & logFC > 0', 
                             'FDR=1% & logFC > 1.5',
                             'FDR=1% & logFC < 0',
                             'FDR=1% & logFC < -1.5'), col = c('pink4', 'red', 'skyblue', 'royalblue4'), pch=20, bty='n', cex=1.2)
names(ldfData)
## select significant genes
round(quantile(abs(ldfData[[1]]$logFC), 0:10/10), 3)
table(ldfData[[1]]$adj.P.Val < 0.01)
dfContrast1.sub = ldfData[[1]][ldfData[[1]]$adj.P.Val < 0.01,]

