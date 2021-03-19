# File: 11_vennDiagram.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: using the results for each contrast, create venn diagrams and some plots
# Date: 19/3/2021

source('header.R')

cvFile = 'dataID51/results/Tcells/DEAnalysis_SC_T0vsHC_T0.xls'
dfData = read.csv(cvFile, header=T, stringsAsFactors = F, row.names=1)

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
  #g = which(dfGenes$adj.P.Val < p.adj.cut)
  #g = c(g, which(dfGenes$logFC < -1.5 & dfGenes$adj.P.Val < p.adj.cut))
  #g = c(g, which(dfGenes$logFC > 1.5 & dfGenes$adj.P.Val < p.adj.cut))
  lab = dfGenes[g, 'SYMBOL']
  text(dfGenes$logFC[g], y = p.val[g], labels = lab, pos=2, cex=0.8)
}

range(dfData$logFC)
cvFile
f_plotVolcano(dfData, main = 'T Cells: SCT0 vs HCT0', 0.01, fc.lim=c(-4, 8.5))

legend('topleft', legend = c('FDR=10% & logFC > 0', 
                             'FDR=10% & logFC > 1.5',
                             'FDR=10% & logFC < 0',
                             'FDR=10% & logFC < -1.5'), col = c('pink4', 'red', 'skyblue', 'royalblue4'), pch=20, bty='n', cex=0.8)


