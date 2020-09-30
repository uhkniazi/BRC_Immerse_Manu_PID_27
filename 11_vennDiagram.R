# File: 11_vennDiagram.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: using the results for each contrast, create venn diagrams and some plots
# Date: 30/9/2020

source('header.R')

lFiles = list.files('results/', pattern='DEAnalysis*', full.names = T, ignore.case = T)

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

cvTitle = gsub('results//DEAnalysis_', '', names(ldfData))
cvTitle = gsub('.xls', '', cvTitle)

## redefine plot volcano
f_plotVolcano = function(dfGenes, main, p.adj.cut = 0.1, fc.lim = c(-3, 3)){
  p.val = -1 * log10(dfGenes$P.Value)
  fc = dfGenes$logFC
  # cutoff for p.value y.axis
  y.cut = -1 * log10(0.01)
  col = rep('lightgrey', times=length(p.val))
  c = which(dfGenes$adj.P.Val < p.adj.cut)
  col[c] = 'red'
  plot(fc, p.val, pch=20, xlab='Fold Change', ylab='-log10 P.Value', col=col, main=main, xlim=fc.lim, cex=0.5)
  abline(v = 0, col='grey', lty=2)
  abline(h = y.cut, col='red', lty=2)
  # second cutoff for adjusted p-values
  y.cut = quantile(p.val[c], probs=0.997)
  #abline(h = y.cut, col='red')
  # identify these genes
  g = which(p.val > y.cut)
  lab = dfGenes[g, 'SYMBOL']
  text(dfGenes$logFC[g], y = p.val[g], labels = lab, pos=2, cex=0.6)
}


for (i in seq_along(cvTitle)){
  df = ldfData[[i]]
  hist(df$logFC, xlab='Log Fold Change', ylab='', main=paste('Fold Changes ', cvTitle[i], sep=''))
  f_plotVolcano(df, cvTitle[i], 0.01, fc.lim=c(-5, 5))
}

names(ldfData)
## select significant genes
dfContrast1.sub = ldfData[[1]][ldfData[[1]]$adj.P.Val < 0.01,]

# # library(VennDiagram)
# 
# # create a list for overlaps
# lVenn = list(rownames(dfContrast1.sub), rownames(dfContrast2.sub), rownames(dfContrast3.sub),
#              rownames(dfContrast4.sub), rownames(dfContrast5.sub), rownames(dfContrast6.sub),
#              rownames(dfContrast7.sub), rownames(dfContrast8.sub), rownames(dfContrast9.sub),
#              rownames(dfContrast10.sub)
# )
# names(ldfData)
# names(lVenn) = cvTitle
# # calculate overlaps
# #lVenn.overlap = calculate.overlap(lVenn)
# cvTitle[c(1, 3, 5)]
# venn.diagram(lVenn[c(1, 3, 5)], filename = 'results/venn_all_contrasts_DiffMutVsDiffWT.tif', margin=0.1)
# cvTitle[c(2, 4, 6, 7)]
# venn.diagram(lVenn[c(2, 4, 6, 7)], filename = 'results/venn_all_contrasts_DiffVsUndiff.tif', margin=0.1)
# cvTitle[c(8, 9, 10)]
# venn.diagram(lVenn[c(8, 9, 10)], filename = 'results/venn_all_contrasts_UndiffMutVsUndiffWT.tif', margin=0.1)

## repeat the analysis but separate up and down regulated genes
# select significant genes
dfContrast1.up = dfContrast1.sub[dfContrast1.sub$logFC > 0, ]
dfContrast1.down = dfContrast1.sub[dfContrast1.sub$logFC < 0, ]

# create a list for overlaps
lVenn = list(rownames(dfContrast1.up), rownames(dfContrast1.down)
)

#cvTitle = gsub('results//DEAnalysis(\\w+)VsControl.xls', '\\1', names(ldfData))
cvTitle.up = paste(cvTitle, 'up', sep='-')
cvTitle.down = paste(cvTitle, 'down', sep='-')
o = c(1, 2)
## sanity check
matrix(c(cvTitle.up, cvTitle.down)[o], ncol = 2, byrow = T)
cvTitle = c(cvTitle.up, cvTitle.down)[o]
names(lVenn) = cvTitle

## create a binary matrix
cvCommonGenes = unique(do.call(c, lVenn))
mCommonGenes = matrix(NA, nrow=length(cvCommonGenes), ncol=length(lVenn))
for (i in 1:ncol(mCommonGenes)){
  mCommonGenes[,i] = cvCommonGenes %in% lVenn[[i]]
}
rownames(mCommonGenes) = cvCommonGenes
colnames(mCommonGenes) = names(lVenn)

# create groups in the data based on permutation with repetition: 2 ^ ncol 
# https://www.mathsisfun.com/combinatorics/combinations-permutations.html
mCommonGenes.grp = mCommonGenes
set.seed(123)
dm = dist(mCommonGenes.grp, method='binary')
hc = hclust(dm)
plot(hc)
# cut the tree at the bottom to create groups
cp = cutree(hc, h = 0.1)
# sanity checks
table(cp)
length(cp)
length(unique(cp))

mCommonGenes.grp = cbind(mCommonGenes.grp, cp)

### print and observe this table and select the groups you are interested in
temp = mCommonGenes.grp
temp = (temp[!duplicated(cp),])
temp2 = cbind(temp, table(cp))
rownames(temp2) = NULL
print(temp2)

## put these results together
dfCommonGenes = data.frame(mCommonGenes, sig.pvals=rowSums(mCommonGenes), groups=cp, Symbol=ldfData[[1]][cvCommonGenes, 'SYMBOL'])
## gene names have an X before them remove those if present
##rownames(dfCommonGenes) = (gsub(pattern = '^X', replacement = '', x = rownames(dfCommonGenes)))
head(dfCommonGenes)

write.csv(dfCommonGenes, file='results/commonDEGenes.xls')
