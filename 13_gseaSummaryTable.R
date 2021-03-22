# File: 13_gseaSummaryTable.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: merge the gsea results for all contrasts in one table
# Date: 19/3/2020

lFiles = list.files('dataID51/results/Bcells/', pattern='*pathways_mSigDb_c2_*', full.names = T, ignore.case = F)

# load the files
ldfData = lapply(lFiles, function(x) read.csv(x, row.names=1))
names(ldfData) = lFiles
sapply(ldfData, nrow)
# put everything in one order by row names
rn = rownames(ldfData[[1]])
head(rn)

ldfData = lapply(ldfData, function(df){
  df = df[rn,]
})

# get the upregulated/downregulated in order
ldfData.up = ldfData[grepl(lFiles, pattern = 'upregulated')]
ldfData.down = ldfData[grepl(lFiles, pattern = 'downregulated')]

## set the names for each contrast
sn = gsub('dataID51/results/Bcells//(.+vs.+)_up.+', '\\1', names(ldfData.up))
names(ldfData.up) = sn

sn = gsub('dataID51/results/Bcells//(.+vs.+)_down.+', '\\1', names(ldfData.down))
names(ldfData.down) = sn

## create a table/matrix of p-values
mMerged.up = sapply(ldfData.up, function(x) x$p.val)
rownames(mMerged.up) = rownames(ldfData.up[[1]])
## create similar table for downregulated
mMerged.down = sapply(ldfData.down, function(x) x$p.val)
rownames(mMerged.down) = rownames(ldfData.down[[1]])

# sanity check
identical(rownames(mMerged.up), rownames(mMerged.down))
colnames(mMerged.up) = paste(colnames(mMerged.up), 'up', sep='-')
colnames(mMerged.down) = paste(colnames(mMerged.down), 'down', sep='-')

mMerged.c2 = cbind(mMerged.up, mMerged.down)
# reorder the columns
colnames(mMerged.c2)
o = c(1, 3, 2, 4)
## sanity check
matrix(colnames(mMerged.c2)[o], ncol = 2, byrow = T)
colnames(mMerged.c2)[o]
mMerged.c2 = mMerged.c2[,o]

# remove na sections
dim(mMerged.c2)
mMerged.c2 = na.omit(mMerged.c2)
dim(mMerged.c2)
head(mMerged.c2)

### create a binary matrix based on cutoffs
getBinaryMatrix = function(mat, cutoff=0.01){
  mat2 = apply(mat, 2, function(x) round(x, 3) <= cutoff)
}

mMerged.c2.bin = getBinaryMatrix(mMerged.c2)

## group this matrix into combinations, permutations with repetition
mMerged.c2.bin.grp = mMerged.c2.bin
set.seed(123)
dm = dist(mMerged.c2.bin.grp, method='binary')
hc = hclust(dm)
plot(hc, labels=F)
# cut the tree at the bottom to create groups
cp = cutree(hc, h = 0.2)
# sanity checks
table(cp)
length(cp)
length(unique(cp))

mMerged.c2.bin.grp = cbind(mMerged.c2.bin.grp, cp)

### print and observe this table and select the groups you are interested in
temp = mMerged.c2.bin.grp
temp = (temp[!duplicated(cp),])
temp2 = cbind(temp, table(cp))
rownames(temp2) = NULL
print(temp2)

## map these names to the cp
groups = cp
sig.pvals = rowSums(mMerged.c2.bin)
dfMerged.c2 = data.frame(round(mMerged.c2, 3), sig.pvals, groups, DB='msigdb-c2')
str(dfMerged.c2)
head(dfMerged.c2)
tail(dfMerged.c2)

########
write.csv(dfMerged.c2, file='dataID51/results/Bcells/gsea_msigdb_c2_merged.xls')

###################### repeat for other databases
################### C3
rm(lFiles)
lFiles = list.files('dataID51/results/Bcells/', pattern='*pathways_mSigDb_c3_*', full.names = T, ignore.case = F)

# load the files
ldfData = lapply(lFiles, function(x) read.csv(x, row.names=1))
names(ldfData) = lFiles
sapply(ldfData, nrow)
# put everything in one order by row names
rn = rownames(ldfData[[1]])
head(rn)

ldfData = lapply(ldfData, function(df){
  df = df[rn,]
})

# get the upregulated/downregulated in order
ldfData.up = ldfData[grepl(lFiles, pattern = 'upregulated')]
ldfData.down = ldfData[grepl(lFiles, pattern = 'downregulated')]

## set the names for each contrast
sn = gsub('dataID51/results/Bcells//(.+vs.+)_up.+', '\\1', names(ldfData.up))
names(ldfData.up) = sn

sn = gsub('dataID51/results/Bcells//(.+vs.+)_down.+', '\\1', names(ldfData.down))
names(ldfData.down) = sn

## create a table/matrix of p-values
mMerged.up = sapply(ldfData.up, function(x) x$p.val)
rownames(mMerged.up) = rownames(ldfData.up[[1]])
## create similar table for downregulated
mMerged.down = sapply(ldfData.down, function(x) x$p.val)
rownames(mMerged.down) = rownames(ldfData.down[[1]])

# sanity check
identical(rownames(mMerged.up), rownames(mMerged.down))
colnames(mMerged.up) = paste(colnames(mMerged.up), 'up', sep='-')
colnames(mMerged.down) = paste(colnames(mMerged.down), 'down', sep='-')

mMerged.c3 = cbind(mMerged.up, mMerged.down)
# reorder the columns
colnames(mMerged.c3)
o = c(1, 3, 2, 4)
## sanity check
matrix(colnames(mMerged.c3)[o], ncol = 2, byrow = T)
colnames(mMerged.c3)[o]
mMerged.c3 = mMerged.c3[,o]

# remove na sections
dim(mMerged.c3)
mMerged.c3 = na.omit(mMerged.c3)
dim(mMerged.c3)
head(mMerged.c3)

### create a binary matrix based on cutoffs
getBinaryMatrix = function(mat, cutoff=0.01){
  mat2 = apply(mat, 2, function(x) round(x, 3) <= cutoff)
}

mMerged.c3.bin = getBinaryMatrix(mMerged.c3)

## group this matrix into combinations, permutations with repetition
mMerged.c3.bin.grp = mMerged.c3.bin
set.seed(123)
dm = dist(mMerged.c3.bin.grp, method='binary')
hc = hclust(dm)
plot(hc, labels=F)
# cut the tree at the bottom to create groups
cp = cutree(hc, h = 0.2)
# sanity checks
table(cp)
length(cp)
length(unique(cp))

mMerged.c3.bin.grp = cbind(mMerged.c3.bin.grp, cp)

### print and observe this table and select the groups you are interested in
temp = mMerged.c3.bin.grp
temp = (temp[!duplicated(cp),])
temp2 = cbind(temp, table(cp))
rownames(temp2) = NULL
print(temp2)

## map these names to the cp
groups = cp
sig.pvals = rowSums(mMerged.c3.bin)
dfMerged.c3 = data.frame(round(mMerged.c3, 3), sig.pvals, groups, DB='msigdb-c3')
str(dfMerged.c3)
head(dfMerged.c3)
tail(dfMerged.c3)

########
write.csv(dfMerged.c3, file='dataID51/results/Bcells/gsea_msigdb_c3_merged.xls')

########################
####### C5
rm(lFiles)
lFiles = list.files('dataID51/results/Bcells/', pattern='*pathways_mSigDb_c5_*', full.names = T, ignore.case = F)

# load the files
ldfData = lapply(lFiles, function(x) read.csv(x, row.names=1))
names(ldfData) = lFiles
sapply(ldfData, nrow)
# put everything in one order by row names
rn = rownames(ldfData[[1]])
head(rn)

ldfData = lapply(ldfData, function(df){
  df = df[rn,]
})

# get the upregulated/downregulated in order
ldfData.up = ldfData[grepl(lFiles, pattern = 'upregulated')]
ldfData.down = ldfData[grepl(lFiles, pattern = 'downregulated')]

## set the names for each contrast
sn = gsub('dataID51/results/Bcells//(.+vs.+)_up.+', '\\1', names(ldfData.up))
names(ldfData.up) = sn

sn = gsub('dataID51/results/Bcells//(.+vs.+)_down.+', '\\1', names(ldfData.down))
names(ldfData.down) = sn

## create a table/matrix of p-values
mMerged.up = sapply(ldfData.up, function(x) x$p.val)
rownames(mMerged.up) = rownames(ldfData.up[[1]])
## create similar table for downregulated
mMerged.down = sapply(ldfData.down, function(x) x$p.val)
rownames(mMerged.down) = rownames(ldfData.down[[1]])

# sanity check
identical(rownames(mMerged.up), rownames(mMerged.down))
colnames(mMerged.up) = paste(colnames(mMerged.up), 'up', sep='-')
colnames(mMerged.down) = paste(colnames(mMerged.down), 'down', sep='-')

mMerged.c5 = cbind(mMerged.up, mMerged.down)
# reorder the columns
colnames(mMerged.c5)
o = c(1, 3, 2, 4)
## sanity check
matrix(colnames(mMerged.c5)[o], ncol = 2, byrow = T)
colnames(mMerged.c5)[o]
mMerged.c5 = mMerged.c5[,o]

# remove na sections
dim(mMerged.c5)
mMerged.c5 = na.omit(mMerged.c5)
dim(mMerged.c5)
head(mMerged.c5)

### create a binary matrix based on cutoffs
getBinaryMatrix = function(mat, cutoff=0.01){
  mat2 = apply(mat, 2, function(x) round(x, 3) <= cutoff)
}

mMerged.c5.bin = getBinaryMatrix(mMerged.c5)

## group this matrix into combinations, permutations with repetition
mMerged.c5.bin.grp = mMerged.c5.bin
set.seed(123)
dm = dist(mMerged.c5.bin.grp, method='binary')
hc = hclust(dm)
plot(hc, labels=F)
# cut the tree at the bottom to create groups
cp = cutree(hc, h = 0.2)
# sanity checks
table(cp)
length(cp)
length(unique(cp))

mMerged.c5.bin.grp = cbind(mMerged.c5.bin.grp, cp)

### print and observe this table and select the groups you are interested in
temp = mMerged.c5.bin.grp
temp = (temp[!duplicated(cp),])
temp2 = cbind(temp, table(cp))
rownames(temp2) = NULL
print(temp2)

## map these names to the cp
groups = cp
sig.pvals = rowSums(mMerged.c5.bin)
dfMerged.c5 = data.frame(round(mMerged.c5, 3), sig.pvals, groups, DB='msigdb-c5')
str(dfMerged.c5)
head(dfMerged.c5)
tail(dfMerged.c5)

########
write.csv(dfMerged.c5, file='dataID51/results/Bcells/gsea_msigdb_c5_merged.xls')

########################################################
############ C7
rm(lFiles)
lFiles = list.files('dataID51/results/Bcells/', pattern='*pathways_mSigDb_c7_*', full.names = T, ignore.case = F)

# load the files
ldfData = lapply(lFiles, function(x) read.csv(x, row.names=1))
names(ldfData) = lFiles
sapply(ldfData, nrow)
# put everything in one order by row names
rn = rownames(ldfData[[1]])
head(rn)

ldfData = lapply(ldfData, function(df){
  df = df[rn,]
})

# get the upregulated/downregulated in order
ldfData.up = ldfData[grepl(lFiles, pattern = 'upregulated')]
ldfData.down = ldfData[grepl(lFiles, pattern = 'downregulated')]

## set the names for each contrast
sn = gsub('dataID51/results/Bcells//(.+vs.+)_up.+', '\\1', names(ldfData.up))
names(ldfData.up) = sn

sn = gsub('dataID51/results/Bcells//(.+vs.+)_down.+', '\\1', names(ldfData.down))
names(ldfData.down) = sn

## create a table/matrix of p-values
mMerged.up = sapply(ldfData.up, function(x) x$p.val)
rownames(mMerged.up) = rownames(ldfData.up[[1]])
## create similar table for downregulated
mMerged.down = sapply(ldfData.down, function(x) x$p.val)
rownames(mMerged.down) = rownames(ldfData.down[[1]])

# sanity check
identical(rownames(mMerged.up), rownames(mMerged.down))
colnames(mMerged.up) = paste(colnames(mMerged.up), 'up', sep='-')
colnames(mMerged.down) = paste(colnames(mMerged.down), 'down', sep='-')

mMerged.c7 = cbind(mMerged.up, mMerged.down)
# reorder the columns
colnames(mMerged.c7)
o = c(1, 3, 2, 4)
## sanity check
matrix(colnames(mMerged.c7)[o], ncol = 2, byrow = T)
colnames(mMerged.c7)[o]
mMerged.c7 = mMerged.c7[,o]

# remove na sections
dim(mMerged.c7)
mMerged.c7 = na.omit(mMerged.c7)
dim(mMerged.c7)
head(mMerged.c7)

### create a binary matrix based on cutoffs
getBinaryMatrix = function(mat, cutoff=0.01){
  mat2 = apply(mat, 2, function(x) round(x, 3) <= cutoff)
}

mMerged.c7.bin = getBinaryMatrix(mMerged.c7)

## group this matrix into combinations, permutations with repetition
mMerged.c7.bin.grp = mMerged.c7.bin
set.seed(123)
dm = dist(mMerged.c7.bin.grp, method='binary')
hc = hclust(dm)
plot(hc, labels=F)
# cut the tree at the bottom to create groups
cp = cutree(hc, h = 0.2)
# sanity checks
table(cp)
length(cp)
length(unique(cp))

mMerged.c7.bin.grp = cbind(mMerged.c7.bin.grp, cp)

### print and observe this table and select the groups you are interested in
temp = mMerged.c7.bin.grp
temp = (temp[!duplicated(cp),])
temp2 = cbind(temp, table(cp))
rownames(temp2) = NULL
print(temp2)

## map these names to the cp
groups = cp
sig.pvals = rowSums(mMerged.c7.bin)
dfMerged.c7 = data.frame(round(mMerged.c7, 3), sig.pvals, groups, DB='msigdb-c7')
str(dfMerged.c7)
head(dfMerged.c7)
tail(dfMerged.c7)

########
write.csv(dfMerged.c7, file='dataID51/results/Bcells/gsea_msigdb_c7_merged.xls')

#######################################################
### C8
rm(lFiles)
lFiles = list.files('dataID51/results/Bcells/', pattern='*pathways_mSigDb_c8_*', full.names = T, ignore.case = F)

# load the files
ldfData = lapply(lFiles, function(x) read.csv(x, row.names=1))
names(ldfData) = lFiles
sapply(ldfData, nrow)
# put everything in one order by row names
rn = rownames(ldfData[[1]])
head(rn)

ldfData = lapply(ldfData, function(df){
  df = df[rn,]
})

# get the upregulated/downregulated in order
ldfData.up = ldfData[grepl(lFiles, pattern = 'upregulated')]
ldfData.down = ldfData[grepl(lFiles, pattern = 'downregulated')]

## set the names for each contrast
sn = gsub('dataID51/results/Bcells//(.+vs.+)_up.+', '\\1', names(ldfData.up))
names(ldfData.up) = sn

sn = gsub('dataID51/results/Bcells//(.+vs.+)_down.+', '\\1', names(ldfData.down))
names(ldfData.down) = sn

## create a table/matrix of p-values
mMerged.up = sapply(ldfData.up, function(x) x$p.val)
rownames(mMerged.up) = rownames(ldfData.up[[1]])
## create similar table for downregulated
mMerged.down = sapply(ldfData.down, function(x) x$p.val)
rownames(mMerged.down) = rownames(ldfData.down[[1]])

# sanity check
identical(rownames(mMerged.up), rownames(mMerged.down))
colnames(mMerged.up) = paste(colnames(mMerged.up), 'up', sep='-')
colnames(mMerged.down) = paste(colnames(mMerged.down), 'down', sep='-')

mMerged.c8 = cbind(mMerged.up, mMerged.down)
# reorder the columns
colnames(mMerged.c8)
o = c(1, 3, 2, 4)
## sanity check
matrix(colnames(mMerged.c8)[o], ncol = 2, byrow = T)
colnames(mMerged.c8)[o]
mMerged.c8 = mMerged.c8[,o]

# remove na sections
dim(mMerged.c8)
mMerged.c8 = na.omit(mMerged.c8)
dim(mMerged.c8)
head(mMerged.c8)

### create a binary matrix based on cutoffs
getBinaryMatrix = function(mat, cutoff=0.01){
  mat2 = apply(mat, 2, function(x) round(x, 3) <= cutoff)
}

mMerged.c8.bin = getBinaryMatrix(mMerged.c8)

## group this matrix into combinations, permutations with repetition
mMerged.c8.bin.grp = mMerged.c8.bin
set.seed(123)
dm = dist(mMerged.c8.bin.grp, method='binary')
hc = hclust(dm)
plot(hc, labels=F)
# cut the tree at the bottom to create groups
cp = cutree(hc, h = 0.2)
# sanity checks
table(cp)
length(cp)
length(unique(cp))

mMerged.c8.bin.grp = cbind(mMerged.c8.bin.grp, cp)

### print and observe this table and select the groups you are interested in
temp = mMerged.c8.bin.grp
temp = (temp[!duplicated(cp),])
temp2 = cbind(temp, table(cp))
rownames(temp2) = NULL
print(temp2)

## map these names to the cp
groups = cp
sig.pvals = rowSums(mMerged.c8.bin)
dfMerged.c8 = data.frame(round(mMerged.c8, 3), sig.pvals, groups, DB='msigdb-c8')
str(dfMerged.c8)
head(dfMerged.c8)
tail(dfMerged.c8)

########
write.csv(dfMerged.c8, file='dataID51/results/Bcells/gsea_msigdb_c8_merged.xls')

######################################################
######### hallmark
rm(lFiles)
lFiles = list.files('dataID51/results/Bcells/', pattern='*pathways_mSigDb_hm_*', full.names = T, ignore.case = F)

# load the files
ldfData = lapply(lFiles, function(x) read.csv(x, row.names=1))
names(ldfData) = lFiles
sapply(ldfData, nrow)
# put everything in one order by row names
rn = rownames(ldfData[[1]])
head(rn)

ldfData = lapply(ldfData, function(df){
  df = df[rn,]
})

# get the upregulated/downregulated in order
ldfData.up = ldfData[grepl(lFiles, pattern = 'upregulated')]
ldfData.down = ldfData[grepl(lFiles, pattern = 'downregulated')]

## set the names for each contrast
sn = gsub('dataID51/results/Bcells//(.+vs.+)_up.+', '\\1', names(ldfData.up))
names(ldfData.up) = sn

sn = gsub('dataID51/results/Bcells//(.+vs.+)_down.+', '\\1', names(ldfData.down))
names(ldfData.down) = sn

## create a table/matrix of p-values
mMerged.up = sapply(ldfData.up, function(x) x$p.val)
rownames(mMerged.up) = rownames(ldfData.up[[1]])
## create similar table for downregulated
mMerged.down = sapply(ldfData.down, function(x) x$p.val)
rownames(mMerged.down) = rownames(ldfData.down[[1]])

# sanity check
identical(rownames(mMerged.up), rownames(mMerged.down))
colnames(mMerged.up) = paste(colnames(mMerged.up), 'up', sep='-')
colnames(mMerged.down) = paste(colnames(mMerged.down), 'down', sep='-')

mMerged.hm = cbind(mMerged.up, mMerged.down)
# reorder the columns
colnames(mMerged.hm)
o = c(1, 3, 2, 4)
## sanity check
matrix(colnames(mMerged.hm)[o], ncol = 2, byrow = T)
colnames(mMerged.hm)[o]
mMerged.hm = mMerged.hm[,o]

# remove na sections
dim(mMerged.hm)
mMerged.hm = na.omit(mMerged.hm)
dim(mMerged.hm)
head(mMerged.hm)

### create a binary matrix based on cutoffs
getBinaryMatrix = function(mat, cutoff=0.01){
  mat2 = apply(mat, 2, function(x) round(x, 3) <= cutoff)
}

mMerged.hm.bin = getBinaryMatrix(mMerged.hm)

## group this matrix into combinations, permutations with repetition
mMerged.hm.bin.grp = mMerged.hm.bin
set.seed(123)
dm = dist(mMerged.hm.bin.grp, method='binary')
hc = hclust(dm)
plot(hc, labels=F)
# cut the tree at the bottom to create groups
cp = cutree(hc, h = 0.2)
# sanity checks
table(cp)
length(cp)
length(unique(cp))

mMerged.hm.bin.grp = cbind(mMerged.hm.bin.grp, cp)

### print and observe this table and select the groups you are interested in
temp = mMerged.hm.bin.grp
temp = (temp[!duplicated(cp),])
temp2 = cbind(temp, table(cp))
rownames(temp2) = NULL
print(temp2)

## map these names to the cp
groups = cp
sig.pvals = rowSums(mMerged.hm.bin)
dfMerged.hm = data.frame(round(mMerged.hm, 3), sig.pvals, groups, DB='msigdb-hm')
str(dfMerged.hm)
head(dfMerged.hm)
tail(dfMerged.hm)

########
write.csv(dfMerged.hm, file='dataID51/results/Bcells/gsea_msigdb_hm_merged.xls')

####################################################
########## merge all the results
####################################################
## merge together into one dataframe
# drop the group with most zeros
table(dfMerged.c2$groups)
t = rowSums(mMerged.c2.bin)
table(t, dfMerged.c2$groups)
dfMerged.c2.sub = dfMerged.c2[dfMerged.c2$groups != 4,]

table(dfMerged.c3$groups)
t = rowSums(mMerged.c3.bin)
table(t, dfMerged.c3$groups)
dfMerged.c3.sub = dfMerged.c3[dfMerged.c3$groups != 2,]

table(dfMerged.c5$groups)
t = rowSums(mMerged.c5.bin)
table(t, dfMerged.c5$groups)
dfMerged.c5.sub = dfMerged.c5[dfMerged.c5$groups != 4,]

table(dfMerged.c7$groups)
t = rowSums(mMerged.c7.bin)
table(t, dfMerged.c7$groups)
dfMerged.c7.sub = dfMerged.c7[dfMerged.c7$groups != 4,]

table(dfMerged.c8$groups)
t = rowSums(mMerged.c8.bin)
table(t, dfMerged.c8$groups)
dfMerged.c8.sub = dfMerged.c8[dfMerged.c8$groups != 3,]

table(dfMerged.hm$groups)
t = rowSums(mMerged.hm.bin)
table(t, dfMerged.hm$groups)
dfMerged.hm.sub = dfMerged.hm[dfMerged.hm$groups != 1,]

dfMerged = rbind(dfMerged.c2.sub, dfMerged.c3.sub, dfMerged.c5.sub,
                 dfMerged.c7.sub, dfMerged.c8.sub, dfMerged.hm.sub)
dfMerged = droplevels.data.frame(dfMerged)
dim(dfMerged)
str(dfMerged)

write.csv(dfMerged, file='dataID51/results/Bcells/gsea_msigdb_significant_6db_merged.xls')

### heatmaps
### just for a quick visual check, do not use for results
df = dfMerged
head(df)
dim(df)
mMat = as.matrix(df[,c(1:4)])
head(mMat)
mMat = -1*log(mMat+1e-16)
g1 = df[,'groups']
g1 = factor(as.character(g1))
levels(g1)
g2 = df[,'DB']
g2 = factor(as.character(g2))
levels(g2)

ann = data.frame(DB=g2 )
#ann = data.frame(Group=g1 )
range(mMat)
quantile(as.vector(mMat), 0:20/20)
#mMat[mMat < 15] = 15 
mMat[mMat > 15] = 15

library(NMF)
library(RColorBrewer)

aheatmap(mMat, annRow = ann, scale = 'none', Rowv = order(g2:g1), Colv=NA, cexRow=5, cexCol = 0.6, #labCol=c('C2vC1', 'K1vC1', 'K2vC2', 'K2vK1'), 
         col=c('white', brewer.pal(9, 'YlOrRd')))

pdf('dataID51/results/Bcells/gsea_msigdb_significant_merged.pdf')
aheatmap(mMat, annRow = ann, scale = 'none', Rowv = order(g2:g1), Colv=NA, cexRow=5, cexCol = 0.6, #labCol=c('C2vC1', 'K1vC1', 'K2vC2', 'K2vK1'), 
         col=c('white', brewer.pal(9, 'YlOrRd')))
dev.off(dev.cur())

f = ann$DB == 'msigdb-c7'
f = grep('t_cell|tcell', rownames(mMat), ignore.case = T)
ann = data.frame(DB=droplevels(g2[f]))
aheatmap(mMat[f,], annRow = ann,  scale = 'none', Rowv = order(g2[f]:g1[f]), Colv=NA, cexRow=5, cexCol = 0.6, #labCol=c('C2vC1', 'K1vC1', 'K2vC2', 'K2vK1'), 
         col=c('white', brewer.pal(9, 'YlOrRd')), main='T cell')


