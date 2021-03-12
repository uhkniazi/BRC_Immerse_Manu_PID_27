# File: 01_createDbEntry.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: list the samples and create appropriate db entries
# Date: 12/3/2021


## set variables and source libraries
source('header.R')
setwd('dataID51/')
## connect to mysql database 
library('RMySQL')

##### connect to mysql database to get samples
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)

# sample and file table
dbGetQuery(db, paste('describe Sample;'))
cSampleCol = dbGetQuery(db, paste('describe Sample;'))$Field[-1]

dbGetQuery(db, paste('describe File;'))
cFileCol = dbGetQuery(db, paste('describe File;'))$Field[-1]

setwd(gcswd)
setwd('dataExternal/remote/dataID51/raw/S486/raw/')

# list the files
cvFiles = list.files(pattern = 'fastq.gz', recursive = T)
cvFiles = cvFiles[-c(69, 70)]

# each sample has 2 files 
fSplit = gsub('_R[1|2]_.+', '', cvFiles)

lFiles = split(cvFiles, fSplit)

setwd(gcswd)
setwd('dataID51/')
## load the metadata file
dfMeta = read.csv('dataExternal/metadata.csv', header=T,
                  stringsAsFactors = F)
str(dfMeta)
cn = colnames(dfMeta)
# remove any white space
for (i in seq_along(1:ncol(dfMeta))) dfMeta[,cn[i]] = gsub(' ', '', dfMeta[,cn[i]])
dim(dfMeta)
dfMeta = na.omit(dfMeta)
# sanity check
table(dfMeta$Sample_id %in% names(lFiles))
# sample names are written differently so don't match
table(c(dfMeta$fastq_R1, dfMeta$fastq_R2) %in% cvFiles)
table(cvFiles %in% c(dfMeta$fastq_R1, dfMeta$fastq_R2))
## 2 files mismatch
cvFiles[!(cvFiles %in% c(dfMeta$fastq_R1, dfMeta$fastq_R2))]
# these don't have a time label and were not in the 
# metadata sheet
# [1] "SC06-B_S26_L001_R1_001.fastq.gz" "SC06-B_S26_L001_R2_001.fastq.gz"
## keep only the matching files
cvFiles = cvFiles[(cvFiles %in% c(dfMeta$fastq_R1, dfMeta$fastq_R2))]
fSplit = gsub('_R[1|2]_.+', '', cvFiles)
lFiles = split(cvFiles, fSplit)

fSubject = gsub('_T\\d_(B|T)', '', dfMeta$Sample_id)
dfMeta$fSubject = factor(fSubject)
## structure of the data and covariates
str(dfMeta)
xtabs( ~ fSubject + Cohort_time, data=dfMeta[dfMeta$Cell_pop == 'B',])
xtabs( ~ Sample_id + fSubject, data=dfMeta[dfMeta$Cell_pop == 'B', ])

## create the entry for samples
cSampleCol

dfSamples = data.frame(idProject=g_pid, idData=g_did, title=dfMeta$Sample_id,
                       location='rosalind scratch and a copy with manu',
                       description= paste('group1 is subject id for pairing',
                                          'group2 is interaction of cohort and time',
                                          'group3 is cell population B or T',
                                          sep=';'),
                       group1 = dfMeta$fSubject,
                       group2= dfMeta$Cohort_time,
                       group3= dfMeta$Cell_pop)
# write this data to the database
rownames(dfSamples) = NULL

### NOTE: Do not execute this anymore as entry created
# write this table to database
#dbWriteTable(db, name='Sample', value=dfSamples, append=T, row.names=F)
# get this table again from database with ids added
g_did
dfSamples = dbGetQuery(db, paste('select * from Sample where Sample.idData = 51;'))

# create entries for these files in the database
dbListTables(db)
cn = dbListFields(db, 'File')[-1]
cn
## put file names in same order as sample ids
lFiles = lapply(1:nrow(dfMeta), function(x){
  return(c(dfMeta[x,'fastq_R1'], dfMeta[x,'fastq_R2']))
})
names(lFiles) = as.character(dfMeta$Sample_id)
table(unlist(lFiles) %in% cvFiles)

table(names(lFiles) %in% dfSamples$title)
identical(names(lFiles), dfSamples$title)
names(lFiles) = dfSamples$id

# get the names of the samples
temp = lapply(as.character(dfSamples$id), function(x){
  # get the file names
  df = data.frame(name=lFiles[[x]], type='fastq', idSample=dfSamples[as.character(dfSamples$id) == x, 'id'])
  return(df)
})

dfFiles = do.call(rbind, temp)
rownames(dfFiles) = NULL

table(dfFiles$name %in% c(dfMeta$fastq_R1, dfMeta$fastq_R2))

# write this table to database
## note: do not execute as it is already done
#dbWriteTable(db, name='File', value=dfFiles, append=T, row.names=F)

dbDisconnect(db)
