# File: 01_createDbEntry.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: list the samples and create appropriate db entries
# Date: 7/9/2020


## set variables and source libraries
source('header.R')

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

setwd('dataExternal/remote/raw/')

# list the files
cvFiles = list.files(pattern = 'fastq.gz', recursive = T)

# each sample has 2 files 
fSplit = gsub('_R[1|2]_.+', '', cvFiles)

lFiles = split(cvFiles, fSplit)

setwd(gcswd)
## load the metadata file
dfMeta = read.csv('dataExternal/cardiac_metadata_RNA.csv', header=T,
                  stringsAsFactors = F, row.names = 1)
str(dfMeta)
cn = colnames(dfMeta)
# remove any white space
for (i in seq_along(1:ncol(dfMeta))) dfMeta[,cn[i]] = gsub(' ', '', dfMeta[,cn[i]])
dim(dfMeta)
# sanity check
table(dfMeta$Sample_ID %in% names(lFiles))
table(names(lFiles) %in% dfMeta$Sample_ID)
## keep only the matching files
f = names(lFiles) %in% dfMeta$Sample_ID
lFiles = lFiles[f]
table(names(lFiles) %in% dfMeta$Sample_ID)
# order the files list and table in same sequence
identical(names(lFiles), dfMeta$Sample_ID)

## structure of the data and covariates
str(dfMeta)
xtabs( ~ Timepoint + study_id, data=dfMeta)
xtabs( ~ Timepoint + outcome_icu_survival, data=dfMeta)
xtabs( ~ Timepoint + RNA_extraction_run_num, data=dfMeta)

## create the entry for samples
cSampleCol

dfSamples = data.frame(idProject=g_pid, idData=g_did, title=dfMeta$Sample_ID,
                       location='rosalind scratch and a copy with manu',
                       description= paste('age', dfMeta$age_years,
                                          'height', dfMeta$height_cms,
                                          'sex', dfMeta$sex,
                                          'group1 is Time point',
                                          'group2 is Lane',
                                          sep=';'),
                       group1 = dfMeta$Timepoint,
                       group2= 0,
                       group3=0)
# write this data to the database
rownames(dfSamples) = NULL

### NOTE: Do not execute this anymore as entry created
# write this table to database
#dbWriteTable(db, name='Sample', value=dfSamples, append=T, row.names=F)
# get this table again from database with ids added
g_did
dfSamples = dbGetQuery(db, paste('select * from Sample where Sample.idData = 49;'))

# create entries for these files in the database
dbListTables(db)
cn = dbListFields(db, 'File')[-1]
cn
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

# sanity check
fn = gsub('_R[1|2]_.+', '', dfFiles$name)
sapply(dfSamples$id, function(x) {
  table(fn[dfFiles$idSample == x] %in% dfSamples$title[dfSamples$id == x])
})
table(dfFiles$name %in% c(dfMeta$File_name_R1, dfMeta$File_name_R2))

# write this table to database
## note: do not execute as it is already done
#dbWriteTable(db, name='File', value=dfFiles, append=T, row.names=F)

dbDisconnect(db)
