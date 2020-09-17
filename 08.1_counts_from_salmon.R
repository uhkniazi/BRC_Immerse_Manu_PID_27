# File: 08_counts_from_salmon.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: generate count tables for transcripts from bam files
# Date: 16/9/2020


## set variables and source libraries
source('header.R')

## following workflow to import data 
## https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html#Salmon_with_inferential_replicates
library(tximport)

## create the file list from database
## connect to mysql database to get sample information
library('RMySQL')

##### connect to mysql database to get samples
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)
dbListFields(db, 'File')
# get the query
g_did
q = paste0('select Sample.id as sid, Sample.group1, Sample.group2, Sample.title, File.* from Sample, File
           where (Sample.idData = 49) AND (File.idSample = Sample.id AND File.type like "%fastq%")')
dfSample = dbGetQuery(db, q)
nrow(dfSample)
dfSample
# close connection after getting data
dbDisconnect(db)
# remove any whitespace from the names
dfSample$name = gsub(" ", "", dfSample$name, fixed = T)
dfSample$title = gsub(" ", "", dfSample$title, fixed = T)
dfSample$group1 = gsub(" ", "", dfSample$group1, fixed = T)
i = grep('_R1_', dfSample$name)
dfSample = dfSample[i,]

# create a new file path variable 
dfSample$fp = paste0('dataExternal/Salmon/', dfSample$name, '_quant/quant.sf')

csFiles = dfSample$fp
names(csFiles) = dfSample$sid
all(file.exists(csFiles))

oTxImport = tximport(csFiles, type='salmon', txOut=T)

# sanity checks
identical(colnames(oTxImport$counts), as.character(dfSample$sid))

## save the summarized experiment object
setwd(gcswd)
n = make.names(paste('TxImport transcript count object did 49 manu rds'))
n2 = paste0('~/Data/MetaData/', n)
save(oTxImport, file=n2)

## comment out after first time execution
library('RMySQL')
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)
dbListFields(db, 'MetaFile')
df = data.frame(idData=g_did, name=n, type='rds', location='~/Data/MetaData/',
                comment='Tx Import object of transcript counts using salmon from manu immerse human project')
# dbWriteTable(db, name = 'MetaFile', value=df, append=T, row.names=F)
dbDisconnect(db)
