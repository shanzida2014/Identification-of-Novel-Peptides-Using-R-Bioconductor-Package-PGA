#Title: Making a Fasta file of all novel exons for use with PGA
# Author: Shanzida Jahan Siddique  (ss51@hood.edu)
# Date: created on 6/14/19
# PGA paper http://www.bioconductor.org/packages/release/bioc/vignettes/PGA/inst/doc/PGA.pdf
###########################################

### Install BiocManager. Only needs to be done once.
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

### Install PGA. Also only needs to be done once
BiocManager::install("PGA")

library(PGA)

###Read in novel exons and create Fasta File
-----------------------------------------------
setwd("~/Capstone Dropbox/Shanzida Siddique/Shanzida/")
#install.packages("seqRFLP")
library("seqRFLP")
newExons <- read.csv("~/Documents/Capstone/ShanzidaSiddiqueCapstoneReformat.csv")

#Create FASTA headers
newExHeader <- paste0(">",newExons$putativeExon,"_",newExons$strand,"::::",newExons$GENE_ID)
#Create vector of novel exon sequences (DNA sequence)
sequences = newEx$sequence....revComp.
#install.packages("seqRFLP") for converting peptide sequences to FASTA format
#install library seqRFLP for that 
library("seqRFLP")
names=newExHeader 
sequences = newEx$sequence....revComp.S
#convert sequences and header to the dataframe
df <- data.frame(names,sequences)
#convert dataframe in to FASTA format
df.fasta = dataframe2fas(df, file="df.fasta")
write.fasta(df.fasta, file="~/Documents/Capstone/df.fasta")
##Followed this tutorialof PGA
#--------------------------------------
#https://bioconductor.org/packages/release/bioc/vignettes/PGA/inst/doc/PGA.pdf
#STEP 1:Construction of customized protein database based on RNA-Seq data
--------------------------------------------------------------
###STEP 1.1:Based on the result from analysis of RNA-Seq data with a reference genome
#---------------------------------------------------------------------------------------
#I didn't follow this because didn't have
#1.VCF format file which contains SNV or INDEL information;
#2. a BED format file which contains splice junctions information;
#3. a GTF format file which contains novel transcripts information.
# STEP1.2:Based on the result from de novo assembly of RNASeq data without a reference genome
#------------------------------------------------------------------
# I followed this step insted 1.1.
library(PGA)


transcript_seq_file <- "/home/ss51/Capstone/df.fasta"

outdb <- createProDB4DenovoRNASeq(infa=transcript_seq_file,outfile_name = "denovo")
# 
# ....... Proteomic database construction ! ............
# 
# Write protein sequences to file: denovo_txFinder.fasta
# 
#  cat(outdb,"\n")
# 
# denovo_txFinder.fasta 

#msfile <- "/home/ss51/Capstone/Alz_P01_A01_097_26Apr12_Roc_12-03-15.mgf"
msfile <- "/BIODATA/PRIDE_DATA/PXD006537/Alz_P07_D06_042_16Apr12_Roc_12-03-26.mgf"
outdir <- "/home/ss51/Capstone/result/Alz_P07_D06_042_16Apr12_Roc_12-03-26"
idfile <- runTandem(spectra = msfile, fasta = outdb,outdir ="/home/ss51/Capstone/result",cpu = 6,enzyme = "[KR]|[X]", varmod = "15.994915@M",itol = 0.05,fixmod ="57.021464@C", tol = 10, tolu = "ppm",itolu = "Daltons", miss = 2, maxCharge = 8, ti = FALSE)

# 2019-07-12 14:06:00 
# 
# Loading spectra
# 
# (mgf).... loaded.
# 
# Spectra matching criteria = 6114
# 
# Starting threads ...... started.
# 
# Computing models:
#   
#   testin
# 
# sequences modelled = 6 ks
# 
# Model refinement:
#   
#   Merging results:
#   
#   from 2.3.4.5.6.
# 
# 
# Creating report:
#   
#   initial calculations  ..... done.
# 
# sorting  ..... done.
# 
# finding repeats ... done.
# 
# evaluating results ... done.
# 
# calculating expectations ... done.
# 
# writing results ..... done.
# 
# 
# Valid models = 2
# 
# Unique models = 2
# 
# Estimated false positives = 0 +/- 1
## Post-processing
#--------------------------------------------------------------
##After the MS/MS data searching, the function parserGear can be used to parse the search result. It calculates the q-value for each peptide spectrum
#matches (PSMs) and then utilizes the Occam’s razor approach [15] to deal
#with degenerated wild peptides by finding a minimum subset of proteins that
#covered all of the identified wild peptides.
parserGear(file = idfile, db = outdb,decoyPrefix="#REV#",xmx=1,thread=8,outdir ="/home/ss51/Capstone/PXD010271/result/parser_outdir")

#It exports some tab-delimited files containing the peptide identification result
#and protein identification result. The annotated spectra for the identified novel
#peptides which pass the threshold are exported.

##HTML-based report generation
------------------------------------------
##The results are then summarised and compiled into an interactive HTML report.
reportGear(parser_dir = "/home/ss51/Capstone/PRIDE_BRAIN/PXD006537/result/parser_outdir",tab_dir = "/home/ss51/Capstone/outfile_name/denovo",report_dir="/home/ss51/Capstone/PRIDE_BRAIN/PXD006537/result/report")
#SNV(DB) didn't exist!
##After the analysis has completed, the file ‘index.html’ in the output directory
#can be opened in a web browser to access report generated. In general, this
#report will show the identification result for four kinds of novel peptides, such
#as SNV-caused peptides, INDEL-caused peptides, alternative splicing caused
#peptides and novel transcripts codeing peptides.
#If the RefSeq annotation is used in the RNA-Seq data analysis, the report will
#show the gene name for each protein.
##I didn't get the result instead got :
#SNV(DB) didn't exist!
#Error in 1:dim(data)[1] : argument of length 0

## Search in a control :

msfile <- "/home/ss51/Capstone/Alz_P07_D06_042_16Apr12_Roc_12-03-26.mgf"
idfile <- runTandem(spectra = msfile, fasta = outdb,outdir ="/home/ss51/Capstone/PRIDE_BRAIN/PXD0065372019-11-05 04:16:21 yme = "[KR]|[X]", varmod = "15.994915@M",itol = 0.05,fixmod ="57.021464@C", tol = 10Loading spectraitolu = "Daltons", miss = 2, maxCharge = 8, ti = FALSE)
 (mgf).... loaded.
Spectra matching criteria = 5944
Starting threads ...... started.
Computing models:
	testin
		sequences modelled = 6 ks
Model refinement:
Merging results:
	from 23456

Creating report:
	initial calculations  ..... done.
	sorting  ..... done.
	finding repeats . done.
	evaluating results . done.
	calculating expectations . done.
	writing results ..... done.

Valid models = 1
Unique models = 1
Estimated false positives = 0 +/- 1

list.files("/BIODATA/PRIDE_DATA/PXD006537/")
data <- read.delim("/BIODATA/PRIDE_DATA/PXD006537/PNNL_Dataset_List.txt")
data <- as.data.frame(data)
samples <- data$SAMPLE
sampletype <- sapply(strsplit(samples,"_"),"[[",1)
controlSamples <- which(sampletype == "Control")
controlFiles <- data$DatasetName[controlSamples]
controlFilesMGF <- paste(controlFiles,".mgf",sep="")


