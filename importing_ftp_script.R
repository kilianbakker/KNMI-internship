# set library paths:
libpath = "/usr/people/bakker/R/x86_64-redhat-linux-gnu-library/3.4"
.libPaths(c(.libPaths()))
library(tidyverse)
library(RCurl)

# The variables/constants
monthNumbers <- c("01","02","03","04","05","06","07","08","09","10","11","12")
indir          <- "ftp://csr_knmi:C5r-Knm1!@bbc.knmi.nl/MEMBERS/knmi/bsrn/"

yearNumber <- 2016
for (i in 4:12){
  url <- paste0(indir,yearNumber,"/",monthNumbers[i],"/")
  filenames <- getURL(url, ftp.use.epsv = T ,dirlistonly = T, verbose = T) 
  filenames <- paste(url, strsplit(filenames, "\r*\n")[[1]], sep = "")
  filenames <- filenames[gsub("_[0-9]{8}.csv.gz","",gsub("ftp://csr_knmi:C5r-Knm1!@bbc.knmi.nl/MEMBERS/knmi/bsrn/[0-9]{4}/[0-9]{2}/","",filenames)) == "BSRN"]
  for (j in 1:length(filenames)){
    fileDate <- gsub(".csv.gz","",gsub("ftp://csr_knmi:C5r-Knm1!@bbc.knmi.nl/MEMBERS/knmi/bsrn/[0-9]{4}/[0-9]{2}/BSRN_","",filenames[j]))
    system(paste0("wget ",filenames[j]," -O /nobackup/users/bakker/Data/observations_min/",fileDate,".csv.gz"))
    system(paste0("gunzip /nobackup/users/bakker/Data/observations_min/",fileDate,".csv.gz"))
  }
}

yearNumber <- 2017
for (i in 1:12){
  url <- paste0(indir,yearNumber,"/",monthNumbers[i],"/")
  filenames <- getURL(url, ftp.use.epsv = T ,dirlistonly = T, verbose = T) 
  filenames <- paste(url, strsplit(filenames, "\r*\n")[[1]], sep = "")
  filenames <- filenames[gsub("_[0-9]{8}.csv.gz","",gsub("ftp://csr_knmi:C5r-Knm1!@bbc.knmi.nl/MEMBERS/knmi/bsrn/[0-9]{4}/[0-9]{2}/","",filenames)) == "BSRN"]
  for (j in 1:length(filenames)){
    fileDate <- gsub(".csv.gz","",gsub("ftp://csr_knmi:C5r-Knm1!@bbc.knmi.nl/MEMBERS/knmi/bsrn/[0-9]{4}/[0-9]{2}/BSRN_","",filenames[j]))
    system(paste0("wget ",filenames[j]," -O /nobackup/users/bakker/Data/observations_min/",fileDate,".csv.gz"))
    system(paste0("gunzip /nobackup/users/bakker/Data/observations_min/",fileDate,".csv.gz"))
  }
}

yearNumber <- 2018
for (i in 1:7){
  url <- paste0(indir,yearNumber,"/",monthNumbers[i],"/")
  filenames <- getURL(url, ftp.use.epsv = T ,dirlistonly = T, verbose = T) 
  filenames <- paste(url, strsplit(filenames, "\r*\n")[[1]], sep = "")
  filenames <- filenames[gsub("_[0-9]{8}.csv.gz","",gsub("ftp://csr_knmi:C5r-Knm1!@bbc.knmi.nl/MEMBERS/knmi/bsrn/[0-9]{4}/[0-9]{2}/","",filenames)) == "BSRN"]
  for (j in 1:length(filenames)){
    fileDate <- gsub(".csv.gz","",gsub("ftp://csr_knmi:C5r-Knm1!@bbc.knmi.nl/MEMBERS/knmi/bsrn/[0-9]{4}/[0-9]{2}/BSRN_","",filenames[j]))
    system(paste0("wget ",filenames[j]," -O /nobackup/users/bakker/Data/observations_min/",fileDate,".csv.gz"))
    system(paste0("gunzip /nobackup/users/bakker/Data/observations_min/",fileDate,".csv.gz"))
  }
}