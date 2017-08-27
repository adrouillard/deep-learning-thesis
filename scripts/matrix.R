# Author : Drouillard Antoine
# Date : 10/04/2017
# Create a input matrix
# Deploy packages-----------------------------------------------------------------------------------------------------------------------
library("DBI")
library("RMySQL")
library("proto")
library("gsubfn")
library("mxnet")
library("doParallel")
library("snow")
library("foreach")
# Connect to database and import data---------------------------------------------------------------------------------------------------
mydb = dbConnect(MySQL(), user='root', password='Macrobio21!', dbname='mass_spectra', host = '127.0.0.1')
#peak = dbReadTable(mydb, "peaks")
mol = dbReadTable(mydb, "molecule")



#import a selected numkber of peaks ----------------------------------------------------------------------------------------------------
sel = dbSendQuery(mydb, "select * from peaks where mz < 300")
res_sel= dbFetch(sel, n = -1)

#Progress bar---------------------------------------------------------------------------------------------------------------------------
pb <- txtProgressBar(min = 0, max = dim(table(res_sel$idmol)), style = 3)

#Create matrices of labeled values bin 0.1----------------------------------------------------------------------------------------------
accession = unique(res_sel$idmol)
digit = 3000
beg <- proc.time()
M = matrix(data = 0, nrow = dim(table(res_sel$idmol)), ncol = digit +2)
for (i in 1:dim(M)[1]) {
  peaks = subset(res_sel,res_sel$idmol ==accession[i])
  molecule = subset(mol,mol$access == accession[i])
  M[i,2] = as.numeric(molecule$mass)
  M[i,1] = molecule$access
  for (j in 1:dim(peaks)[1]) {
    if (peaks$mz[j] < 0.1){
      M[i,3] = as.numeric(M[i,3]) + peaks$int[j]
    }else {
      M[i,round(peaks$mz[j],1)*10+2] =  as.numeric(M[i,round(peaks$mz[j],1)*10+2]) + peaks$int[j]
    }
  }
  setTxtProgressBar(pb, i)
}

print(proc.time() - beg)
#Create a vector of labels and matrices of values --------------------------------------------------------------------------------------
dat = data.frame(M)
dat = subset(dat, (!is.na(dat$X2)))
write.csv(dat,"martix_bin01.csv")

#Create matrices of labeled values bin 0.5----------------------------------------------------------------------------------------------
accession = unique(res_sel$idmol)
digit = 600
beg <- proc.time()
M = matrix(data = 0, nrow = dim(table(res_sel$idmol)), ncol = digit +2)
for (i in 1:dim(M)[1]) {
  peaks = subset(res_sel,res_sel$idmol ==accession[i])
  molecule = subset(mol,mol$access == accession[i])
  M[i,2] = as.numeric(molecule$mass)
  M[i,1] = molecule$access
  for (j in 1:dim(peaks)[1]) {
    if (peaks$mz[j] < 0.1){
      M[i,3] = as.numeric(M[i,3]) + peaks$int[j]
    }else {
      M[i,ceiling(peaks$mz[j]/0.5)+2] =  as.numeric(M[i,ceiling(peaks$mz[j]/0.5)+2]) + peaks$int[j]
    }
  }
  setTxtProgressBar(pb, i)
}

print(proc.time() - beg)
#Create a vector of labels and matrices of values --------------------------------------------------------------------------------------
dat = data.frame(M)
dat = subset(dat, (!is.na(dat$X2)))
write.csv(dat,"martix_bin05.csv")

#Create matrices of labeled values bin 1----------------------------------------------------------------------------------------------
accession = unique(res_sel$idmol)
digit = 300
beg <- proc.time()
M = matrix(data = 0, nrow = dim(table(res_sel$idmol)), ncol = digit +2)
for (i in 1:dim(M)[1]) {
  peaks = subset(res_sel,res_sel$idmol ==accession[i])
  molecule = subset(mol,mol$access == accession[i])
  M[i,2] = paste("'",molecule$bitstring,"'",sep = '')
  M[i,1] = molecule$access
  for (j in 1:dim(peaks)[1]) {
    if (peaks$mz[j] < 0.1){
      M[i,3] = as.numeric(M[i,3]) + peaks$int[j]
    }else {
      M[i,round(peaks$mz[j])+2] =  as.numeric(M[i,round(peaks$mz[j])+2]) + peaks$int[j]
    }
  }
  setTxtProgressBar(pb, i)
}

print(proc.time() - beg)
#Create a vector of labels and matrices of values --------------------------------------------------------------------------------------
dat = data.frame(M)
dat = subset(dat, (!is.na(dat$X2)))
write.csv(dat,"martix_bin1.csv")


#Create matrices of labeled values bin 2------------------------------------------------------------------------------------------------
accession = unique(res_sel$idmol)
digit = 150
beg <- proc.time()
M = matrix(data = 0, nrow = dim(table(res_sel$idmol)), ncol = digit +2)
for (i in 1:dim(M)[1]) {
  peaks = subset(res_sel,res_sel$idmol ==accession[i])
  molecule = subset(mol,mol$access == accession[i])
  M[i,2] = as.numeric(molecule$mass)
  M[i,1] = molecule$access
  for (j in 1:dim(peaks)[1]) {
    if (peaks$mz[j] < 0.1){
      M[i,3] = as.numeric(M[i,3]) + peaks$int[j]
    }else {
        M[i,ceiling(peaks$mz[j]/2)+2] =  as.numeric(M[i,ceiling(peaks$mz[j]/2)+2]) + peaks$int[j]
    }
  }
  setTxtProgressBar(pb, i)
}

print(proc.time() - beg)
#Create a vector of labels and matrices of values --------------------------------------------------------------------------------------
dat = data.frame(M)
dat = subset(dat, (!is.na(dat$X2)))
write.csv(dat,"martix_bin2.csv")

#Create matrices of labeled values bin 2------------------------------------------------------------------------------------------------
accession = unique(res_sel$idmol)
digit = 60
beg <- proc.time()
M = matrix(data = 0, nrow = dim(table(res_sel$idmol)), ncol = digit +2)
for (i in 1:dim(M)[1]) {
  peaks = subset(res_sel,res_sel$idmol ==accession[i])
  molecule = subset(mol,mol$access == accession[i])
  M[i,2] = as.numeric(molecule$mass)
  M[i,1] = molecule$access
  for (j in 1:dim(peaks)[1]) {
    if (peaks$mz[j] < 0.1){
      M[i,3] = as.numeric(M[i,3]) + peaks$int[j]
    }else {
      M[i,ceiling(peaks$mz[j]/5)+2] =  as.numeric(M[i,ceiling(peaks$mz[j]/5)+2]) + peaks$int[j]
    }
  }
  setTxtProgressBar(pb, i)
}

print(proc.time() - beg)
#Create a vector of labels and matrices of values --------------------------------------------------------------------------------------
dat = data.frame(M)
dat = subset(dat, (!is.na(dat$X2)))
write.csv(dat,"martix_bin5.csv")

