# Author : Drouillard Antoine
# Date : 10/05/2017
# Convert information into bitstring and fill the database

library("DBI")
library("RMySQL")
library("fingerprint")
library("rcdk")

mydb = dbConnect(MySQL(), user='root', password='Macrobio21!', dbname='mass_spectra', host = '127.0.0.1')
mol = dbReadTable(mydb, "molecule")

#Progress bar---------------------------------------------------------------------------------------------------------------------------
pb <- txtProgressBar(min = 1, max = dim(mol)[1], style = 3)
IA <- parse.smiles(mol$smiles, kekulise = TRUE)
#core ----------------------------------------------------------------------------------------------------------------------------------
for (el in 1:dim(mol)[1]){
  if (is.na(mol$bitstring[el])){
    if (!is.na(IA[[el]])){
      finpr <- get.fingerprint(IA[[el]], type ='maccs',verbose=FALSE)
      arr <- array(finpr@bits)  
      bitstr <- array(data = 0, dim = 166)
      for (i in 1:dim(arr)){
        bitstr[arr[i]] = 1 
      }
      str = ''
      for (j in 1:166){
        str = paste(str, bitstr[j], sep='')
      }
      
      #*****Update the database*****
      mol$bitstring[el] = str
      query = paste("UPDATE `mass_spectra`.`molecule` SET `bitstring`=", " '",str, "' WHERE `access`=", " '",mol$access[el], "' ;", sep = '')
      dbSendQuery(mydb, query)
    }
  }  
  setTxtProgressBar(pb, el)
}

#Evaluate the-- repartition of values between all bits-----------------------------------------------------------------------------------
M = matrix(data = 0, nrow = dim(mol)[1] , ncol = 166)
for (x in 1:dim(mol)[1]){
  for (y in 1:166)
  M[x,y] = unlist(strsplit(mol$bitstring[x], ''))[y]
}

for (z in 1:166){
  if (table(M[,z])[1]/43694 < 0.6){
    if (table(M[,z])[1]/43694 >0.4){
      print(z) 
      print(table(M[,z])[1]/43694)
    }
  }
}


