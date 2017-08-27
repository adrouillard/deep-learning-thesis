# Author : Drouillard Antoine
# Date : 05/04/2017
# Import and prepare date for the Deep Learning Application

# Deploy packages-----------------------------------------------------------------------------------------------------------------------
library("DBI", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.3")
library("RMySQL", lib.loc="/usr/lib/R/site-library")
library("gsubfn", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.3")
library("imager", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.3")
# Connect to database and import data---------------------------------------------------------------------------------------------------
mydb = dbConnect(MySQL(), user='root', password='pass', dbname='mass_spectra', host = '127.0.0.1')
peak = dbReadTable(mydb, "peaks")
mol = dbReadTable(mydb, "molecule")

#import data from one choosen molecule--------------------------------------------------------------------------------------------------
acc = 'AU100701'
test = fn$dbSendQuery(mydb, "select * from peaks where idmol ='$acc'")
result = dbFetch(test)

#draw plots-----------------------------------------------------------------------------------------------------------------------------
plot(result$mz,result$int, pch = 3, type = "h" )
plot(result$mz,log(result$int), pch = 3, type = "h" )
plot(result$mz,result$int, pch = 3)
plot(result$mz,log(result$int), pch = 3)

dim(plot(result$mz,log(result$int), pch = 3))
