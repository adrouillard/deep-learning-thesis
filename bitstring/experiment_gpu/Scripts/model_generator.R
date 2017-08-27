# Author : Drouillard Antoine
# Date : 10/04/2017
# Create a input matrix
# Deploy packages-----------------------------------------------------------------------------------------------------------------------
library("DBI")
library("RMySQL")
library("proto")
library("gsubfn")
library("mxnet")
library("grid")
library("gplots")
library("pROC")
# Connect to database and import data---------------------------------------------------------------------------------------------------
mydb = dbConnect(MySQL(), user='root', password='Macrobio21!', dbname='mass_spectra', host = '127.0.0.1')
#peak = dbReadTable(mydb, "peaks")
mol = dbReadTable(mydb, "molecule")
tacc = read.csv("../train_access.csv")
test_acc = tacc[,2]
#initialize result matrix---------------------------------------------------------------------------------------------------------------
#results = matrix(data = NA, nrow =  dim(val1)[1]+dim(val2)[1]+dim(val3)[1] , ncol =  2000 + 4)
results = matrix(data = NA, nrow =  (166-8)*5, ncol =  6000 + 4)
constant_bit = c(1,2,3,4,5,6,7,10,166)
complex_bit = c(99,112, 118,123,126,127,131,132,137,139,140,143,146,147,148,149,150,151,152,153,156,158,162)
seed = c(0,10,50,100)

short = array(data = 1:166, dim = 166)
short = short[!(short %in% constant_bit)]
short= short[!(short %in% complex_bit)]
miss = array(data = 0, dim = 5950)
#import a selected numkber of peaks ----------------------------------------------------------------------------------------------------
sel = dbSendQuery(mydb, "select * from peaks where mz < 300")
res_sel= dbFetch(sel, n = -1)

#Code for GPU and CPU------------------------------------------------------------------------------------------------------------------- 
n.gpu <- 1
device.cpu <- mx.cpu()
device.gpu <- lapply(0:(n.gpu-1), function(i) {
  mx.gpu(i)
})

row = 0
# Main loop testing --------------------------------------------------------------------------------------------------------------------

#"***************Creation of models for the "easy bits"***************
for (bit in short) {
  #create matricies for s----------------------------------------------------------------------------------------------------------
  data = read.table("martix_bin1.csv", sep = ',', header = TRUE, row.names = 1)
  data$X2 <- as.character(data$X2)
  #data = data[,-1]
  for (x in 1:dim(data)[1]){
    data[x,2] = unlist(strsplit(as.character(data[x,2]),''))[bit+1]
  }
  data = subset.data.frame(data, !(is.na(data$X2)))
  data_test = data[test_acc,]
  data_train = data[!(data$X1 %in% test_acc),]
  
  #********Creation of input matrices*************
  Mtrain = data.matrix(data_train[,-1])
  Mtest = data.matrix(data_test[,-1])
  train.x <- Mtrain[,-1]
  train.y <- Mtrain[,1]
  test.x <- Mtest[,-1]
  
  #Clean ---------------------------------------------------------------------------------------------------------------------------------
  remove(data_test,data_train, mol, para, tacc, res_sel, data)
  # Deep NN-------------------------------------------------------------------------------------------------------------------------------
  for (s in seed) {
    row = row +1
    #********Creation of the architechture************
    data <- mx.symbol.Variable("data")
    fc1 <- mx.symbol.FullyConnected(data, name="fc1", num_hidden=300)
    act1 <- mx.symbol.Activation(fc1, name="sig1", act_type="sigmoid")
    drop1 <- mx.symbol.Dropout(act1, p = 0.5)
    fc2 <- mx.symbol.FullyConnected(data, name="fc1", num_hidden=200)
    act2 <- mx.symbol.Activation(fc1, name="sig2", act_type="sigmoid")
    drop2 <- mx.symbol.Dropout(act2, p = 0.5)
    fc3 <- mx.symbol.FullyConnected(drop2, name="fc3", num_hidden=1)
    out <- mx.symbol.LogisticRegressionOutput(fc3, name="sm")
    #Learning function----------------------------------------------------------------------------------------------------------------------
    mx.set.seed(s)
    beg <- proc.time()
    log <- mx.metric.logger$new()
    model <- mx.model.FeedForward.create(out,
                                         X=train.x,
                                         y=train.y,
                                         ctx=device.gpu,
                                         num.round=50,
                                         array.batch.size=200,
                                         learning.rate=5.5e-4,
                                         momentum=0.1, 
                                         eval.metric=mx.metric.rmsle,
                                         batch.end.callback = mx.callback.log.train.metric(100, log),
                                         initializer=mx.init.uniform(0.01)
    )
    t <- (proc.time() - beg)
    
    #testing function----------------------------------------------------------------------------------------------------------------------
    beg2 <- proc.time()
    preds <- predict(model,test.x)
    #print(proc.time() - beg2)
    
    beg2 <- proc.time()
    preds_tr <- predict(model,train.x)
    #print(proc.time() - beg2)
    
    #accuracy testing-----------------------------------------------------------------------------------------------------------------------
    Q2 <- 1-( sum((Mtest[,1]-preds)^2) / sum((Mtest[,1]-mean(Mtest[,1]))^2) )
    Q2
    
    Q2_tr <- 1-( sum((train.y-preds_tr)^2) / sum((train.y-mean(train.y))^2) )
    Q2_tr
    auroc = 1
    try({auroc = auc(Mtest[,1],preds)
    auroc})
    
    auroc_tr = 1
    try({auroc_tr = auc(train.y,preds_tr)
    auroc_tr})
    #Labels---------------------------------------------------------------------------------------------------------------------------------
    bin = "LogR_bin1"
    lab = paste(bin,"_",bit,"_",s, sep = '')
    Glab = paste(lab,"_results.png", sep = "")
    Modellab = paste("model/",lab, sep = "")
    
    #Save model-----------------------------------------------------------------------------------------------------------------------------
    mx.model.save(model = model, Modellab, 50)
    res = c(auroc, auroc_tr,t[1], log$train)
    rep_res = data.frame(round(res[1:3], 4))
    rownames(rep_res)[1:3] <- c("AUC test", "AUC train","training time")
    colnames(rep_res) <- paste("results",lab, sep = ' ')
    
    res = c(lab,res, miss)
    results[row,] = res
    #Graph ---------------------------------------------------------------------------------------------------------------------------------
    png(Glab, width = 1600, height = 1200)
    par(mfrow=c(2,2))
    
    plot(log$train, type = 'l', main = "Learning curve", ylab = "Accuracy")
    
    textplot(rep_res)
    
    try(plot(roc(Mtest[,1],preds),main = "Prediction results Test set", ylab = "True positive rate", xlab = "False positive rate" ))
    
    try(plot(roc(train.y,preds_tr),main = "Prediction results Train set", ylab = "True positive rate", xlab = "False positive rate"))
    dev.off()
    
    remove(preds,preds_tr,obs, obs_tr, log, model)
  }
}

#"***************Creation of models for the "complex bits"***************
for (bit in complex_bit) {
  #create matricies for analysis----------------------------------------------------------------------------------------------------------
  data = read.table("martix_bin1.csv", sep = ',', header = TRUE, row.names = 1)
  data$X2 <- as.character(data$X2)
  #data = data[,-1]
  for (x in 1:dim(data)[1]){
    data[x,2] = unlist(strsplit(as.character(data[x,2]),''))[bit+1]
  }
  data = subset.data.frame(data, !(is.na(data$X2)))
  data_test = data[test_acc,]
  data_train = data[!(data$X1 %in% test_acc),]
  
  
  Mtrain = data.matrix(data_train[,-1])
  Mtest = data.matrix(data_test[,-1])
  train.x <- Mtrain[,-1]
  train.y <- Mtrain[,1]
  test.x <- Mtest[,-1]
  
  #Clean ---------------------------------------------------------------------------------------------------------------------------------
  remove(data_test,data_train, mol, para, tacc, res_sel, data)
  # Deep NN-------------------------------------------------------------------------------------------------------------------------------
  
  row = row +1
  data <- mx.symbol.Variable("data")
  fc1 <- mx.symbol.FullyConnected(data, name="fc1", num_hidden=300)
  act1 <- mx.symbol.Activation(fc1, name="sig1", act_type="sigmoid")
  drop1 <- mx.symbol.Dropout(act1, p = 0.5)
  fc2 <- mx.symbol.FullyConnected(data, name="fc1", num_hidden=200)
  act2 <- mx.symbol.Activation(fc1, name="sig2", act_type="sigmoid")
  drop2 <- mx.symbol.Dropout(act2, p = 0.5)
  fc3 <- mx.symbol.FullyConnected(drop2, name="fc3", num_hidden=1)
  out <- mx.symbol.LogisticRegressionOutput(fc3, name="sm")
  
  #Learning function----------------------------------------------------------------------------------------------------------------------
  mx.set.seed(0)
  beg <- proc.time()
  log <- mx.metric.logger$new()
  model <- mx.model.FeedForward.create(out,
                                       X=train.x,
                                       y=train.y,
                                       ctx=device.gpu,
                                       num.round=6000,
                                       array.batch.size=200,
                                       learning.rate=5.5e-4,
                                       momentum=0.1, 
                                       eval.metric=mx.metric.rmsle,
                                       batch.end.callback = mx.callback.log.train.metric(100, log),
                                       initializer=mx.init.uniform(0.01)
  )
  t <- (proc.time() - beg)
  
  #testing function----------------------------------------------------------------------------------------------------------------------
  beg2 <- proc.time()
  preds <- predict(model,test.x)
  #print(proc.time() - beg2)
  
  beg2 <- proc.time()
  preds_tr <- predict(model,train.x)
  #print(proc.time() - beg2)
  
  #accuracy testing-----------------------------------------------------------------------------------------------------------------------
  Q2 <- 1-( sum((Mtest[,1]-preds)^2) / sum((Mtest[,1]-mean(Mtest[,1]))^2) )
  Q2
  
  Q2_tr <- 1-( sum((train.y-preds_tr)^2) / sum((train.y-mean(train.y))^2) )
  Q2_tr
  auroc = 1
  try({auroc = auc(Mtest[,1],preds)
  auroc})
  
  auroc_tr = 1
  try({auroc_tr = auc(train.y,preds_tr)
  auroc_tr})
  #Labels---------------------------------------------------------------------------------------------------------------------------------
  bin = "LogR_bin1"
  lab = paste(bin,"_",bit,"_",0, sep = '')
  Glab = paste(lab,"_results.png", sep = "")
  Modellab = paste("model/",lab, sep = "")
  
  #Save model-----------------------------------------------------------------------------------------------------------------------------
  mx.model.save(model, Modellab,1)
  
  #save results---------------------------------------------------------------------------------------------------------------------------
  res = c(auroc, auroc_tr,t[1], log$train)
  rep_res = data.frame(round(res[1:3], 4))
  rownames(rep_res)[1:3] <- c("AUC test", "AUC train","training time")
  colnames(rep_res) <- paste("results",lab, sep = ' ')
  
  res = c(lab,res)
  results[row,] = res
  #Graph ---------------------------------------------------------------------------------------------------------------------------------
  png(Glab, width = 1600, height = 1200)
  par(mfrow=c(2,2))
  
  plot(log$train, type = 'l', main = "Learning curve", ylab = "Accuracy")
  
  textplot(rep_res)
  
  try(plot(roc(Mtest[,1],preds),main = "Prediction results Test set", ylab = "True positive rate", xlab = "False positive rate" ))
  
  try(plot(roc(train.y,preds_tr),main = "Prediction results Train set", ylab = "True positive rate", xlab = "False positive rate"))
  
  dev.off()
  
  remove(preds,preds_tr,obs, obs_tr, log, model)
}
 


#Creation of the result file into a CSV format--------------------------------------------------------------------------------------------
results = data.frame(results)
colnames(results)[1:4] <- c("label","AUC test","AUC train","training time")
colnames(results)[5:6004] <- c(1:6000)
write.csv(results,"results.csv")

