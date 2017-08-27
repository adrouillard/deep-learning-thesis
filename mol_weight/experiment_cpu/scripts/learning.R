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
# Connect to database and import data---------------------------------------------------------------------------------------------------
mydb = dbConnect(MySQL(), user='root', password='Macrobio21!', dbname='mass_spectra', host = '127.0.0.1')
#peak = dbReadTable(mydb, "peaks")
mol = dbReadTable(mydb, "molecule")
tacc = read.csv("../train_access.csv")
test_acc = tacc[,2]
#import parameters----------------------------------------------------------------------------------------------------------------------
para = read.csv("../parameters.csv", header = TRUE)
val1 = na.omit(as.matrix(expand.grid(s=para$seed, l1 = para$neurons)))
val2 = na.omit(as.matrix(expand.grid(s=para$seed, l1 = para$neurons, l2 = para$neurons)))
val3 = na.omit(as.matrix(expand.grid(s=para$seed, l1 = para$neurons, l2 = para$neurons, l3 = para$neurons)))

#*********Selection to creat a reduce the number of neurons at each layer
val2 = subset(val2, val2[,2] > val2[,3])
val3 = subset(val3, val3[,2] > val3[,3])
val3 =subset(val3, val3[,3] > val3[,4])
value = list(val1=val1,val2=val2,val3=val3)

results = matrix(data = NA, nrow =  dim(val1)[1]+dim(val2)[1]+dim(val3)[1] , ncol =  1000 + 4)

#import a selected numkber of peaks ----------------------------------------------------------------------------------------------------
sel = dbSendQuery(mydb, "select * from peaks where mz < 300")
res_sel= dbFetch(sel, n = -1)

#create matricies for analysis----------------------------------------------------------------------------------------------------------
data = read.csv("martix_bin1.csv")
data = data[,-1]
data = subset.data.frame(data, data$X2 < 1000)
data_test = data[test_acc,]
data_train = data[!(data$X1 %in% test_acc),]

#********Creation of input matrices*************
Mtrain = data.matrix(data_train[,-1])
Mtest = data.matrix(data_test[,-1])
train.x <- Mtrain[,-1]
train.y <- Mtrain[,1]
test.x <- Mtest[,-1]

#Code for GPU and CPU------------------------------------------------------------------------------------------------------------------- 
n.gpu <- 1
device.cpu <- mx.cpu()
device.gpu <- lapply(0:(n.gpu-1), function(i) {
  mx.gpu(i)
})

# Main loop testing --------------------------------------------------------------------------------------------------------------------
for (lay in 1:3){
  for (repet in 1:dim(value[[lay]])[1]){
    # Deep NN-------------------------------------------------------------------------------------------------------------------------------
    if (lay == 1) {
      l1 = value[[lay]][repet,][2]
      l2 = 1
      l3 = 1
      x = "act1"
    } else if (lay == 2) {
      l1 = value[[lay]][repet,][2]
      l2 = value[[lay]][repet,][3]
      l3 = 1
      x = "act2"
    } else {
      l1 = value[[lay]][repet,][2]
      l2 = value[[lay]][repet,][3]
      l3 = value[[lay]][repet,][4]
      x = "act3"
    }
    #********Creation of the architechture************
    data <- mx.symbol.Variable("data")
    fc1 <- mx.symbol.FullyConnected(data, name="fc1", num_hidden=l1)
    act1 <- mx.symbol.Activation(fc1, name="sig1", act_type="sigmoid")
    fc2 <- mx.symbol.FullyConnected(act1, name="fc2", num_hidden=l2)
    act2 <- mx.symbol.Activation(fc2, name="sig2", act_type="sigmoid")
    fc3 <-  mx.symbol.FullyConnected(act2, name="fc3", num_hidden=l3)
    act3 <- mx.symbol.Activation(fc3, name="sig3", act_type="sigmoid")
    fc4 <-  mx.symbol.FullyConnected(get(x), name="fc4", num_hidden=1)
    out <- mx.symbol.LinearRegressionOutput(fc4, name="sm")
    
    
    #Learning function----------------------------------------------------------------------------------------------------------------------
    mx.set.seed(value[[lay]][repet,][1])
    beg <- proc.time()
    log <- mx.metric.logger$new()
    model <- mx.model.FeedForward.create(out,
                                         X=train.x,
                                         y=train.y,
                                         ctx=device.cpu,
                                         num.round=1000,
                                         array.batch.size=200,
                                         learning.rate=5.5e-6,
                                         momentum=0.9, 
                                         eval.metric=mx.metric.rmse,
                                         batch.end.callback = mx.callback.log.train.metric(100, log),
                                         initializer=mx.init.uniform(0.07)
    )
    t <- (proc.time() - beg)
    
    #testing function----------------------------------------------------------------------------------------------------------------------
    beg2 <- proc.time()
    preds <- predict(model,test.x)
    print(proc.time() - beg2)
    
    beg2 <- proc.time()
    preds_tr <- predict(model,train.x)
    print(proc.time() - beg2)
    
    #accuracy testing-----------------------------------------------------------------------------------------------------------------------
    Q2 <- 1-( sum((Mtest[,1]-preds)^2) / sum((Mtest[,1]-mean(Mtest[,1]))^2) )
    Q2
    
    Q2_tr <- 1-( sum((train.y-preds_tr)^2) / sum((train.y-mean(train.y))^2) )
    Q2_tr
    
    #Labels---------------------------------------------------------------------------------------------------------------------------------
    bin = "bin1"
    lab = paste(bin,"_",value[[lay]][repet,][1],"_",l1,"-",l2,"-",l3, sep = '')
    Glab = paste(lab,"_results.png", sep = "")
    Modellab = paste("model/",lab,"_model", sep = "")
    
    #Save model-----------------------------------------------------------------------------------------------------------------------------
    mx.model.save(model, Modellab,1)
    res = c(Q2, Q2_tr,t[1], log$train)
    rep_res = data.frame(round(res[1:3], 4))
    rownames(rep_res)[1:3] <- c("Q2 test", "Q2 train","training time")
    colnames(rep_res) <- paste("results",lab, sep = ' ')
    
    res = c(lab,res)
    results[repet,] = res
    #Graph ---------------------------------------------------------------------------------------------------------------------------------
    png(Glab, width = 1600, height = 1200)
    par(mfrow=c(2,2))
    
    plot(log$train, type = 'l', main = "Learning curve", ylab = "Root mean square error")
    
    textplot(rep_res)
    
    obs = array(Mtest[,1])
    plot(array(preds[1,])~obs, main = "Prediction results Test set", ylab = "Prediction", xlab = "Observation" )
    abline(a = 0, b = 1, col = "red")
    
    obs_tr = array(train.y)
    plot(array(preds_tr[1,])~obs_tr, main = "Prediction results Train set", ylab = "Prediction", xlab = "Observation" )
    abline(a = 0, b = 1, col = "red")
    dev.off()
  }
}

#Creation of the result file into a CSV format--------------------------------------------------------------------------------------------
results = data.frame(results)
colnames(results)[1:4] <- c("label","Q2 test","Q2 train","training time")
colnames(results)[5:1004] <- c(1:1000)
write.csv(results,"results.csv")

