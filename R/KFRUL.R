library(kernlab)

#' Create MLP and SVM-RBF Models from the inputted data.
#' @title MLP, SVM-RBF models (model.ens.h2o.rbf)
#' @param data A matrix or vector of the independent data for the models.
#' @param rul A vector of the inputted dependent data for the models.
#' @return A list containing a MLP model created using h2o and a SVM-RBF model created using kernlab
#' @export
model.ens.h2o.rbf <- function(data,rul){
  model.mlp = h2o.deeplearning(x = 1:ncol(data), y = (ncol(data)+1), training_frame = as.h2o(cbind(data,rul)),hidden = c(200,200,200))
  model.rbf = ksvm(x = data, y  = rul, kernal = "rbfdot")
  return(list(mlp = model.mlp, rbf = model.rbf))
}

#' Create a model for combining the results of the MLP and SVM-RBF models. This will be called the Ensemble Model.
#' @title Ensemble Model (ensemble.ens.h2o)
#' @param data A matrix or vector of the independent data for the model.
#' @param rul A vector of the inputted dependent data for the models.
#' @param model The MLP and SVM-RBF models from the 'model.ens.h2o.rbf' function.
#' @return An ensemble model created using MLP from h2o
#' @export
ensemble.ens.h2o <- function(data, rul, model){
  pred.mlp = as.data.frame(h2o.predict(model$mlp,as.h2o(data)))$predict
  pred.rbf = predict(model$rbf, data)
  pred = cbind(pred.mlp,pred.rbf)
  ensemble.mlp = h2o.deeplearning(x = 1:(ncol(pred)), y = (ncol(pred)+1), training_frame = as.h2o(cbind(pred,rul)),hidden = c(200,200,200))
  return(ensemble.mlp)
}

#' Validate created model with a subset of data from the training set
#' @title Validate Ensemble Model (validate.ens.h2o)
#' @param data The validation data as a vector
#' @param model The models created by the 'model.ens.h2o.rbf' function
#' @param ensembleModel The model created by the 'ensemble.ens.h2o' function
#' @return The validated prediction as a vector
validate.ens.h2o <- function(data, model, ensembleModel){
  pred.mlp = as.data.frame(h2o.predict(model$mlp,as.h2o(data)))$predict
  pred.rbf = predict(model$rbf, data)
  pred = cbind(pred.mlp,pred.rbf)
  validatePred = as.data.frame(h2o.predict(ensembleModel,as.h2o(pred)))$predict
  return(validatePred)
}

#' Create Ensemble Model for RUL prediction
#' @title Create Ensemble Model (createModel.ens.h2o)
#' @param data The training data for the model
#' @param trainingSplit The percentage to split the data by for model training, ensemble training and validation
#' @return An ensemble model
createModel.ens.h2o <- function(data, trainingSplit = c(0.75,0.15,0.1)){
  len = nrow(data)
  # len = max(data[,1])
  inputCol = ncol(data)

  # rul = max(data[,1]) - data[,1]
  # rul = nrow(data):1

  rul = foreach(i = 1:max(data[,1]), .combine = 'c', .inorder = TRUE) %do%
    rev(data[which(data[,1]==i),2])

  data = data[,-1]

  data = scale(data)
  scaleMean = attr(data,"scaled:center")
  scaleStd = attr(data,"scaled:scale")
  # dataS = data

  unitId = 1:len

  modelIndex = sample(1:len, as.integer(len*trainingSplit[1]), FALSE)
  # unitId = unitId[-modelIndex]
  # modelIndex = as.numeric(unlist(sapply(modelIndex, function(x){which(data[,1]==x)})))
  dataModel = data[modelIndex,]
  rulModel = rul[modelIndex]
  data = data[-modelIndex,]
  rul = rul[-modelIndex]
  # dataModel = dataModel[,-1]

  ensembleIndex = sample(1:nrow(data), as.integer(len*trainingSplit[2]), FALSE)
  # ensembleIndex = as.numeric(unlist(sapply(ensembleIndex, function(x){which(data[,1]==unitId[x])})))
  dataEnsemble = data[ensembleIndex,]
  rulEnsemble = rul[ensembleIndex]
  data = data[-ensembleIndex,]
  rul = rul[-ensembleIndex]
  # dataEnsemble = dataEnsemble[,-1]

  removeCol = numeric()
  for(i in 1:(ncol(data))){
    if(length(which(is.na(dataModel[,i])))>0 || all(dataModel[,i]==dataModel[1,i]))
      removeCol = c(removeCol,i)
    else if(length(which(is.na(dataEnsemble[,i])))>0 || all(dataEnsemble[,i]==dataEnsemble[1,i]))
      removeCol = c(removeCol,i)
    else if(length(which(is.na(data[,i])))>0 || all(data[,i]==data[1,i]))
      removeCol = c(removeCol,i)
  }
  if(length(removeCol)!=0){
    dataModel = dataModel[,-removeCol]
    dataEnsemble = dataEnsemble[,-removeCol]
    data = data[,-removeCol]
    scaleMean = scaleMean[-removeCol]
    scaleStd = scaleStd[-removeCol]
  }

  models = model.ens.h2o.rbf(dataModel, rulModel)

  ensembleModel = ensemble.ens.h2o(dataEnsemble, rulEnsemble, models)

  validPred = validate.ens.h2o(data, models, ensembleModel)

  return(list(models = models,ensembleModel = ensembleModel, inputCol = inputCol, removeCol = removeCol, scaleMean = scaleMean, scaleStd = scaleStd, validate = validPred, validateRUL = rul))
}

#' Predict Using Ensemble Model
#' @title Predicting with Ensemble Model (predict.ens.h2o)
#' @param data The data that will have its RUL predicted
#' @param model An Emsemble model outputed by the 'createModel.ens.h2o' function
#' @return The predicted RUL of the inputted data
predict.ens.h2o <- function(data, model){
  if(ncol(data)!=model$inputCol){
    warning("Not the same amount of columns as training set")
    return(0)
  }
  data = data[,-1]
  # data = scale(data)
  if(length(model$removeCol!=0))
    data = data[,-model$removeCol]
  for(i in 1:ncol(data))
    data[,i] = (data[,i]-model$scaleMean[i])/model$scaleStd[i]

  pred.mlp = as.data.frame(h2o.predict(model$models$mlp,as.h2o(data)))
  pred.mlp = pred.mlp$predict
  pred.rbf = predict(model$models$rbf, data)

  pred.ens = h2o.predict(model$ensembleModel, as.h2o(cbind(pred.mlp,pred.rbf)))
  return(as.data.frame(pred.ens)$predict)
}

phmScoring <- function(predRUL, actualRUL){
  score = foreach(i = 1:length(predRUL),.combine = 'c', .inorder = TRUE) %do% {
    diff = predRUL[i]-actualRUL[i]
    if(diff<0)
      exp(-1*diff/10)-1
    else
      exp(diff/13)-1
  }
  return(sum(score))
}
#' library(ipred)
#'
#' #' @export
#' createModel.h2o.kf <- function(data, groupSize, nbaggs){
#'   boundSeq = as.integer(seq(from = 0, to = nrow(data), by = groupSize))
#'   rul = boundSeq[length(boundSeq)]:1
#'   registerDoMC()
#'   dataReps = foreach(i = 1:(length(boundSeq)-1), .combine = 'cbind', .inorder = TRUE) %dopar% {
#'     bag = ipredbagg(y = rul[(boundSeq[i]+1):(boundSeq[i+1])],X = data[(boundSeq[i]+1):(boundSeq[i+1]),],nbagg = nbaggs)
#'     reps = foreach(ii = 1:nbaggs, .combine = 'cbind') %do% {
#'       bag$mtrees[[ii]]$bindx
#'       print(paste("ii = ",ii,sep = ""))
#'     }
#'     print(paste("----------------------i = ",i,sep = ""))
#'     reps
#'   }
#'   registerDoMC()
#'   models = foreach(i = 1:ncol(dataReps)) %dopar% {
#'     group = (i/nbaggs)
#'     df = data.frame(unname(dataReps[,i]),rul = rul[((groupSize*group)+1):(groupSize*(group+1))])
#'     h2o.deeplearning(x = 1, y = 2, training_frame = as.h2o(df))
#'     print(paste("model = ",i,sep = ""))
#'   }
#'   return(models)
#' }
#'
#' predict.h2o.kf <- function(data, model){
#'   bag = ipredbagg()
#' }
#'
# model <- SSModel(as.numeric(pred.rbf)[1:10000] ~ SSMtrend(1, Q = list(matrix(NA))), H = matrix(NA))
# model <- fitSSM(model, c(log(var(as.numeric(pred.rbf)[1:10000])), log(var(as.numeric(pred.rbf)[1:10000]))),method = "BFGS")$model
# out <- KFS(model, filtering = "state", smoothing = "state")
#
#
# model <- SSModel(as.numeric(pred.rbf) ~ SSMtrend(1, Q = list(matrix(NA))),H = matrix(NA))
# out_NileNA <- KFS(model_NileNA, "mean", "mean")
#
# Expected = seq(from = (as.integer(length(data)/exp2)), to = (as.integer(length(data)/(exp2*length(RUL)))), length.out = length(RUL))
# RUL = cbind(RUL,Expected)
# Time = seq

# FD1 = read.table("~/Downloads/CMAPSSData/train_FD001.txt")
# FD2 = read.table("~/Downloads/CMAPSSData/train_FD002.txt")
# FD3 = read.table("~/Downloads/CMAPSSData/train_FD003.txt")
# FD4 = read.table("~/Downloads/CMAPSSData/train_FD004.txt")
# test1 = read.table("~/Downloads/CMAPSSData/test_FD001.txt")
# test2 = read.table("~/Downloads/CMAPSSData/test_FD002.txt")
# test3 = read.table("~/Downloads/CMAPSSData/test_FD003.txt")
# test4 = read.table("~/Downloads/CMAPSSData/test_FD004.txt")
# rul1 = read.table("~/Downloads/CMAPSSData/RUL_FD001.txt")
# rul2 = read.table("~/Downloads/CMAPSSData/RUL_FD002.txt")
# rul3 = read.table("~/Downloads/CMAPSSData/RUL_FD003.txt")
# rul4 = read.table("~/Downloads/CMAPSSData/RUL_FD004.txt")
# ensModel1 = createModel.ens.h2o(FD1)
# ensModel2 = createModel.ens.h2o(FD2)
# ensModel3 = createModel.ens.h2o(FD3)
# ensModel4 = createModel.ens.h2o(FD4)
# predTr1_Te1 = predict.ens.h2o(test1,ensModel1)
# predTr2_Te1 = predict.ens.h2o(test1,ensModel2)
# predTr3_Te1 = predict.ens.h2o(test1,ensModel3)
# predTr4_Te1 = predict.ens.h2o(test1,ensModel4)
#
# predTr1_Te2 = predict.ens.h2o(test2,ensModel1)
# predTr2_Te2 = predict.ens.h2o(test2,ensModel2)
# predTr3_Te2 = predict.ens.h2o(test2,ensModel3)
# predTr4_Te2 = predict.ens.h2o(test2,ensModel4)
#
# predTr1_Te3 = predict.ens.h2o(test3,ensModel1)
# predTr2_Te3 = predict.ens.h2o(test3,ensModel2)
# predTr3_Te3 = predict.ens.h2o(test3,ensModel3)
# predTr4_Te3 = predict.ens.h2o(test3,ensModel4)
#
# predTr1_Te4 = predict.ens.h2o(test4,ensModel1)
# predTr2_Te4 = predict.ens.h2o(test4,ensModel2)
# predTr3_Te4 = predict.ens.h2o(test4,ensModel3)
# predTr4_Te4 = predict.ens.h2o(test4,ensModel4)

# predTr1_Te1 = predict.ens.h2o(test1,ensModel1)
# predTr1_Te2 = predict.ens.h2o(test2,ensModel1)
# predTr1_Te3 = predict.ens.h2o(test3,ensModel1)
# predTr1_Te4 = predict.ens.h2o(test4,ensModel1)


# i = sapply(1:max(test2[,1]),function(x){
#   w = which(test2[,1]==x)
#   max(w)
# })
#
# fullRUL <- function(testdata,ruldata){
#   newRUL = numeric()
#   for(i in 1:nrow(ruldata)){
#     tempRUL = ruldata[i,] + length(which(testdata[,1]==i)) - 1
#     newRUL = c(newRUL,tempRUL:ruldata[i,])
#   }
#   return(newRUL)
# }

