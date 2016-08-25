#' Create deep learning model for nodal energies of wpd using h2o.
#' @title Create Model WPD energy RUL (createModel.h2o.wpd.rul)
#' @param data The time series data.
#' @param exp2 Two to the power of exp2 will be the size of each group fed into the WPD of each group fed into the WPD.
#' @param levels The level of the WPD to compute the energy for.
#' @param degradStartPercent The percentage at which we assume the system begins degradation.
#' @param nPast Number of past energy states to look at when predicting RUL. Default is 0
#' @param cores The number of cores to use when creating the model. Default is NULL, which will use half of the detected number of cores.
#' @return A deeplearning model from h2o of the wpd nodal energies
#' @export
createModel.h2o.wpd.rul <- function(data,exp2,levels,degradStartPercent,nPast = 0,cores=NULL){
  energy = getEnergies.wpd(data,exp2,levels)
  time = 1:nrow(energy)
  registerDoMC(cores = cores)
  energy = foreach(i = 1:ncol(energy),.combine = 'cbind',.inorder = TRUE) %dopar% {
    predict(loess(energy[,i]~time))
  }
  # for(i in 1:ncol(energy))
  #   energy[,i] = predict(loess(energy[,i]~time))
  if(nPast != 0){
    rowNum = nrow(energy)-nPast
    registerDoMC(cores = cores)
    energy = foreach(i = 0:rowNum, .combine = 'rbind', .inorder = TRUE) %dopar% {
      lb = i+1
      ub = nPast+i
      as.vector(energy[lb:ub,])
    }
    energy = unname(energy)
  }
  rul = c(rep(nrow(energy),as.integer(nrow(energy)*degradStartPercent)),seq(from=(nrow(energy)), to = 1, length.out = (nrow(energy))-as.integer(nrow(energy)*degradStartPercent)))
  energy.df = data.frame(energy,rul = rul)
  energy.hex = as.h2o(energy.df)
  energy.dl = h2o.deeplearning(x = 1:(2^levels),y = (ncol(energy.df)),training_frame = energy.hex)
  return(list(model = energy.dl,exp2 = exp2, levels = levels,nPast = nPast))
}

#' Compute the energy and RUL for the inputted data and return it as a data frame.
#' @title Compute Energy and RUL (getEnergy.df.rul)
#' @param data The time series data.
#' @param exp2 Two to the power of exp2 will be the size of each group fed into the WPD of each group fed into the WPD.
#' @param levels The level of the WPD to compute the energy for.
#' @param degradStartPercent The percentage at which we assume the system begins degradation.
#' @param nPast Number of past energy states to look at when predicting RUL. Default is 0
#' @param cores The number of cores to use when creating the model. Default is NULL, which will use half of the detected number of cores.
#' @return A data frame containing the nodal energies and RUL, used for h2o.deeplearning training.
#' @export
getEnergy.df.rul <- function(data,exp2,levels,degradStartPercent,nPast = 0,cores=NULL){
  energy = getEnergies.wpd(data,exp2,levels)
  time = 1:nrow(energy)
  registerDoMC(cores = cores)
  energy = foreach(i = 1:ncol(energy),.combine = 'cbind',.inorder = TRUE) %dopar% {
    predict(loess(energy[,i]~time))
  }
  if(nPast != 0){
    rowNum = nrow(energy)-nPast
    registerDoMC(cores = cores)
    energy = foreach(i = 0:rowNum, .combine = 'rbind', .inorder = TRUE) %dopar% {
      lb = i+1
      ub = nPast+i
      as.vector(energy[lb:ub,])
    }
    energy = unname(energy)
  }
  rul = c(rep(nrow(energy),as.integer(nrow(energy)*degradStartPercent)),seq(from=(nrow(energy)), to = 1, length.out = (nrow(energy))-as.integer(nrow(energy)*degradStartPercent)))
  energy.df = data.frame(energy,rul = rul)
  return(energy.df)
}

#' Create deep learning model for nodal energies of wpd using h2o using several sets of data
#' @title Create Model WPD energy RUL (createModel.h2o.wpd.rul)
#' @param data The A list of time series data.
#' @param exp2 Two to the power of exp2 will be the size of each group fed into the WPD of each group fed into the WPD.
#' @param levels The level of the WPD to compute the energy for.
#' @param degradStartPercent The percentage at which we assume the system begins degradation.
#' @param nPast Number of past energy states to look at when predicting RUL. Default is 0
#' @param cores The number of cores to use when creating the model. Default is NULL, which will use half of the detected number of cores.
#' @return A deeplearning model from h2o of the wpd nodal energies
#' @export
createModel.h2o.wpd.group.rul <- function(data,exp2,levels,degradStartPercent,hidden,nPast = 0,cores=NULL){
  energy.df = foreach(i = 1:length(data), .combine = 'rbind') %do% {
    getEnergy.df.rul(data[[i]],exp2,levels,degradStartPercent,nPast,cores)
  }
  energy.hex = as.h2o(energy.df)
  energy.dl = h2o.deeplearning(x = 1:(2^levels),y = (ncol(energy.df)),training_frame = energy.hex,hidden = hidden)
  return(list(model = energy.dl,exp2 = exp2, levels = levels,nPast = nPast))
}

#' Predict the RUL by using the nodal energies of the WPD
#' @title WPD RUL (predict.h2o.wpd.rul)
#' @param data The time series data
#' @param model The deep learning model created by the 'createModel.h2o.wpd.rul' function
#' @param cores The number of cores to use when predicting the RUL Default is NULL, which will use half of the detected number of cores.
#' @return A vector of the predicted RUL
#' @export
predict.h2o.wpd.rul <- function(data,model,cores = NULL){
  levels = model$levels
  exp2 = model$exp2
  nPast = model$nPast
  model = model$model
  energy = getEnergies.wpd(data,exp2,levels)
  time = 1:nrow(energy)
  registerDoMC(cores = cores)
  energy = foreach(i = 1:ncol(energy),.combine = 'cbind',.inorder = TRUE) %dopar% {
    predict(loess(energy[,i]~time))
  }
  if(nPast != 0){
    rowNum = nrow(energy)-nPast
    registerDoMC(cores = cores)
    energy = foreach(i = 0:rowNum, .combine = 'rbind', .inorder = TRUE) %dopar% {
      lb = i+1
      ub = nPast+i
      as.vector(energy[lb:ub,])
    }
    energy = unname(energy)
  }
  # for(i in 1:ncol(energy))
  #   energy[,i] = predict(loess(energy[,i]~time))
  #rul = rep(0,nrow(energy))
  energy.df = data.frame(energy)
  energy.hex = as.h2o(energy.df)
  pred = h2o.predict(model,energy.hex)
  pred = as.data.frame(pred)
  return(pred$predict)
}
# plot.ts(predict.h2o.wpd.rul(h1_1,model1_1))
# plot.ts(predict.h2o.wpd.rul(h1_2,model1_2))
# plot.ts(predict.h2o.wpd.rul(h1_3,model1_3))
# plot.ts(predict.h2o.wpd.rul(h1_4,model1_4))
# plot.ts(predict.h2o.wpd.rul(h1_5,model1_5))
# plot.ts(predict.h2o.wpd.rul(h1_6,model1_6))
# plot.ts(predict.h2o.wpd.rul(h1_7,model1_7)
#
# data = list(h1_1,h1_2,h1_3,h1_4,h1_5,h1_6,h1_7)
# model1_1 = createModel.h2o.wpd.group.rul(data[-1],8,3,0)
# model1_2 = createModel.h2o.wpd.group.rul(data[-2],8,3,0)
# model1_3 = createModel.h2o.wpd.group.rul(data[-3],8,3,0)
# model1_4 = createModel.h2o.wpd.group.rul(data[-4],8,3,0)
# model1_5 = createModel.h2o.wpd.group.rul(data[-5],8,3,0)
# model1_6 = createModel.h2o.wpd.group.rul(data[-6],8,3,0)
# model1_7 = createModel.h2o.wpd.group.rul(data[-7],8,3,0)
#
# plot.ts(predict.h2o.wpd.rul(h1_1,model1_1))
# plot.ts(predict.h2o.wpd.rul(h1_2,model1_2))
# plot.ts(predict.h2o.wpd.rul(h1_3,model1_3))
# plot.ts(predict.h2o.wpd.rul(h1_4,model1_4))
# plot.ts(predict.h2o.wpd.rul(h1_5,model1_5))
# plot.ts(predict.h2o.wpd.rul(h1_6,model1_6))
# plot.ts(predict.h2o.wpd.rul(h1_7,model1_7))
#
#
# model1_1 = h2o.deeplearning(x=1:8,y=9,training_frame = as.h2o(all.df1_1))
# model1_2 = h2o.deeplearning(x=1:8,y=9,training_frame = as.h2o(all.df1_2))
# model1_3 = h2o.deeplearning(x=1:8,y=9,training_frame = as.h2o(all.df1_3))
# model1_4 = h2o.deeplearning(x=1:8,y=9,training_frame = as.h2o(all.df1_4))
# model1_5 = h2o.deeplearning(x=1:8,y=9,training_frame = as.h2o(all.df1_5))
# model1_6 = h2o.deeplearning(x=1:8,y=9,training_frame = as.h2o(all.df1_6))
# model1_7 = h2o.deeplearning(x=1:8,y=9,training_frame = as.h2o(all.df1_7))
#
# model1_1 = list(model = model1_1,exp2 = 8,levels = 3,nPast = 0)
# model1_2 = list(model = model1_2,exp2 = 8,levels = 3,nPast = 0)
# model1_3 = list(model = model1_3,exp2 = 8,levels = 3,nPast = 0)
# model1_4 = list(model = model1_4,exp2 = 8,levels = 3,nPast = 0)
# model1_5 = list(model = model1_5,exp2 = 8,levels = 3,nPast = 0)
# model1_6 = list(model = model1_6,exp2 = 8,levels = 3,nPast = 0)
# model1_7 = list(model = model1_7,exp2 = 8,levels = 3,nPast = 0)

# modelH1_3p1 = createModel.h2o.wpd.group.rul(data = list(h1_1,h1_2,h1_4,h1_5,h1_6,h1_7),exp2 = 8,levels = 3,degradStartPercent = 0,hidden = c(200,200,200,200),nPast = 2,cores = 4)
# plot.ts(predict.h2o.wpd.rul(h1_3,modelH1_3p1))
# modelH1_3p2 = createModel.h2o.wpd.group.rul(data = list(h1_1,h1_2,h1_4,h1_5,h1_6,h1_7),exp2 = 8,levels = 3,degradStartPercent = 0,hidden = c(200,200,200,200),nPast = 2,cores = 4)
# plot.ts(predict.h2o.wpd.rul(h1_3,modelH1_3p2))
# modelH1_3p3 = createModel.h2o.wpd.group.rul(data = list(h1_1,h1_2,h1_4,h1_5,h1_6,h1_7),exp2 = 8,levels = 3,degradStartPercent = 0,hidden = c(200,200,200,200),nPast = 2,cores = 4)
# plot.ts(predict.h2o.wpd.rul(h1_3,modelH1_3p3))


