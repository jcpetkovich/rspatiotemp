#' Create deep learning model for nodal energies or wpd using h2o.
#' @title Create Model WPD energy RUL (createModel.h2o.wpd.rul)
#' @param data The time series data.
#' @param exp2 Two to the power of exp2 will be the size of each group fed into the WPD of each group fed into the WPD.
#' @param levels The level of the WPD to compute the energy for.
#' @param degradStartPercent The percentage at which we assume the system begins degradation.
#' @param nPast Number of past energy states to look at when predicting RUL. Default is 0
#' @param cores The number of cores to use when creating the model. Default is NULL, which will use half of the detect number of cores.
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
  return(energy.df)
}
#   energy.hex = as.h2o(energy.df)
#   energy.dl = h2o.deeplearning(x = 1:(2^levels),y = (ncol(energy.df)),training_frame = energy.hex)
#   return(list(model = energy.dl,exp2 = exp2, levels = levels,nPast = nPast))
# }

#' Predict the RUL by using the nodal energies of the WPD
#' @title WPD RUL (predict.h2o.wpd.rul)
#' @param data The time series data
#' @param model The deep learning model created by the 'createModel.h2o.wpd.rul' function
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
plot.ts(predict.h2o.wpd.rul(h1_1,modelH1_3_p0))
plot.ts(predict.h2o.wpd.rul(h1_2,modelH1_3_p0))
plot.ts(predict.h2o.wpd.rul(h1_3,modelH1_3_p0))
plot.ts(predict.h2o.wpd.rul(h1_4,modelH1_3_p0))
plot.ts(predict.h2o.wpd.rul(h1_5,modelH1_3_p0))
plot.ts(predict.h2o.wpd.rul(h1_6,modelH1_3_p0))
plot.ts(predict.h2o.wpd.rul(h1_7,modelH1_3_p0))

df1 = createModel.h2o.wpd.rul(h1_1,8,3,0)
df2 = createModel.h2o.wpd.rul(h1_2,8,3,0)
df3 = createModel.h2o.wpd.rul(h1_3,8,3,0)
df4 = createModel.h2o.wpd.rul(h1_4,8,3,0)
df5 = createModel.h2o.wpd.rul(h1_5,8,3,0)
df6 = createModel.h2o.wpd.rul(h1_6,8,3,0)
df7 = createModel.h2o.wpd.rul(h1_7,8,3,0)
