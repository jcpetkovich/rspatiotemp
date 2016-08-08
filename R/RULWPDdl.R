
#' Create deep learning model for nodal energies or wpd using h2o.
#' @title Create Model WPD energy RUL (createModel.h2o.wpd.rul)
#' @param data The time series data.
#' @param exp2 Two to the power of exp2 will be the size of each group fed into the WPD of each group fed into the WPD.
#' @param levels The level of the WPD to compute the energy for.
#' @param degradStartPercent The percentage at which we assume the system begins degradation.
#' @param scale TRUE or FALSE whether or not to scale the output
#' @return A deeplearning model from h2o of the wpd nodal energies
#' @export
createModel.h2o.wpd.rul <- function(data,exp2,levels,degradStartPercent,scale){
  energy = getEnergies.wpd(data,exp2,levels)
  time = 1:nrow(energy)
  registerDoMC()
  energy = foreach(i = 1:ncol(energy),.combine = 'cbind',.inorder = TRUE) %dopar% {
    predict(loess(energy[,i]~time))
  }
  # for(i in 1:ncol(energy))
  #   energy[,i] = predict(loess(energy[,i]~time))

  rul = c(rep(nrow(energy),as.integer(nrow(energy)*degradStartPercent)),seq(from=(nrow(energy)), to = 1, length.out = (nrow(energy))-as.integer(nrow(energy)*degradStartPercent)))
  energy.df = data.frame(energy,rul = rul)
  energy.hex = as.h2o(energy.df)
  energy.dl = h2o.deeplearning(x = 1:(2^levels),y = (ncol(energy.df)),training_frame = energy.hex)
  return(list(model = energy.dl,exp2 = exp2, levels = levels,scale = scale))
}

#' Predict the RUL by using the nodal energies of the WPD
#' @title WPD RUL (predict.h2o.wpd.rul)
#' @param data The time series data
#' @param model The deep learning model created by the 'createModel.h2o.wpd.rul' function
#' @return A vector of the predicted RUL
#' @export
predict.h2o.wpd.rul <- function(data,model){
  levels = model$levels
  exp2 = model$exp2
  scale = model$scale
  model = model$model
  energy = getEnergies.wpd(data,exp2,levels)
  time = 1:nrow(energy)
  registerDoMC()
  energy = foreach(i = 1:ncol(energy),.combine = 'cbind',.inorder = TRUE) %dopar% {
    predict(loess(energy[,i]~time))
  }
  # for(i in 1:ncol(energy))
  #   energy[,i] = predict(loess(energy[,i]~time))
  rul = rep(0,nrow(energy))
  energy.df = data.frame(energy,rul = rul)
  energy.hex = as.h2o(energy.df)
  pred = h2o.predict(model,energy.hex)
  pred = as.data.frame(pred)
  return(pred$predict)
}
