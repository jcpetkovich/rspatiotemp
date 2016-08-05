
#' @export
createModel.h2o.wpd.rul <- function(data,exp2,levels,degradStartPercent){
  energy = getEnergies.wpd(data,exp2,levels)
  rul = c(rep(nrow(energy),as.integer(nrow(energy)*degradStartPercent)),seq(from=(nrow(energy)), to = 1, length.out = (nrow(energy))-as.integer(nrow(energy)*degradStartPercent)))
  energy.df = data.frame(energy,rul = rul)
  energy.hex = as.h2o(energy.df)
  energy.dl = h2o.deeplearning(x = 1:8,y = 9,training_frame = energy.hex)
  return(list(model = energy.dl,exp2 = exp2, levels = levels))
}

#' @export
predict.h2o.wpd.rul <- function(data,model){
  levels = model$levels
  exp2 = model$exp2
  model = model$model
  energy = getEnergies.wpd(data,exp2,levels)
  rul = rep(0,nrow(energy))
  energy.df = data.frame(energy,rul = rul)
  energy.hex = as.h2o(energy.df)
  pred = h2o.predict(model,energy.hex)
  pred = as.data.frame(pred)
  return(pred$predict)
}
