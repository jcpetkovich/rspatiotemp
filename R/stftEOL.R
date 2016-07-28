library(h2o)
h2o.init()
#' @export
get.fftloess <- function(data, groupSize){
  boundSeq = seq(from=1,to = (length(data)+1), by = groupSize)
  len = length(boundSeq) - 1
  time = 1:groupSize
  lowerBound = boundSeq[1]
  upperBound = boundSeq[2] - 1
  fft = abs(fft(data[lowerBound:upperBound]))
  lo = predict(loess(fft~time))

  for(i in 2:len){
    lowerBound = boundSeq[i]
    upperBound = boundSeq[i+1] - 1
    fft = abs(fft(data[lowerBound:upperBound]))
    lo1 = predict(loess(fft~time))
    lo = rbind(lo,lo1)
  }

  return(lo)
}

#' @export
get.fftloess.skip <- function(data, groupSize,smoothSampleNum){
  boundSeq = seq(from=1,to = (length(data)+1), by = groupSize)
  len = length(boundSeq) - 1
  boundSeq2 = seq(from = 1, to = length(data),length.out = smoothSampleNum)
  time = 1:groupSize
  lowerBound = boundSeq[1]
  upperBound = boundSeq[2] - 1
  fft = abs(fft(data[lowerBound:upperBound]))
  lo = predict(loess(fft~time))

  for(i in 2:len){
    lowerBound = boundSeq[i]
    upperBound = boundSeq[i+1] - 1
    fft = abs(fft(data[lowerBound:upperBound]))
    lo1 = predict(loess(fft~time))[boundSeq2]
    lo = rbind(lo,lo1)
  }
  return(lo)
}

#' @export
get.fftloess.d10 <- function(data, groupSize){
  boundSeq = seq(from=1,to = (length(data)+1), by = groupSize)
  len = length(boundSeq) - 1
  winSize = as.integer(groupSize/20)
  lb = as.integer(groupSize/2-winSize)
  ub = as.integer(groupSize/2+winSize)
  time = 1:groupSize
  lowerBound = boundSeq[1]
  upperBound = boundSeq[2] - 1
  fft = abs(fft(data[lowerBound:upperBound]))
  lo = predict(loess(fft~time))[lb:ub]

  for(i in 2:len){
    lowerBound = boundSeq[i]
    upperBound = boundSeq[i+1] - 1
    fft = abs(fft(data[lowerBound:upperBound]))
    lo1 = predict(loess(fft~time))[lb:ub]
    lo = rbind(lo,lo1)
  }
  return(lo)
}

#' @export
createModel.h2o.logic <- function(data,groupSize,degradStartPercent){
  fftloess = get.fftloess(data,groupSize)
  rul = c(rep(TRUE,as.integer(nrow(fftloess)*degradStartPercent)),rep(FALSE,(nrow(fftloess))-as.integer(nrow(fftloess)*degradStartPercent)))
  fftloess.df = data.frame(fftloess, rul = rul)
  fftloess.hex = as.h2o(fftloess.df)
  fftloess.dl = h2o.deeplearning(x = 1:(ncol(fftloess)),y = ncol(fftloess.df),training_frame = fftloess.hex)
  return(fftloess.dl)
}

#' @export
createModel.h2o.logic.save <- function(data,groupSize,degradStartPercent,bearingNum){
  fftloess = get.fftloess(data,groupSize)
  rul = c(rep(TRUE,as.integer(nrow(fftloess)*degradStartPercent)),rep(FALSE,(nrow(fftloess))-as.integer(nrow(fftloess)*degradStartPercent)))
  fftloess.df = data.frame(fftloess, rul = rul)
  save(fftloess.df,file = paste("fftloess",bearingNum,".Rd",sep = ""))
  fftloess.hex = as.h2o(fftloess.df)
  fftloess.dl = h2o.deeplearning(x = 1:(ncol(fftloess)),y = ncol(fftloess.df),training_frame = fftloess.hex)
  save(fftloess.dl, file = paste("model",bearingNum,".Rd",sep=""))
}

#' @export
predict.h2o.logic <- function(data,groupSize, model){
  testfft = get.fftloess(data,groupSize)
  testfft.df = data.frame(testfft)
  testfft.hex = as.h2o(testfft.df)
  predictions = h2o.predict(model,testfft.hex)
  predictions = as.data.frame(predictions)
  return(predictions$predict)
}

#' @export
status.h2o.logic <- function(dataH,dataV,modelH,modelV,groupSize,n){
  if(length(predH) > (n*groupSize)){
    up = length(predH)
    lo = length(predH)-n*groupSize+1
    dataH = dataH[lo:up]
    dataV = dataV[lo:up]
  }
  predH = predict.h2o.logic(dataH,groupSize,modelH)
  predV = predict.h2o.logic(dataV,groupSize,modelV)

  if(all(predH == TRUE)){
    if(all(predV == FALSE)){
      return(-1)
    }
    else{
      return(1)
    }
  }
  else if(all(predH == FALSE)){
    return(-1)
  }
  else{
    if(all(predV == TRUE)){
      return(1)
    }
    else if(all(predV == FALSE)){
      return(-1)
    }
    else{
      return(0)
    }
  }
}

#' @export
predict.status.h2o.logic <- function(dataH,dataV,groupSize, modelH,modelV,n){
  predictionsH = predict.h2o.logic(dataH,groupSize,modelH)
  predictionsV = predict.h2o.logic(dataV,groupSize,modelV)
  boundSeq = c(1:(n-1),seq(from=n,to = (length(data)+1), by = n))
  status = numeric()
  for(i in boundSeq){
    lo = numeric()
    if(i < n)
      lo = 1
    else
      lo = i - n + 1

    predH = predictionsH[lo:i]
    predV = predictionsV[lo:i]

    if(all(predH == TRUE)){
      if(all(predV == FALSE)){
        status = c(status,-1)
      }
      else{
        status = c(status,1)
      }
    }
    else if(all(predH == FALSE)){
      status = c(status,-1)
    }
    else{
      if(all(predV == TRUE)){
        status = c(status,1)
      }
      else if(all(predV == FALSE)){
        status = c(status,-1)
      }
      else{
        status = c(status,0)
      }
    }
  }
  return(status)
}

#' @export
selfRUL.h2o.fftloess <- function(data,groupSize,degradStartPercent){
  fftloess = get.fftloess(data,groupSize)
  rul = c(rep(nrow(fftloess),as.integer(nrow(fftloess)*degradStartPercent)),seq(from=(nrow(fftloess)), to = 1, length.out = (nrow(fftloess))-as.integer(nrow(fftloess)*degradStartPercent)))
  fftloess.df = data.frame(fftloess,rul = rul)
  fftloess.hex = as.h2o(fftloess.df)
  fftloess.dl = h2o.deeplearning(x = 1:(ncol(fftloess)),y = ncol(fftloess.df),training_frame = fftloess.hex)
  predictions <- h2o.predict(fftloess.dl, fftloess.hex)
  predictions = as.data.frame(predictions)
  plot.ts(predictions$predict)
}

#' @export
selfRUL.h2o.fftloess.logic <- function(data,groupSize,degradStartPercent){
  fftloess = get.fftloess(data,groupSize)
  rul = c(rep(TRUE,as.integer(nrow(fftloess)*degradStartPercent)),rep(FALSE,(nrow(fftloess))-as.integer(nrow(fftloess)*degradStartPercent)))
  fftloess.df = data.frame(fftloess,rul = rul)
  fftloess.hex = as.h2o(fftloess.df)
  fftloess.dl = h2o.deeplearning(x = 1:(ncol(fftloess)),y = ncol(fftloess.df),training_frame = fftloess.hex)
  predictions <- h2o.predict(fftloess.dl, fftloess.hex)
  predictions = as.data.frame(predictions)
  plot.ts(predictions$predict)
}

#' @export
selfRUL.h2o.fftloess.logic.skip <- function(data,groupSize,degradStartPercent,smoothSampleNum){
  fftloess = get.fftloess.skip(data,groupSize,smoothSampleNum)
  #rul = c(rep(nrow(fftloess),as.integer(nrow(fftloess)*degradStartPercent)),seq(from=(nrow(fftloess)), to = 1, length.out = (nrow(fftloess))-as.integer(nrow(fftloess)*degradStartPercent)))
  rul = c(rep(TRUE,as.integer(nrow(fftloess)*degradStartPercent)),rep(FALSE,(nrow(fftloess))-as.integer(nrow(fftloess)*degradStartPercent)))
  fftloess.df = data.frame(fftloess,rul = rul)
  fftloess.hex = as.h2o(fftloess.df)
  fftloess.dl = h2o.deeplearning(x = 1:(ncol(fftloess)-1),y = ncol(fftloess),training_frame = fftloess.hex)
  predictions <- h2o.predict(fftloess.dl, fftloess.hex)
  predictions = as.data.frame(predictions)
  plot.ts(predictions$predict)
}

#' @export
selfRUL.h2o.fftloess.d10 <- function(data,groupSize,degradStartPercent){
  fftloess = get.fftloess.d10(data,groupSize)
  rul = c(rep(nrow(fftloess),as.integer(nrow(fftloess)*degradStartPercent)),seq(from=(nrow(fftloess)), to = 1, length.out = (nrow(fftloess))-as.integer(nrow(fftloess)*degradStartPercent)))
  fftloess.df = data.frame(rul = rul, fftloess)
  fftloess.hex = as.h2o(fftloess.df)
  fftloess.dl = h2o.deeplearning(x = 2:(ncol(fftloess)),y = 1,training_frame = fftloess.hex)
  predictions <- h2o.predict(fftloess.dl, fftloess.hex)
  predictions = as.data.frame(predictions)
  plot.ts(predictions$predict)
}

#' @export
plot.fftloess <- function(data, groupSize){

  fft1 = abs(fft(data[1:groupSize]))
  fft3 = abs(fft(data[(length(data)-groupSize+1):(length(data))]))
  if(groupSize%%2 == 0){
    fft2 = abs(fft(data[(length(data)/2 - groupSize/2 + 1):(length(data)/2 + groupSize/2)]))
  }
  else{
    fft2 = abs(fft(data[(length(data)/2 - groupSize/2 + 0.5):(length(data)/2) + groupSize/2 - 0.5]))
  }
  time = 1:groupSize
  lo1 = predict(loess(fft1~time))
  lo2 = predict(loess(fft2~time))
  lo3 = predict(loess(fft3~time))
  plot.ts(fft1)
  lines(lo1,col= "red")
  lines(lo2,col= "blue")
  lines(lo3,col= "green")
}

#' @export
plot.fftloess.d10 <- function(data, groupSize){
  fft1 = abs(fft(data[1:groupSize]))
  fft3 = abs(fft(data[(length(data)-groupSize+1):(length(data))]))
  if(groupSize%%2 == 0){
    fft2 = abs(fft(data[(length(data)/2 - groupSize/2 + 1):(length(data)/2 + groupSize/2)]))
  }
  else{
    fft2 = abs(fft(data[(length(data)/2 - groupSize/2 + 0.5):(length(data)/2) + groupSize/2 - 0.5]))
  }
  winSize = as.integer(groupSize/20)
  lb = as.integer(groupSize/2-winSize)
  ub = as.integer(groupSize/2+winSize)
  time = 1:groupSize
  lo1 = predict(loess(fft1~time))[lb:ub]
  lo2 = predict(loess(fft2~time))[lb:ub]
  lo3 = predict(loess(fft3~time))[lb:ub]
  plot.ts(fft1)
  lines(lb:ub,lo1,col= "red")
  lines(lb:ub,lo2,col= "blue")
  lines(lb:ub,lo3,col= "green")
}

#' @export
saveModel.fftloess <- function(data.path){
  load(paste(data.path,"bearingDataH1_1.Rd",sep=""))
  h1_1=horizontal
  h1_1=as.numeric(unlist(h1_1))
  load(paste(data.path,"bearingDataH1_2.Rd",sep=""))
  h1_2=horizontal
  h1_2=as.numeric(unlist(h1_2))
  load(paste(data.path,"bearingDataH1_3.Rd",sep=""))
  h1_3=horizontal
  h1_3=as.numeric(unlist(h1_3))
  load(paste(data.path,"bearingDataH1_4.Rd",sep=""))
  h1_4=horizontal
  h1_4=as.numeric(unlist(h1_4))
  load(paste(data.path,"bearingDataH1_5.Rd",sep=""))
  h1_5=horizontal
  h1_5=as.numeric(unlist(h1_5))
  load(paste(data.path,"bearingDataH1_6.Rd",sep=""))
  h1_6=horizontal
  h1_6=as.numeric(unlist(h1_6))
  load(paste(data.path,"bearingDataH1_7.Rd",sep=""))
  h1_7=horizontal
  h1_7=as.numeric(unlist(h1_7))
  load(paste(data.path,"bearingDataH2_1.Rd",sep=""))
  h2_1=horizontal
  h2_1=as.numeric(unlist(h1_1))
  load(paste(data.path,"bearingDataH2_2.Rd",sep=""))
  h2_2=horizontal
  h2_2=as.numeric(unlist(h1_2))
  load(paste(data.path,"bearingDataH2_3.Rd",sep=""))
  h2_3=horizontal
  h2_3=as.numeric(unlist(h2_3))
  load(paste(data.path,"bearingDataH2_4.Rd",sep=""))
  h2_4=horizontal
  h2_4=as.numeric(unlist(h2_4))
  load(paste(data.path,"bearingDataH2_5.Rd",sep=""))
  h2_5=horizontal
  h2_5=as.numeric(unlist(h2_5))
  load(paste(data.path,"bearingDataH2_6.Rd",sep=""))
  h2_6=horizontal
  h2_6=as.numeric(unlist(h2_6))
  load(paste(data.path,"bearingDataH2_7.Rd",sep=""))
  h2_7=horizontal
  h2_7=as.numeric(unlist(h2_7))
  load(paste(data.path,"bearingDataH3_1.Rd",sep=""))
  h3_1=horizontal
  h3_1=as.numeric(unlist(h3_1))
  load(paste(data.path,"bearingDataH3_2.Rd",sep=""))
  h3_2=horizontal
  h3_2=as.numeric(unlist(h3_2))
  load(paste(data.path,"bearingDataH3_3.Rd",sep=""))
  h3_3=horizontal
  h3_3=as.numeric(unlist(h3_3))
  load(paste(data.path,"bearingDataV1_1.Rd",sep=""))
  v1_1=vertical
  v1_1=as.numeric(unlist(v1_1))
  load(paste(data.path,"bearingDataV1_2.Rd",sep=""))
  v1_2=horizontal
  v1_2=as.numeric(unlist(v1_2))
  load(paste(data.path,"bearingDataV1_3.Rd",sep=""))
  v1_3=vertical
  v1_3=as.numeric(unlist(v1_4))
  load(paste(data.path,"bearingDataV1_4.Rd",sep=""))
  v1_4=vertical
  v1_4=as.numeric(unlist(v1_4))
  load(paste(data.path,"bearingDataV1_5.Rd",sep=""))
  v1_5=vertical
  v1_5=as.numeric(unlist(v1_5))
  load(paste(data.path,"bearingDataV1_6.Rd",sep=""))
  v1_6=vertical
  v1_6=as.numeric(unlist(v1_6))
  load(paste(data.path,"bearingDataV1_7.Rd",sep=""))
  v1_7=vertical
  v1_7=as.numeric(unlist(v1_7))
  load(paste(data.path,"bearingDataV2_1.Rd",sep=""))
  v2_1=vertical
  v2_1=as.numeric(unlist(v2_1))
  load(paste(data.path,"bearingDataV2_2.Rd",sep=""))
  v2_2=horizontal
  v2_2=as.numeric(unlist(v2_2))
  load(paste(data.path,"bearingDataV2_3.Rd",sep=""))
  v2_3=vertical
  v2_3=as.numeric(unlist(v2_3))
  load(paste(data.path,"bearingDataV2_4.Rd",sep=""))
  v2_4=vertical
  v2_4=as.numeric(unlist(v2_4))
  load(paste(data.path,"bearingDataV2_5.Rd",sep=""))
  v2_5=vertical
  v2_5=as.numeric(unlist(v2_5))
  load(paste(data.path,"bearingDataV2_6.Rd",sep=""))
  v2_6=vertical
  v2_6=as.numeric(unlist(v2_6))
  load(paste(data.path,"bearingDataV2_7.Rd",sep=""))
  v2_7=vertical
  v2_7=as.numeric(unlist(v2_7))
  load(paste(data.path,"bearingDataV1_1.Rd",sep=""))
  v3_1=vertical
  v3_1=as.numeric(unlist(v1_1))
  load(paste(data.path,"bearingDataV1_2.Rd",sep=""))
  v3_2=horizontal
  v3_2=as.numeric(unlist(v1_2))
  load(paste(data.path,"bearingDataV3_3.Rd",sep=""))
  v3_3=vertical
  v3_3=as.numeric(unlist(v3_3))
  createModel.h2o.logic.save(h1_1,2048,0.8,"1_1")
  createModel.h2o.logic.save(h1_2,2048,0.8,"1_2")
  createModel.h2o.logic.save(h1_3,2048,0.8,"1_3")
  createModel.h2o.logic.save(h1_4,2048,0.8,"1_4")
  createModel.h2o.logic.save(h1_5,2048,0.8,"1_5")
  createModel.h2o.logic.save(h1_6,2048,0.8,"1_6")
  createModel.h2o.logic.save(h1_7,2048,0.8,"1_7")
  createModel.h2o.logic.save(h2_1,2048,0.8,"2_1")
  createModel.h2o.logic.save(h2_2,2048,0.8,"2_2")
  createModel.h2o.logic.save(h2_3,2048,0.8,"2_3")
  createModel.h2o.logic.save(h2_4,2048,0.8,"2_4")
  createModel.h2o.logic.save(h2_5,2048,0.8,"2_5")
  createModel.h2o.logic.save(h2_6,2048,0.8,"2_6")
  createModel.h2o.logic.save(h2_7,2048,0.8,"2_7")
  createModel.h2o.logic.save(h3_1,2048,0.8,"1_1")
  createModel.h2o.logic.save(h3_2,2048,0.8,"1_2")
  createModel.h2o.logic.save(h3_3,2048,0.8,"1_3")
  createModel.h2o.logic.save(v1_1,2048,0.8,"1_1")
  createModel.h2o.logic.save(v1_2,2048,0.8,"1_2")
  createModel.h2o.logic.save(v1_3,2048,0.8,"1_3")
  createModel.h2o.logic.save(v1_4,2048,0.8,"1_4")
  createModel.h2o.logic.save(v1_5,2048,0.8,"1_5")
  createModel.h2o.logic.save(v1_6,2048,0.8,"1_6")
  createModel.h2o.logic.save(v1_7,2048,0.8,"1_7")
  createModel.h2o.logic.save(v2_1,2048,0.8,"2_1")
  createModel.h2o.logic.save(v2_2,2048,0.8,"2_2")
  createModel.h2o.logic.save(v2_3,2048,0.8,"2_3")
  createModel.h2o.logic.save(v2_4,2048,0.8,"2_4")
  createModel.h2o.logic.save(v2_5,2048,0.8,"2_5")
  createModel.h2o.logic.save(v2_6,2048,0.8,"2_6")
  createModel.h2o.logic.save(v2_7,2048,0.8,"2_7")
  createModel.h2o.logic.save(v3_1,2048,0.8,"1_1")
  createModel.h2o.logic.save(v3_2,2048,0.8,"1_2")
  createModel.h2o.logic.save(v3_3,2048,0.8,"1_3")
}
