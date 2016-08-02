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

#' #' @export
#' predict.status.h2o.logic <- function(dataH,dataV,groupSize, modelH,modelV,n){
#'   predictionsH = predict.h2o.logic(dataH,groupSize,modelH)
#'   predictionsV = predict.h2o.logic(dataV,groupSize,modelV)
#'   boundSeq = c(1:(n-1),seq(from=n,to = (length(predictionsH)+1), by = n))
#'   status = numeric()
#'   for(i in boundSeq){
#'     lo = numeric()
#'     if(i < n)
#'       lo = 1
#'     else
#'       lo = i - n + 1
#'
#'     predH = predictionsH[lo:i]
#'     predV = predictionsV[lo:i]
#'
#'     if(all(predH == TRUE)){
#'       if(all(predV == FALSE)){
#'         status = c(status,-1)
#'       }
#'       else{
#'         status = c(status,1)
#'       }
#'     }
#'     else if(all(predH == FALSE)){
#'       status = c(status,-1)
#'     }
#'     else{
#'       if(all(predV == TRUE)){
#'         status = c(status,1)
#'       }
#'       else if(all(predV == FALSE)){
#'         status = c(status,-1)
#'       }
#'       else{
#'         status = c(status,0)
#'       }
#'     }
#'   }
#'   return(status)
#' }

#' @export
predict.status.h2o.logic <- function(predictH,predictV,n){
  boundSeq = c(1:(n-1),seq(from=n,to = (length(predictH)), by = n))
  status = numeric()
  for(i in boundSeq){
    lo = numeric()
    if(i < n)
      lo = 1
    else
      lo = i - n + 1

    predH = predictH[lo:i]
    predV = predictV[lo:i]

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

#' #' @export
#' saveModel.fftloess <- function(data.path){
#'
#' }
#'
#' #' @export
#' predict.all.logic <- function(data.path,model.path){
#'   load(paste(data.path,"bearingDataH1_1.Rd",sep=""))
#'   h1_1=horizontal
#'   h1_1=as.numeric(unlist(h1_1))
#'   load(paste(data.path,"bearingDataH1_2.Rd",sep=""))
#'   h1_2=horizontal
#'   h1_2=as.numeric(unlist(h1_2))
#'   load(paste(data.path,"bearingDataH1_3.Rd",sep=""))
#'   h1_3=horizontal
#'   h1_3=as.numeric(unlist(h1_3))
#'   load(paste(data.path,"bearingDataH1_4.Rd",sep=""))
#'   h1_4=horizontal
#'   h1_4=as.numeric(unlist(h1_4))
#'   load(paste(data.path,"bearingDataH1_5.Rd",sep=""))
#'   h1_5=horizontal
#'   h1_5=as.numeric(unlist(h1_5))
#'   load(paste(data.path,"bearingDataH1_6.Rd",sep=""))
#'   h1_6=horizontal
#'   h1_6=as.numeric(unlist(h1_6))
#'   load(paste(data.path,"bearingDataH1_7.Rd",sep=""))
#'   h1_7=horizontal
#'   h1_7=as.numeric(unlist(h1_7))
#'   load(paste(data.path,"bearingDataH2_1.Rd",sep=""))
#'   h2_1=horizontal
#'   h2_1=as.numeric(unlist(h1_1))
#'   load(paste(data.path,"bearingDataH2_2.Rd",sep=""))
#'   h2_2=horizontal
#'   h2_2=as.numeric(unlist(h1_2))
#'   load(paste(data.path,"bearingDataH2_3.Rd",sep=""))
#'   h2_3=horizontal
#'   h2_3=as.numeric(unlist(h2_3))
#'   load(paste(data.path,"bearingDataH2_4.Rd",sep=""))
#'   h2_4=horizontal
#'   h2_4=as.numeric(unlist(h2_4))
#'   load(paste(data.path,"bearingDataH2_5.Rd",sep=""))
#'   h2_5=horizontal
#'   h2_5=as.numeric(unlist(h2_5))
#'   load(paste(data.path,"bearingDataH2_6.Rd",sep=""))
#'   h2_6=horizontal
#'   h2_6=as.numeric(unlist(h2_6))
#'   load(paste(data.path,"bearingDataH2_7.Rd",sep=""))
#'   h2_7=horizontal
#'   h2_7=as.numeric(unlist(h2_7))
#'   load(paste(data.path,"bearingDataH3_1.Rd",sep=""))
#'   h3_1=horizontal
#'   h3_1=as.numeric(unlist(h3_1))
#'   load(paste(data.path,"bearingDataH3_2.Rd",sep=""))
#'   h3_2=horizontal
#'   h3_2=as.numeric(unlist(h3_2))
#'   load(paste(data.path,"bearingDataH3_3.Rd",sep=""))
#'   h3_3=horizontal
#'   h3_3=as.numeric(unlist(h3_3))
#'   load(paste(data.path,"bearingDataV1_1.Rd",sep=""))
#'   v1_1=vertical
#'   v1_1=as.numeric(unlist(v1_1))
#'   load(paste(data.path,"bearingDataV1_2.Rd",sep=""))
#'   v1_2=horizontal
#'   v1_2=as.numeric(unlist(v1_2))
#'   load(paste(data.path,"bearingDataV1_3.Rd",sep=""))
#'   v1_3=vertical
#'   v1_3=as.numeric(unlist(v1_3))
#'   load(paste(data.path,"bearingDataV1_4.Rd",sep=""))
#'   v1_4=vertical
#'   v1_4=as.numeric(unlist(v1_4))
#'   load(paste(data.path,"bearingDataV1_5.Rd",sep=""))
#'   v1_5=vertical
#'   v1_5=as.numeric(unlist(v1_5))
#'   load(paste(data.path,"bearingDataV1_6.Rd",sep=""))
#'   v1_6=vertical
#'   v1_6=as.numeric(unlist(v1_6))
#'   load(paste(data.path,"bearingDataV1_7.Rd",sep=""))
#'   v1_7=vertical
#'   v1_7=as.numeric(unlist(v1_7))
#'   load(paste(data.path,"bearingDataV2_1.Rd",sep=""))
#'   v2_1=vertical
#'   v2_1=as.numeric(unlist(v2_1))
#'   load(paste(data.path,"bearingDataV2_2.Rd",sep=""))
#'   v2_2=horizontal
#'   v2_2=as.numeric(unlist(v2_2))
#'   load(paste(data.path,"bearingDataV2_3.Rd",sep=""))
#'   v2_3=vertical
#'   v2_3=as.numeric(unlist(v2_3))
#'   load(paste(data.path,"bearingDataV2_4.Rd",sep=""))
#'   v2_4=vertical
#'   v2_4=as.numeric(unlist(v2_4))
#'   load(paste(data.path,"bearingDataV2_5.Rd",sep=""))
#'   v2_5=vertical
#'   v2_5=as.numeric(unlist(v2_5))
#'   load(paste(data.path,"bearingDataV2_6.Rd",sep=""))
#'   v2_6=vertical
#'   v2_6=as.numeric(unlist(v2_6))
#'   load(paste(data.path,"bearingDataV2_7.Rd",sep=""))
#'   v2_7=vertical
#'   v2_7=as.numeric(unlist(v2_7))
#'   load(paste(data.path,"bearingDataV1_1.Rd",sep=""))
#'   v3_1=vertical
#'   v3_1=as.numeric(unlist(v1_1))
#'   load(paste(data.path,"bearingDataV1_2.Rd",sep=""))
#'   v3_2=horizontal
#'   v3_2=as.numeric(unlist(v1_2))
#'   load(paste(data.path,"bearingDataV3_3.Rd",sep=""))
#'   v3_3=vertical
#'   v3_3=as.numeric(unlist(v3_3))
#'
#'   modelH1_1 = createModel.h2o.logic(h1_1,2048,0.8)
#'   save(modelH1_1, file = "modelH1_1.Rd")
#'   modelH1_2 = createModel.h2o.logic(h1_2,2048,0.8)
#'   save(modelH1_2, file = "modelH1_2.Rd")
#'   modelH1_3 = createModel.h2o.logic(h1_3,2048,0.8)
#'   save(modelH1_3, file = "modelH1_3.Rd")
#'   modelH1_4 = createModel.h2o.logic(h1_4,2048,0.8)
#'   save(modelH1_4, file = "modelH1_4.Rd")
#'   modelH1_5 = createModel.h2o.logic(h1_5,2048,0.8)
#'   save(modelH1_5, file = "modelH1_5.Rd")
#'   modelH1_6 = createModel.h2o.logic(h1_6,2048,0.8)
#'   save(modelH1_6, file = "modelH1_6.Rd")
#'   modelH1_7 = createModel.h2o.logic(h1_7,2048,0.8)
#'   save(modelH1_7, file = "modelH1_7.Rd")
#'   modelH2_1 = createModel.h2o.logic(h2_1,2048,0.8)
#'   save(modelH2_1, file = "modelH2_1.Rd")
#'   modelH2_2 = createModel.h2o.logic(h2_2,2048,0.8)
#'   save(modelH2_2, file = "modelH2_2.Rd")
#'   modelH2_3 = createModel.h2o.logic(h2_3,2048,0.8)
#'   save(modelH2_3, file = "modelH2_3.Rd")
#'   modelH2_4 = createModel.h2o.logic(h2_4,2048,0.8)
#'   save(modelH2_4, file = "modelH2_4.Rd")
#'   modelH2_5 = createModel.h2o.logic(h2_5,2048,0.8)
#'   save(modelH2_5, file = "modelH2_5.Rd")
#'   modelH2_6 = createModel.h2o.logic(h2_6,2048,0.8)
#'   save(modelH2_6, file = "modelH2_6.Rd")
#'   modelH2_7 = createModel.h2o.logic(h2_7,2048,0.8)
#'   save(modelH2_7, file = "modelH2_7.Rd")
#'   modelH3_1 = createModel.h2o.logic(h3_1,2048,0.8)
#'   save(modelH3_1, file = "modelH3_1.Rd")
#'   modelH3_2 = createModel.h2o.logic(h3_2,2048,0.8)
#'   save(modelH3_2, file = "modelH3_2.Rd")
#'   modelH3_3 = createModel.h2o.logic(h3_3,2048,0.8)
#'   save(modelH3_2, file = "modelH3_2.Rd")
#'
#'   modelV1_1 = createModel.h2o.logic(v1_1,2048,0.8)
#'   save(modelV1_1, file = "modelV1_1.Rd")
#'   modelV1_2 = createModel.h2o.logic(v1_2,2048,0.8)
#'   save(modelV1_2, file = "modelV1_2.Rd")
#'   modelV1_3 = createModel.h2o.logic(v1_3,2048,0.8)
#'   save(modelV1_3, file = "modelV1_3.Rd")
#'   modelV1_4 = createModel.h2o.logic(v1_4,2048,0.8)
#'   save(modelV1_4, file = "modelV1_4.Rd")
#'   modelV1_5 = createModel.h2o.logic(v1_5,2048,0.8)
#'   save(modelV1_5, file = "modelV1_5.Rd")
#'   modelV1_6 = createModel.h2o.logic(v1_6,2048,0.8)
#'   save(modelV1_6, file = "modelV1_6.Rd")
#'   modelV1_7 = createModel.h2o.logic(v1_7,2048,0.8)
#'   save(modelV1_7, file = "modelV1_7.Rd")
#'   modelV2_1 = createModel.h2o.logic(v2_1,2048,0.8)
#'   save(modelV2_1, file = "modelV2_1.Rd")
#'   modelV2_2 = createModel.h2o.logic(v2_2,2048,0.8)
#'   save(modelV2_2, file = "modelV2_2.Rd")
#'   modelV2_3 = createModel.h2o.logic(v2_3,2048,0.8)
#'   save(modelV2_3, file = "modelV2_3.Rd")
#'   modelV2_4 = createModel.h2o.logic(v2_4,2048,0.8)
#'   save(modelV2_4, file = "modelV2_4.Rd")
#'   modelV2_5 = createModel.h2o.logic(v2_5,2048,0.8)
#'   save(modelV2_5, file = "modelV2_5.Rd")
#'   modelV2_6 = createModel.h2o.logic(v2_6,2048,0.8)
#'   save(modelV2_6, file = "modelV2_6.Rd")
#'   modelV2_7 = createModel.h2o.logic(v2_7,2048,0.8)
#'   save(modelV2_7, file = "modelV2_7.Rd")
#'   modelV3_1 = createModel.h2o.logic(v3_1,2048,0.8)
#'   save(modelV3_1, file = "modelV3_1.Rd")
#'   modelV3_2 = createModel.h2o.logic(v3_2,2048,0.8)
#'   save(modelV3_2, file = "modelV3_2.Rd")
#'   modelV3_3 = createModel.h2o.logic(v3_3,2048,0.8)
#'   save(modelV3_2, file = "modelV3_2.Rd")
#'
#'   predictH1_1 = list()
#'   predictH1_2 = list()
#'   predictH1_3 = list()
#'   predictH1_4 = list()
#'   predictH1_5 = list()
#'   predictH1_6 = list()
#'   predictH1_7 = list()
#'
#'   predictH1_1[[1]] = predict.h2o.logic(h1_1,2048,modelH1_1)
#'   predictH1_1[[2]] = predict.h2o.logic(h1_2,2048,modelH1_1)
#'   predictH1_1[[3]] = predict.h2o.logic(h1_3,2048,modelH1_1)
#'   predictH1_1[[4]] = predict.h2o.logic(h1_4,2048,modelH1_1)
#'   predictH1_1[[5]] = predict.h2o.logic(h1_5,2048,modelH1_1)
#'   predictH1_1[[6]] = predict.h2o.logic(h1_6,2048,modelH1_1)
#'   predictH1_1[[7]] = predict.h2o.logic(h1_7,2048,modelH1_1)
#'
#'   predictH1_2[[1]] = predict.h2o.logic(h1_1,2048,modelH1_2)
#'   predictH1_2[[2]] = predict.h2o.logic(h1_2,2048,modelH1_2)
#'   predictH1_2[[3]] = predict.h2o.logic(h1_3,2048,modelH1_2)
#'   predictH1_2[[4]] = predict.h2o.logic(h1_4,2048,modelH1_2)
#'   predictH1_2[[5]] = predict.h2o.logic(h1_5,2048,modelH1_2)
#'   predictH1_2[[6]] = predict.h2o.logic(h1_6,2048,modelH1_2)
#'   predictH1_2[[7]] = predict.h2o.logic(h1_7,2048,modelH1_2)
#'
#'   predictH1_3[[1]] = predict.h2o.logic(h1_1,2048,modelH1_3)
#'   predictH1_3[[2]] = predict.h2o.logic(h1_2,2048,modelH1_3)
#'   predictH1_3[[3]] = predict.h2o.logic(h1_3,2048,modelH1_3)
#'   predictH1_3[[4]] = predict.h2o.logic(h1_4,2048,modelH1_3)
#'   predictH1_3[[5]] = predict.h2o.logic(h1_5,2048,modelH1_3)
#'   predictH1_3[[6]] = predict.h2o.logic(h1_6,2048,modelH1_3)
#'   predictH1_3[[7]] = predict.h2o.logic(h1_7,2048,modelH1_3)
#'
#'   predictH1_4[[1]] = predict.h2o.logic(h1_1,2048,modelH1_4)
#'   predictH1_4[[2]] = predict.h2o.logic(h1_2,2048,modelH1_4)
#'   predictH1_4[[3]] = predict.h2o.logic(h1_3,2048,modelH1_4)
#'   predictH1_4[[4]] = predict.h2o.logic(h1_4,2048,modelH1_4)
#'   predictH1_4[[5]] = predict.h2o.logic(h1_5,2048,modelH1_4)
#'   predictH1_4[[6]] = predict.h2o.logic(h1_6,2048,modelH1_4)
#'   predictH1_4[[7]] = predict.h2o.logic(h1_7,2048,modelH1_4)
#'
#'   predictH1_5[[1]] = predict.h2o.logic(h1_1,2048,modelH1_5)
#'   predictH1_5[[2]] = predict.h2o.logic(h1_2,2048,modelH1_5)
#'   predictH1_5[[3]] = predict.h2o.logic(h1_3,2048,modelH1_5)
#'   predictH1_5[[4]] = predict.h2o.logic(h1_4,2048,modelH1_5)
#'   predictH1_5[[5]] = predict.h2o.logic(h1_5,2048,modelH1_5)
#'   predictH1_5[[6]] = predict.h2o.logic(h1_6,2048,modelH1_5)
#'   predictH1_5[[7]] = predict.h2o.logic(h1_7,2048,modelH1_5)
#'
#'   predictH1_6[[1]] = predict.h2o.logic(h1_1,2048,modelH1_6)
#'   predictH1_6[[2]] = predict.h2o.logic(h1_2,2048,modelH1_6)
#'   predictH1_6[[3]] = predict.h2o.logic(h1_3,2048,modelH1_6)
#'   predictH1_6[[4]] = predict.h2o.logic(h1_4,2048,modelH1_6)
#'   predictH1_6[[5]] = predict.h2o.logic(h1_5,2048,modelH1_6)
#'   predictH1_6[[6]] = predict.h2o.logic(h1_6,2048,modelH1_6)
#'   predictH1_6[[7]] = predict.h2o.logic(h1_7,2048,modelH1_6)
#'
#'   predictH1_7[[1]] = predict.h2o.logic(h1_1,2048,modelH1_7)
#'   predictH1_7[[2]] = predict.h2o.logic(h1_2,2048,modelH1_7)
#'   predictH1_7[[3]] = predict.h2o.logic(h1_3,2048,modelH1_7)
#'   predictH1_7[[4]] = predict.h2o.logic(h1_4,2048,modelH1_7)
#'   predictH1_7[[5]] = predict.h2o.logic(h1_5,2048,modelH1_7)
#'   predictH1_7[[6]] = predict.h2o.logic(h1_6,2048,modelH1_7)
#'   predictH1_7[[7]] = predict.h2o.logic(h1_7,2048,modelH1_7)
#'
#'   save(predictH1_1,file = "predictH1_1.Rd")
#'   save(predictH1_2,file = "predictH1_2.Rd")
#'   save(predictH1_3,file = "predictH1_3.Rd")
#'   save(predictH1_1,file = "predictH1_4.Rd")
#'   save(predictH1_2,file = "predictH1_5.Rd")
#'   save(predictH1_3,file = "predictH1_6.Rd")
#'   save(predictH1_1,file = "predictH1_7.Rd")
#'
#'   predictH2_1 = list()
#'   predictH2_2 = list()
#'   predictH2_3 = list()
#'   predictH2_4 = list()
#'   predictH2_5 = list()
#'   predictH2_6 = list()
#'   predictH2_7 = list()
#'
#'   predictH2_1[[1]] = predict.h2o.logic(h2_1,2048,modelH2_1)
#'   predictH2_1[[2]] = predict.h2o.logic(h2_2,2048,modelH2_1)
#'   predictH2_1[[3]] = predict.h2o.logic(h2_3,2048,modelH2_1)
#'   predictH2_1[[4]] = predict.h2o.logic(h2_4,2048,modelH2_1)
#'   predictH2_1[[5]] = predict.h2o.logic(h2_5,2048,modelH2_1)
#'   predictH2_1[[6]] = predict.h2o.logic(h2_6,2048,modelH2_1)
#'   predictH2_1[[7]] = predict.h2o.logic(h2_7,2048,modelH2_1)
#'
#'   predictH2_2[[1]] = predict.h2o.logic(h2_1,2048,modelH2_2)
#'   predictH2_2[[2]] = predict.h2o.logic(h2_2,2048,modelH2_2)
#'   predictH2_2[[3]] = predict.h2o.logic(h2_3,2048,modelH2_2)
#'   predictH2_2[[4]] = predict.h2o.logic(h2_4,2048,modelH2_2)
#'   predictH2_2[[5]] = predict.h2o.logic(h2_5,2048,modelH2_2)
#'   predictH2_2[[6]] = predict.h2o.logic(h2_6,2048,modelH2_2)
#'   predictH2_2[[7]] = predict.h2o.logic(h2_7,2048,modelH2_2)
#'
#'   predictH2_3[[1]] = predict.h2o.logic(h2_1,2048,modelH2_3)
#'   predictH2_3[[2]] = predict.h2o.logic(h2_2,2048,modelH2_3)
#'   predictH2_3[[3]] = predict.h2o.logic(h2_3,2048,modelH2_3)
#'   predictH2_3[[4]] = predict.h2o.logic(h2_4,2048,modelH2_3)
#'   predictH2_3[[5]] = predict.h2o.logic(h2_5,2048,modelH2_3)
#'   predictH2_3[[6]] = predict.h2o.logic(h2_6,2048,modelH2_3)
#'   predictH2_3[[7]] = predict.h2o.logic(h2_7,2048,modelH2_3)
#'
#'   predictH2_4[[1]] = predict.h2o.logic(h2_1,2048,modelH2_4)
#'   predictH2_4[[2]] = predict.h2o.logic(h2_2,2048,modelH2_4)
#'   predictH2_4[[3]] = predict.h2o.logic(h2_3,2048,modelH2_4)
#'   predictH2_4[[4]] = predict.h2o.logic(h2_4,2048,modelH2_4)
#'   predictH2_4[[5]] = predict.h2o.logic(h2_5,2048,modelH2_4)
#'   predictH2_4[[6]] = predict.h2o.logic(h2_6,2048,modelH2_4)
#'   predictH2_4[[7]] = predict.h2o.logic(h2_7,2048,modelH2_4)
#'
#'   predictH2_5[[1]] = predict.h2o.logic(h2_1,2048,modelH2_5)
#'   predictH2_5[[2]] = predict.h2o.logic(h2_2,2048,modelH2_5)
#'   predictH2_5[[3]] = predict.h2o.logic(h2_3,2048,modelH2_5)
#'   predictH2_5[[4]] = predict.h2o.logic(h2_4,2048,modelH2_5)
#'   predictH2_5[[5]] = predict.h2o.logic(h2_5,2048,modelH2_5)
#'   predictH2_5[[6]] = predict.h2o.logic(h2_6,2048,modelH2_5)
#'   predictH2_5[[7]] = predict.h2o.logic(h2_7,2048,modelH2_5)
#'
#'   predictH2_6[[1]] = predict.h2o.logic(h2_1,2048,modelH2_6)
#'   predictH2_6[[2]] = predict.h2o.logic(h2_2,2048,modelH2_6)
#'   predictH2_6[[3]] = predict.h2o.logic(h2_3,2048,modelH2_6)
#'   predictH2_6[[4]] = predict.h2o.logic(h2_4,2048,modelH2_6)
#'   predictH2_6[[5]] = predict.h2o.logic(h2_5,2048,modelH2_6)
#'   predictH2_6[[6]] = predict.h2o.logic(h2_6,2048,modelH2_6)
#'   predictH2_6[[7]] = predict.h2o.logic(h2_7,2048,modelH2_6)
#'
#'   predictH2_7[[1]] = predict.h2o.logic(h2_1,2048,modelH2_7)
#'   predictH2_7[[2]] = predict.h2o.logic(h2_2,2048,modelH2_7)
#'   predictH2_7[[3]] = predict.h2o.logic(h2_3,2048,modelH2_7)
#'   predictH2_7[[4]] = predict.h2o.logic(h2_4,2048,modelH2_7)
#'   predictH2_7[[5]] = predict.h2o.logic(h2_5,2048,modelH2_7)
#'   predictH2_7[[6]] = predict.h2o.logic(h2_6,2048,modelH2_7)
#'   predictH2_7[[7]] = predict.h2o.logic(h2_7,2048,modelH2_7)
#'
#'   save(predictH2_1,file = "predictH2_1.Rd")
#'   save(predictH2_2,file = "predictH2_2.Rd")
#'   save(predictH2_3,file = "predictH2_3.Rd")
#'   save(predictH2_1,file = "predictH2_4.Rd")
#'   save(predictH2_2,file = "predictH2_5.Rd")
#'   save(predictH2_3,file = "predictH2_6.Rd")
#'   save(predictH2_1,file = "predictH2_7.Rd")
#'
#'   predictH3_1 = list()
#'   predictH3_2 = list()
#'   predictH3_3 = list()
#'
#'   predictH3_1[[1]] = predict.h2o.logic(h3_1,2048,modelH3_1)
#'   predictH3_1[[2]] = predict.h2o.logic(h3_2,2048,modelH3_1)
#'   predictH3_1[[3]] = predict.h2o.logic(h3_3,2048,modelH3_1)
#'
#'   predictH3_2[[1]] = predict.h2o.logic(h3_1,2048,modelH3_2)
#'   predictH3_2[[2]] = predict.h2o.logic(h3_2,2048,modelH3_2)
#'   predictH3_2[[3]] = predict.h2o.logic(h3_3,2048,modelH3_2)
#'
#'   predictH3_3[[1]] = predict.h2o.logic(h3_1,2048,modelH3_3)
#'   predictH3_3[[2]] = predict.h2o.logic(h3_2,2048,modelH3_3)
#'   predictH3_3[[3]] = predict.h2o.logic(h3_3,2048,modelH3_3)
#'
#'   save(predictH3_1,file = "predictH3_1.Rd")
#'   save(predictH3_2,file = "predictH3_2.Rd")
#'   save(predictH3_3,file = "predictH3_3.Rd")
#'
#'   predictV1_1 = list()
#'   predictV1_2 = list()
#'   predictV1_3 = list()
#'   predictV1_4 = list()
#'   predictV1_5 = list()
#'   predictV1_6 = list()
#'   predictV1_7 = list()
#'
#'   predictV1_1[[1]] = predict.h2o.logic(v1_1,2048,modelV1_1)
#'   predictV1_1[[2]] = predict.h2o.logic(v1_2,2048,modelV1_1)
#'   predictV1_1[[3]] = predict.h2o.logic(v1_3,2048,modelV1_1)
#'   predictV1_1[[4]] = predict.h2o.logic(v1_4,2048,modelV1_1)
#'   predictV1_1[[5]] = predict.h2o.logic(v1_5,2048,modelV1_1)
#'   predictV1_1[[6]] = predict.h2o.logic(v1_6,2048,modelV1_1)
#'   predictV1_1[[7]] = predict.h2o.logic(v1_7,2048,modelV1_1)
#'
#'   predictV1_2[[1]] = predict.h2o.logic(v1_1,2048,modelV1_2)
#'   predictV1_2[[2]] = predict.h2o.logic(v1_2,2048,modelV1_2)
#'   predictV1_2[[3]] = predict.h2o.logic(v1_3,2048,modelV1_2)
#'   predictV1_2[[4]] = predict.h2o.logic(v1_4,2048,modelV1_2)
#'   predictV1_2[[5]] = predict.h2o.logic(v1_5,2048,modelV1_2)
#'   predictV1_2[[6]] = predict.h2o.logic(v1_6,2048,modelV1_2)
#'   predictV1_2[[7]] = predict.h2o.logic(v1_7,2048,modelV1_2)
#'
#'   predictV1_3[[1]] = predict.h2o.logic(v1_1,2048,modelV1_3)
#'   predictV1_3[[2]] = predict.h2o.logic(v1_2,2048,modelV1_3)
#'   predictV1_3[[3]] = predict.h2o.logic(v1_3,2048,modelV1_3)
#'   predictV1_3[[4]] = predict.h2o.logic(v1_4,2048,modelV1_3)
#'   predictV1_3[[5]] = predict.h2o.logic(v1_5,2048,modelV1_3)
#'   predictV1_3[[6]] = predict.h2o.logic(v1_6,2048,modelV1_3)
#'   predictV1_3[[7]] = predict.h2o.logic(v1_7,2048,modelV1_3)
#'
#'   predictV1_4[[1]] = predict.h2o.logic(v1_1,2048,modelV1_4)
#'   predictV1_4[[2]] = predict.h2o.logic(v1_2,2048,modelV1_4)
#'   predictV1_4[[3]] = predict.h2o.logic(v1_3,2048,modelV1_4)
#'   predictV1_4[[4]] = predict.h2o.logic(v1_4,2048,modelV1_4)
#'   predictV1_4[[5]] = predict.h2o.logic(v1_5,2048,modelV1_4)
#'   predictV1_4[[6]] = predict.h2o.logic(v1_6,2048,modelV1_4)
#'   predictV1_4[[7]] = predict.h2o.logic(v1_7,2048,modelV1_4)
#'
#'   predictV1_5[[1]] = predict.h2o.logic(v1_1,2048,modelV1_5)
#'   predictV1_5[[2]] = predict.h2o.logic(v1_2,2048,modelV1_5)
#'   predictV1_5[[3]] = predict.h2o.logic(v1_3,2048,modelV1_5)
#'   predictV1_5[[4]] = predict.h2o.logic(v1_4,2048,modelV1_5)
#'   predictV1_5[[5]] = predict.h2o.logic(v1_5,2048,modelV1_5)
#'   predictV1_5[[6]] = predict.h2o.logic(v1_6,2048,modelV1_5)
#'   predictV1_5[[7]] = predict.h2o.logic(v1_7,2048,modelV1_5)
#'
#'   predictV1_6[[1]] = predict.h2o.logic(v1_1,2048,modelV1_6)
#'   predictV1_6[[2]] = predict.h2o.logic(v1_2,2048,modelV1_6)
#'   predictV1_6[[3]] = predict.h2o.logic(v1_3,2048,modelV1_6)
#'   predictV1_6[[4]] = predict.h2o.logic(v1_4,2048,modelV1_6)
#'   predictV1_6[[5]] = predict.h2o.logic(v1_5,2048,modelV1_6)
#'   predictV1_6[[6]] = predict.h2o.logic(v1_6,2048,modelV1_6)
#'   predictV1_6[[7]] = predict.h2o.logic(v1_7,2048,modelV1_6)
#'
#'   predictV1_7[[1]] = predict.h2o.logic(v1_1,2048,modelV1_7)
#'   predictV1_7[[2]] = predict.h2o.logic(v1_2,2048,modelV1_7)
#'   predictV1_7[[3]] = predict.h2o.logic(v1_3,2048,modelV1_7)
#'   predictV1_7[[4]] = predict.h2o.logic(v1_4,2048,modelV1_7)
#'   predictV1_7[[5]] = predict.h2o.logic(v1_5,2048,modelV1_7)
#'   predictV1_7[[6]] = predict.h2o.logic(v1_6,2048,modelV1_7)
#'   predictV1_7[[7]] = predict.h2o.logic(v1_7,2048,modelV1_7)
#'
#'   save(predictV1_1,file = "predictV1_1.Rd")
#'   save(predictV1_2,file = "predictV1_2.Rd")
#'   save(predictV1_3,file = "predictV1_3.Rd")
#'   save(predictV1_1,file = "predictV1_4.Rd")
#'   save(predictV1_2,file = "predictV1_5.Rd")
#'   save(predictV1_3,file = "predictV1_6.Rd")
#'   save(predictV1_1,file = "predictV1_7.Rd")
#'
#'   predictV2_1 = list()
#'   predictV2_2 = list()
#'   predictV2_3 = list()
#'   predictV2_4 = list()
#'   predictV2_5 = list()
#'   predictV2_6 = list()
#'   predictV2_7 = list()
#'
#'   predictV2_1[[1]] = predict.h2o.logic(v2_1,2048,modelV2_1)
#'   predictV2_1[[2]] = predict.h2o.logic(v2_2,2048,modelV2_1)
#'   predictV2_1[[3]] = predict.h2o.logic(v2_3,2048,modelV2_1)
#'   predictV2_1[[4]] = predict.h2o.logic(v2_4,2048,modelV2_1)
#'   predictV2_1[[5]] = predict.h2o.logic(v2_5,2048,modelV2_1)
#'   predictV2_1[[6]] = predict.h2o.logic(v2_6,2048,modelV2_1)
#'   predictV2_1[[7]] = predict.h2o.logic(v2_7,2048,modelV2_1)
#'
#'   predictV2_2[[1]] = predict.h2o.logic(v2_1,2048,modelV2_2)
#'   predictV2_2[[2]] = predict.h2o.logic(v2_2,2048,modelV2_2)
#'   predictV2_2[[3]] = predict.h2o.logic(v2_3,2048,modelV2_2)
#'   predictV2_2[[4]] = predict.h2o.logic(v2_4,2048,modelV2_2)
#'   predictV2_2[[5]] = predict.h2o.logic(v2_5,2048,modelV2_2)
#'   predictV2_2[[6]] = predict.h2o.logic(v2_6,2048,modelV2_2)
#'   predictV2_2[[7]] = predict.h2o.logic(v2_7,2048,modelV2_2)
#'
#'   predictV2_3[[1]] = predict.h2o.logic(v2_1,2048,modelV2_3)
#'   predictV2_3[[2]] = predict.h2o.logic(v2_2,2048,modelV2_3)
#'   predictV2_3[[3]] = predict.h2o.logic(v2_3,2048,modelV2_3)
#'   predictV2_3[[4]] = predict.h2o.logic(v2_4,2048,modelV2_3)
#'   predictV2_3[[5]] = predict.h2o.logic(v2_5,2048,modelV2_3)
#'   predictV2_3[[6]] = predict.h2o.logic(v2_6,2048,modelV2_3)
#'   predictV2_3[[7]] = predict.h2o.logic(v2_7,2048,modelV2_3)
#'
#'   predictV2_4[[1]] = predict.h2o.logic(v2_1,2048,modelV2_4)
#'   predictV2_4[[2]] = predict.h2o.logic(v2_2,2048,modelV2_4)
#'   predictV2_4[[3]] = predict.h2o.logic(v2_3,2048,modelV2_4)
#'   predictV2_4[[4]] = predict.h2o.logic(v2_4,2048,modelV2_4)
#'   predictV2_4[[5]] = predict.h2o.logic(v2_5,2048,modelV2_4)
#'   predictV2_4[[6]] = predict.h2o.logic(v2_6,2048,modelV2_4)
#'   predictV2_4[[7]] = predict.h2o.logic(v2_7,2048,modelV2_4)
#'
#'   predictV2_5[[1]] = predict.h2o.logic(v2_1,2048,modelV2_5)
#'   predictV2_5[[2]] = predict.h2o.logic(v2_2,2048,modelV2_5)
#'   predictV2_5[[3]] = predict.h2o.logic(v2_3,2048,modelV2_5)
#'   predictV2_5[[4]] = predict.h2o.logic(v2_4,2048,modelV2_5)
#'   predictV2_5[[5]] = predict.h2o.logic(v2_5,2048,modelV2_5)
#'   predictV2_5[[6]] = predict.h2o.logic(v2_6,2048,modelV2_5)
#'   predictV2_5[[7]] = predict.h2o.logic(v2_7,2048,modelV2_5)
#'
#'   predictV2_6[[1]] = predict.h2o.logic(v2_1,2048,modelV2_6)
#'   predictV2_6[[2]] = predict.h2o.logic(v2_2,2048,modelV2_6)
#'   predictV2_6[[3]] = predict.h2o.logic(v2_3,2048,modelV2_6)
#'   predictV2_6[[4]] = predict.h2o.logic(v2_4,2048,modelV2_6)
#'   predictV2_6[[5]] = predict.h2o.logic(v2_5,2048,modelV2_6)
#'   predictV2_6[[6]] = predict.h2o.logic(v2_6,2048,modelV2_6)
#'   predictV2_6[[7]] = predict.h2o.logic(v2_7,2048,modelV2_6)
#'
#'   predictV2_7[[1]] = predict.h2o.logic(v2_1,2048,modelV2_7)
#'   predictV2_7[[2]] = predict.h2o.logic(v2_2,2048,modelV2_7)
#'   predictV2_7[[3]] = predict.h2o.logic(v2_3,2048,modelV2_7)
#'   predictV2_7[[4]] = predict.h2o.logic(v2_4,2048,modelV2_7)
#'   predictV2_7[[5]] = predict.h2o.logic(v2_5,2048,modelV2_7)
#'   predictV2_7[[6]] = predict.h2o.logic(v2_6,2048,modelV2_7)
#'   predictV2_7[[7]] = predict.h2o.logic(v2_7,2048,modelV2_7)
#'
#'   save(predictV2_1,file = "predictV2_1.Rd")
#'   save(predictV2_2,file = "predictV2_2.Rd")
#'   save(predictV2_3,file = "predictV2_3.Rd")
#'   save(predictV2_1,file = "predictV2_4.Rd")
#'   save(predictV2_2,file = "predictV2_5.Rd")
#'   save(predictV2_3,file = "predictV2_6.Rd")
#'   save(predictV2_1,file = "predictV2_7.Rd")
#'
#'   predictV3_1 = list()
#'   predictV3_2 = list()
#'   predictV3_3 = list()
#'
#'   predictV3_1[[1]] = predict.h2o.logic(v3_1,2048,modelV3_1)
#'   predictV3_1[[2]] = predict.h2o.logic(v3_2,2048,modelV3_1)
#'   predictV3_1[[3]] = predict.h2o.logic(v3_3,2048,modelV3_1)
#'
#'   predictV3_2[[1]] = predict.h2o.logic(v3_1,2048,modelV3_2)
#'   predictV3_2[[2]] = predict.h2o.logic(v3_2,2048,modelV3_2)
#'   predictV3_2[[3]] = predict.h2o.logic(v3_3,2048,modelV3_2)
#'
#'   predictV3_3[[1]] = predict.h2o.logic(v3_1,2048,modelV3_3)
#'   predictV3_3[[2]] = predict.h2o.logic(v3_2,2048,modelV3_3)
#'   predictV3_3[[3]] = predict.h2o.logic(v3_3,2048,modelV3_3)
#'
#'   save(predictV3_1,file = "predictV3_1.Rd")
#'   save(predictV3_2,file = "predictV3_2.Rd")
#'   save(predictV3_3,file = "predictV3_3.Rd")
#' }
#'
#' load("modelH1_1.Rd")
#' load("modelH1_2.Rd")
#' load("modelH1_3.Rd")
#' load("modelH1_4.Rd")
#' load("modelH1_5.Rd")
#' load("modelH1_6.Rd")
#' load("modelH1_7.Rd")
#' load("modelH2_1.Rd")
#' load("modelH2_2.Rd")
#' load("modelH2_3.Rd")
#' load("modelH2_4.Rd")
#' load("modelH2_5.Rd")
#' load("modelH2_6.Rd")
#' load("modelH2_7.Rd")
#' load("modelH3_1.Rd")
#' load("modelH3_2.Rd")
#' load("modelH3_3.Rd")
#' load("modelV1_1.Rd")
#' load("modelV1_2.Rd")
#' load("modelV1_3.Rd")
#' load("modelV1_4.Rd")
#' load("modelV1_5.Rd")
#' load("modelV1_6.Rd")
#' load("modelV1_7.Rd")
#' load("modelV2_1.Rd")
#' load("modelV2_2.Rd")
#' load("modelV2_3.Rd")
#' load("modelV2_4.Rd")
#' load("modelV2_5.Rd")
#' load("modelV2_6.Rd")
#' load("modelV2_7.Rd")
#' load("modelV3_1.Rd")
#' load("modelV3_2.Rd")
#' load("modelV3_3.Rd")

load(paste(data.path,"predictH1_1.Rd",sep=""))
load(paste(data.path,"predictV1_1.Rd",sep=""))
status1_1 = list()
status1_1[[1]] = predict.status.h2o.logic(predictH1_1[[1]],predictV1_1[[1]],50)
status1_1[[2]] = predict.status.h2o.logic(predictH1_1[[2]],predictV1_1[[2]],50)
status1_1[[3]] = predict.status.h2o.logic(predictH1_1[[3]],predictV1_1[[3]],50)
status1_1[[4]] = predict.status.h2o.logic(predictH1_1[[4]],predictV1_1[[4]],50)
status1_1[[5]] = predict.status.h2o.logic(predictH1_1[[5]],predictV1_1[[5]],50)
status1_1[[6]] = predict.status.h2o.logic(predictH1_1[[6]],predictV1_1[[6]],50)
status1_1[[7]] = predict.status.h2o.logic(predictH1_1[[7]],predictV1_1[[7]],50)
save(status1_1,file = "status1_1.Rd")
load(paste(data.path,"predictH1_2.Rd",sep=""))
load(paste(data.path,"predictV1_2.Rd",sep=""))
status1_2 = list()
status1_2[[1]] = predict.status.h2o.logic(predictH1_2[[1]],predictV1_2[[1]],50)
status1_2[[2]] = predict.status.h2o.logic(predictH1_2[[2]],predictV1_2[[2]],50)
status1_2[[3]] = predict.status.h2o.logic(predictH1_2[[3]],predictV1_2[[3]],50)
status1_2[[4]] = predict.status.h2o.logic(predictH1_2[[4]],predictV1_2[[4]],50)
status1_2[[5]] = predict.status.h2o.logic(predictH1_2[[5]],predictV1_2[[5]],50)
status1_2[[6]] = predict.status.h2o.logic(predictH1_2[[6]],predictV1_2[[6]],50)
status1_2[[7]] = predict.status.h2o.logic(predictH1_2[[7]],predictV1_2[[7]],50)
save(status1_2,file = "status1_2.Rd")
load(paste(data.path,"predictH1_3.Rd",sep=""))
load(paste(data.path,"predictV1_3.Rd",sep=""))
status1_3 = list()
status1_3[[1]] = predict.status.h2o.logic(predictH1_3[[1]],predictV1_3[[1]],50)
status1_3[[2]] = predict.status.h2o.logic(predictH1_3[[2]],predictV1_3[[2]],50)
status1_3[[3]] = predict.status.h2o.logic(predictH1_3[[3]],predictV1_3[[3]],50)
status1_3[[4]] = predict.status.h2o.logic(predictH1_3[[4]],predictV1_3[[4]],50)
status1_3[[5]] = predict.status.h2o.logic(predictH1_3[[5]],predictV1_3[[5]],50)
status1_3[[6]] = predict.status.h2o.logic(predictH1_3[[6]],predictV1_3[[6]],50)
status1_3[[7]] = predict.status.h2o.logic(predictH1_3[[7]],predictV1_3[[7]],50)
save(status1_3,file = "status1_3.Rd")
load(paste(data.path,"predictH1_4.Rd",sep=""))
load(paste(data.path,"predictV1_4.Rd",sep=""))
status1_4 = list()
status1_4[[1]] = predict.status.h2o.logic(predictH1_4[[1]],predictV1_4[[1]],50)
status1_4[[2]] = predict.status.h2o.logic(predictH1_4[[2]],predictV1_4[[2]],50)
status1_4[[3]] = predict.status.h2o.logic(predictH1_4[[3]],predictV1_4[[3]],50)
status1_4[[4]] = predict.status.h2o.logic(predictH1_4[[4]],predictV1_4[[4]],50)
status1_4[[5]] = predict.status.h2o.logic(predictH1_4[[5]],predictV1_4[[5]],50)
status1_4[[6]] = predict.status.h2o.logic(predictH1_4[[6]],predictV1_4[[6]],50)
status1_4[[7]] = predict.status.h2o.logic(predictH1_4[[7]],predictV1_4[[7]],50)
save(status1_4,file = "status1_4.Rd")
load(paste(data.path,"predictH1_5.Rd",sep=""))
load(paste(data.path,"predictV1_5.Rd",sep=""))
status1_5 = list()
status1_5[[1]] = predict.status.h2o.logic(predictH1_5[[1]],predictV1_5[[1]],50)
status1_5[[2]] = predict.status.h2o.logic(predictH1_5[[2]],predictV1_5[[2]],50)
status1_5[[3]] = predict.status.h2o.logic(predictH1_5[[3]],predictV1_5[[3]],50)
status1_5[[4]] = predict.status.h2o.logic(predictH1_5[[4]],predictV1_5[[4]],50)
status1_5[[5]] = predict.status.h2o.logic(predictH1_5[[5]],predictV1_5[[5]],50)
status1_5[[6]] = predict.status.h2o.logic(predictH1_5[[6]],predictV1_5[[6]],50)
status1_5[[7]] = predict.status.h2o.logic(predictH1_5[[7]],predictV1_5[[7]],50)
save(status1_5,file = "status1_5.Rd")
load(paste(data.path,"predictH1_6.Rd",sep=""))
load(paste(data.path,"predictV1_6.Rd",sep=""))
status1_6 = list()
status1_6[[1]] = predict.status.h2o.logic(predictH1_6[[1]],predictV1_6[[1]],50)
status1_6[[2]] = predict.status.h2o.logic(predictH1_6[[2]],predictV1_6[[2]],50)
status1_6[[3]] = predict.status.h2o.logic(predictH1_6[[3]],predictV1_6[[3]],50)
status1_6[[4]] = predict.status.h2o.logic(predictH1_6[[4]],predictV1_6[[4]],50)
status1_6[[5]] = predict.status.h2o.logic(predictH1_6[[5]],predictV1_6[[5]],50)
status1_6[[6]] = predict.status.h2o.logic(predictH1_6[[6]],predictV1_6[[6]],50)
status1_6[[7]] = predict.status.h2o.logic(predictH1_6[[7]],predictV1_6[[7]],50)
save(status1_6,file = "status1_6.Rd")
load(paste(data.path,"predictH1_7.Rd",sep=""))
load(paste(data.path,"predictV1_7.Rd",sep=""))
status1_1 = list()
status1_7[[1]] = predict.status.h2o.logic(predictH1_7[[1]],predictV1_7[[1]],50)
status1_7[[2]] = predict.status.h2o.logic(predictH1_7[[2]],predictV1_7[[2]],50)
status1_7[[3]] = predict.status.h2o.logic(predictH1_7[[3]],predictV1_7[[3]],50)
status1_7[[4]] = predict.status.h2o.logic(predictH1_7[[4]],predictV1_7[[4]],50)
status1_7[[5]] = predict.status.h2o.logic(predictH1_7[[5]],predictV1_7[[5]],50)
status1_7[[6]] = predict.status.h2o.logic(predictH1_7[[6]],predictV1_7[[6]],50)
status1_7[[7]] = predict.status.h2o.logic(predictH1_7[[7]],predictV1_7[[7]],50)
save(status1_7,file = "status1_7.Rd")

load(paste(data.path,"predictH2_1.Rd",sep=""))
load(paste(data.path,"predictV2_1.Rd",sep=""))
status2_1 = list()
status2_1[[1]] = predict.status.h2o.logic(predictH2_1[[1]],predictV2_1[[1]],50)
status2_1[[2]] = predict.status.h2o.logic(predictH2_1[[2]],predictV2_1[[2]],50)
status2_1[[3]] = predict.status.h2o.logic(predictH2_1[[3]],predictV2_1[[3]],50)
status2_1[[4]] = predict.status.h2o.logic(predictH2_1[[4]],predictV2_1[[4]],50)
status2_1[[5]] = predict.status.h2o.logic(predictH2_1[[5]],predictV2_1[[5]],50)
status2_1[[6]] = predict.status.h2o.logic(predictH2_1[[6]],predictV2_1[[6]],50)
status2_1[[7]] = predict.status.h2o.logic(predictH2_1[[7]],predictV2_1[[7]],50)
save(status2_1,file = "status2_1.Rd")
load(paste(data.path,"predictH2_2.Rd",sep=""))
load(paste(data.path,"predictV2_2.Rd",sep=""))
status2_2 = list()
status2_2[[1]] = predict.status.h2o.logic(predictH2_2[[1]],predictV2_2[[1]],50)
status2_2[[2]] = predict.status.h2o.logic(predictH2_2[[2]],predictV2_2[[2]],50)
status2_2[[3]] = predict.status.h2o.logic(predictH2_2[[3]],predictV2_2[[3]],50)
status2_2[[4]] = predict.status.h2o.logic(predictH2_2[[4]],predictV2_2[[4]],50)
status2_2[[5]] = predict.status.h2o.logic(predictH2_2[[5]],predictV2_2[[5]],50)
status2_2[[6]] = predict.status.h2o.logic(predictH2_2[[6]],predictV2_2[[6]],50)
status2_2[[7]] = predict.status.h2o.logic(predictH2_2[[7]],predictV2_2[[7]],50)
save(status2_2,file = "status2_2.Rd")
load(paste(data.path,"predictH2_3.Rd",sep=""))
load(paste(data.path,"predictV2_3.Rd",sep=""))
status2_3 = list()
status2_3[[1]] = predict.status.h2o.logic(predictH2_3[[1]],predictV2_3[[1]],50)
status2_3[[2]] = predict.status.h2o.logic(predictH2_3[[2]],predictV2_3[[2]],50)
status2_3[[3]] = predict.status.h2o.logic(predictH2_3[[3]],predictV2_3[[3]],50)
status2_3[[4]] = predict.status.h2o.logic(predictH2_3[[4]],predictV2_3[[4]],50)
status2_3[[5]] = predict.status.h2o.logic(predictH2_3[[5]],predictV2_3[[5]],50)
status2_3[[6]] = predict.status.h2o.logic(predictH2_3[[6]],predictV2_3[[6]],50)
status2_3[[7]] = predict.status.h2o.logic(predictH2_3[[7]],predictV2_3[[7]],50)
save(status2_3,file = "status2_3.Rd")
load(paste(data.path,"predictH2_4.Rd",sep=""))
load(paste(data.path,"predictV2_4.Rd",sep=""))
status2_4 = list()
status2_4[[1]] = predict.status.h2o.logic(predictH2_4[[1]],predictV2_4[[1]],50)
status2_4[[2]] = predict.status.h2o.logic(predictH2_4[[2]],predictV2_4[[2]],50)
status2_4[[3]] = predict.status.h2o.logic(predictH2_4[[3]],predictV2_4[[3]],50)
status2_4[[4]] = predict.status.h2o.logic(predictH2_4[[4]],predictV2_4[[4]],50)
status2_4[[5]] = predict.status.h2o.logic(predictH2_4[[5]],predictV2_4[[5]],50)
status2_4[[6]] = predict.status.h2o.logic(predictH2_4[[6]],predictV2_4[[6]],50)
status2_4[[7]] = predict.status.h2o.logic(predictH2_4[[7]],predictV2_4[[7]],50)
save(status2_4,file = "status2_4.Rd")
load(paste(data.path,"predictH2_5.Rd",sep=""))
load(paste(data.path,"predictV2_5.Rd",sep=""))
status2_5 = list()
status2_5[[1]] = predict.status.h2o.logic(predictH2_5[[1]],predictV2_5[[1]],50)
status2_5[[2]] = predict.status.h2o.logic(predictH2_5[[2]],predictV2_5[[2]],50)
status2_5[[3]] = predict.status.h2o.logic(predictH2_5[[3]],predictV2_5[[3]],50)
status2_5[[4]] = predict.status.h2o.logic(predictH2_5[[4]],predictV2_5[[4]],50)
status2_5[[5]] = predict.status.h2o.logic(predictH2_5[[5]],predictV2_5[[5]],50)
status2_5[[6]] = predict.status.h2o.logic(predictH2_5[[6]],predictV2_5[[6]],50)
status2_5[[7]] = predict.status.h2o.logic(predictH2_5[[7]],predictV2_5[[7]],50)
save(status2_5,file = "status2_5.Rd")
load(paste(data.path,"predictH2_6.Rd",sep=""))
load(paste(data.path,"predictV2_6.Rd",sep=""))
status2_6 = list()
status2_6[[1]] = predict.status.h2o.logic(predictH2_6[[1]],predictV2_6[[1]],50)
status2_6[[2]] = predict.status.h2o.logic(predictH2_6[[2]],predictV2_6[[2]],50)
status2_6[[3]] = predict.status.h2o.logic(predictH2_6[[3]],predictV2_6[[3]],50)
status2_6[[4]] = predict.status.h2o.logic(predictH2_6[[4]],predictV2_6[[4]],50)
status2_6[[5]] = predict.status.h2o.logic(predictH2_6[[5]],predictV2_6[[5]],50)
status2_6[[6]] = predict.status.h2o.logic(predictH2_6[[6]],predictV2_6[[6]],50)
status2_6[[7]] = predict.status.h2o.logic(predictH2_6[[7]],predictV2_6[[7]],50)
save(status2_6,file = "status2_6.Rd")
load(paste(data.path,"predictH2_7.Rd",sep=""))
load(paste(data.path,"predictV2_7.Rd",sep=""))
status2_1 = list()
status2_7[[1]] = predict.status.h2o.logic(predictH2_7[[1]],predictV2_7[[1]],50)
status2_7[[2]] = predict.status.h2o.logic(predictH2_7[[2]],predictV2_7[[2]],50)
status2_7[[3]] = predict.status.h2o.logic(predictH2_7[[3]],predictV2_7[[3]],50)
status2_7[[4]] = predict.status.h2o.logic(predictH2_7[[4]],predictV2_7[[4]],50)
status2_7[[5]] = predict.status.h2o.logic(predictH2_7[[5]],predictV2_7[[5]],50)
status2_7[[6]] = predict.status.h2o.logic(predictH2_7[[6]],predictV2_7[[6]],50)
status2_7[[7]] = predict.status.h2o.logic(predictH2_7[[7]],predictV2_7[[7]],50)
save(status2_7,file = "status2_7.Rd")

load(paste(data.path,"predictH3_1.Rd",sep=""))
load(paste(data.path,"predictV3_1.Rd",sep=""))
status3_1 = list()
status3_1[[1]] = predict.status.h2o.logic(predictH3_1[[1]],predictV3_1[[1]],50)
status3_1[[2]] = predict.status.h2o.logic(predictH3_1[[2]],predictV3_1[[2]],50)
status3_1[[3]] = predict.status.h2o.logic(predictH3_1[[3]],predictV3_1[[3]],50)
save(status3_1,file = "status3_1.Rd")
load(paste(data.path,"predictH3_2.Rd",sep=""))
load(paste(data.path,"predictV3_2.Rd",sep=""))
status3_2 = list()
status3_2[[1]] = predict.status.h2o.logic(predictH3_2[[1]],predictV3_2[[1]],50)
status3_2[[2]] = predict.status.h2o.logic(predictH3_2[[2]],predictV3_2[[2]],50)
status3_2[[3]] = predict.status.h2o.logic(predictH3_2[[3]],predictV3_2[[3]],50)
save(status3_2,file = "status3_2.Rd")
load(paste(data.path,"predictH3_3.Rd",sep=""))
load(paste(data.path,"predictV3_3.Rd",sep=""))
status3_3 = list()
status3_3[[1]] = predict.status.h2o.logic(predictH3_3[[1]],predictV3_3[[1]],50)
status3_3[[2]] = predict.status.h2o.logic(predictH3_3[[2]],predictV3_3[[2]],50)
status3_3[[3]] = predict.status.h2o.logic(predictH3_3[[3]],predictV3_3[[3]],50)
save(status3_3,file = "status3_3.Rd")
