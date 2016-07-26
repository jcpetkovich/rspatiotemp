#' Compute Mean Arithmetic PSD of the data
#' @title Mean Arithmetic PSD (meanPSD)
#' @param data The data to have its mean arithmetic psd computed
#' @return The Mean Arithmetic PSD of the data
#' @export
meanPSD <- function(data){
  result = sum(abs(fft(data)))/length(data)
  result = 20*log10(result/(10e-5))
  return(result)
}

#' Compute Peak-to-peak value of the data
#' @title Peak-to-Peak (pk2pk)
#' @param data The data to have its peak-to-peak value computed
#' @return The Peak-to-Peak Value of the data
#' @export
pk2pk <- function(data){
  data.diff = diff(data)
  zeros = which(data.diff == 0)
  if (length(zeros) != 0){
    data = data[-zeros]
    data.diff = data.diff[-zeros]
  }
  peakMax = c(0,0)
  peakMin = c(0,0)
  len = length(data.diff)
  isPositive = 0
  if(data.diff[1] > 0)
    isPositive = 1
  for(i in 2:len){
    if(isPositive){
      if(data.diff[i] < 0){
        peakMax[1] = peakMax[1] + data[i]
        peakMax[2] = peakMax[2] + 1
        isPositive = 0
      }
    }
    else{
      if(data.diff[i] > 0){
        peakMin[1] = peakMin[1] + data[i]
        peakMin[2] = peakMin[2] + 1
        isPositive = 1
      }
    }
  }
  if(peakMax[2]!=0)
    peakMax = peakMax[1]/peakMax[2]
  if(peakMin[2]!=0)
    peakMin = peakMin[1]/peakMin[2]
  return(peakMax-peakMin)
}

#' Compute the AutoRegressive model values of the data
#' @title AutoRegressive model values (arModel)
#' @param data The data to have its autoregressive model values computed
#' @param order.max The maximum number of coefficients the ar model can compute
#' @return the AutoRegressive model values of the data
#' @export
arModel <- function(data, order.max){
  len = length(data)
  arCoef = ar(data,order.max = order.max)$ar
  arClen = length(arCoef)
  values = numeric()
  for(i in arClen:len){
    low = i - arClen + 1
    value = sum(data[low:i]*arCoef)
    values = c(values,value)
  }
  return(values)
}

#' Compute the entropy of the data
#' @title Entropy (entropyR)
#' @param data The data to have its entropy computed
#' @param binNum The binNum parameter for entropy::discretize function
#' @return The entropy value of the data
#' @export
entropyR <- function(data, binNum){
  data.dis = entropy::discretize(data,binNum)
  return(entropy::entropy(data.dis))
}

#' Compute Mutual Information value of the data
#' @title Mutual Information (mutInfo)
#' @param data1 The first set of data to have its mutual information computed
#' @param data2 The second set of data to have its mutual information computed
#' @param binNum1 The binNum parameter for entropy::discretize2d function for the first data set
#' @param binNum2 The binNum parameter for entropy::discretize2d function for the second data set
#' @return The Mutual Information value between the two data sets
#' @export
mutInfo <- function(data1, data2, binNum1, binNum2){
  data.dis = entropy::discretize2d(data1,data2,binNum1,binNum2)
  return(entropy::mi.plugin(data.dis))
}

#' An "all in one" function that computes several features together
#' Compute features (featureEx)
#' @param dataH The horizontal vibration data
#' @param dataV The vertical vibration data
#' @param binNumH The binNum parameter for entropy::discretize2d and entropy::discretize function
#' @param binNumV The binNum parameter for entropy::discretize2d and entropy::discretize function
#' @param groupSize The size of each group that will be split from the larger data for feature computation
#' @param max.order The maximum number of coefficients to be computed for the autoregressive model feature
#' @param ar.threshold The maximum percentage of NA values in the autoregressive model feature
#' @return A list of features computed from the data. [1]pk2pkH [2]pk2pkV [3]MaxPeakH [4]MaxPeakV [5]rmsH [6]rmsV [7]KurtosisH [8]KurtosisV [9]SkewnessH [10]SkewnessV [11]MutInfo [12]EntropyH [13]EntropyV [14]amPSDH [15]amPSDV [16]LineIntH [17]LineIntV [18]EnergyH [19]EnergyV [20...]arModelH [21...]arModelV
#' @export
featureEx <- function(dataH, dataV, binNumH, binNumV, groupSize, max.order, ar.threshold){
  boundSeq = seq(from=1,to = (length(dataH)+1), by = groupSize)
  len = length(boundSeq)-1
  peak2peakH = numeric()
  peak2peakV = numeric()
  maxPeakH = numeric()
  maxPeakV = numeric()
  RMSH = numeric()
  RMSV = numeric()
  KurtosisH = numeric()
  KurtosisV = numeric()
  SkewnessH = numeric()
  SkewnessV = numeric()
  amPSDH = numeric()
  amPSDV = numeric()
  EntropyH = numeric()
  EntropyV = numeric()
  LIH = numeric()
  LIV = numeric()
  EnergyH = numeric()
  EnergyV = numeric()
  MI = numeric()
  arMH = numeric()
  arMV = numeric()
  for(i in 1:len){
    lowerBound = boundSeq[i]
    upperBound = boundSeq[i+1] - 1
    dataHTemp = dataH[lowerBound:upperBound]
    dataVTemp = dataV[lowerBound:upperBound]
    partFeatures = allFeatures(dataHTemp,dataVTemp)

    peak2peakH = c(peak2peakH,partFeatures$pk2pkH)
    peak2peakV = c(peak2peakV,partFeatures$pk2pkV)

    maxPeakH = c(maxPeakH,max(dataHTemp))
    maxPeakV = c(maxPeakV,max(dataVTemp))

    RMSH = c(RMSH,partFeatures$rmsH)
    RMSV = c(RMSV,partFeatures$rmsV)

    KurtosisH = c(KurtosisH,partFeatures$kurtosisH)
    KurtosisV = c(KurtosisV,partFeatures$kurtosisV)

    SkewnessH = c(SkewnessH,partFeatures$skewH)
    SkewnessV = c(SkewnessV,partFeatures$skewV)

    MI = c(MI,mutInfo(dataHTemp,dataVTemp,binNumH,binNumV))

    EntropyH = c(EntropyH,entropyR(dataHTemp,binNumH))
    EntropyV = c(EntropyV,entropyR(dataVTemp,binNumV))

    amPSDH = c(amPSDH,meanPSD(dataHTemp))
    amPSDV = c(amPSDV,meanPSD(dataVTemp))

    LIH = c(LIH,partFeatures$lineIntH)
    LIV = c(LIV,partFeatures$lineIntV)

    if(length(arMH) == 0){
      arMH = ar(dataHTemp,order.max = max.order)$ar[1:max.order]
    }
    else{
      arMH = rbind(arMH,ar(dataHTemp,order.max = max.order)$ar[1:max.order])
    }

    if(length(arMV) == 0){
      arMV = ar(dataVTemp,order.max = max.order)$ar[1:max.order]
    }
    else{
      arMV = rbind(arMV,ar(dataVTemp,order.max = max.order)$ar[1:max.order])
    }

    EnergyH = c(EnergyH,partFeatures$energyH)
    EnergyV = c(EnergyV,partFeatures$energyV)
  }
  lenAR = length(arMV[,1])
  done = c(0,0)
  for(i in 1:max.order){
    if(done[1]==0){
      zerosH = length(which(is.na(arMH[,i])))/lenAR
      if(zerosH > ar.threshold){
        arMH = arMH[,-(i:max.order)]
        done[1] = 1
      }
    }
    if(done[2]==0){
      zerosV = length(which(is.na(arMV[,i])))/lenAR
      if(zerosV > ar.threshold){
        arMV = arMV[,-(i:max.order)]
        done[2] = 1
      }
    }
    if((done[1]==1)&&(done[2]==1))
      break
  }
  arMH[is.na(arMH)] = 0
  arMV[is.na(arMV)] = 0
  result = list(peak2peakH,peak2peakV,maxPeakH,maxPeakV,RMSH,RMSV,KurtosisH,KurtosisV,SkewnessH,SkewnessV,MI,EntropyH,EntropyV,amPSDH,amPSDV,LIH,LIV,EnergyH,EnergyV)
  arResultH = as.list(as.data.frame(arMH))
  arResultV = as.list(as.data.frame(arMV))
  result = append(result,arResultH)
  result = append(result,arResultV)
  return(result)
}
# [1]pk2pkH [2]pk2pkV [3]MaxPeakH [4]MaxPeakV [5]rmsH [6]rmsV [7]KurtosisH [8]KurtosisV [9]SkewnessH [10]SkewnessV [11]MutInfo [12]EntropyH [13]EntropyV [14]amPSDH [15]amPSDV [16]LineIntH [17]LineIntV [18]EnergyH [19]EnergyV [20...]arModelH [21...]arModelV

#' Compute the (1 - Symmetrical uncertainty) of the features
#' @title Symmetrical Uncertainty (SU)
#' @param features A list of features, the output of the 'featureEx' function
#' @param binNum The binNum parameter for entropy::discretize2d and entropy::discretize function
#' @return A matrix with the SU values of each feature against each other
#' @export
SU <- function(features,binNum){
  len = length(features)
  entropy = numeric()
  result = numeric()
  #c = 1
  #r = 1
  entropy = c(entropy, entropyR(as.numeric(unlist(features[1])),binNum))
  MI = mutInfo(as.numeric(unlist(features[1])),as.numeric(unlist(features[1])),binNum,binNum)
  res = MI/entropy[1]
  result = c(result,res)
  #c = 1
  for(r in 2:len){
    entropy = c(entropy, entropyR(as.numeric(unlist(features[r])),binNum))
    MI = mutInfo(as.numeric(unlist(features[r])),as.numeric(unlist(features[1])),binNum,binNum)
    res = 2*MI/(entropy[r]+entropy[1])
    result = c(result,res)
  }
  for(c in 2:len){
    for(r in 1:len){
      MI = mutInfo(as.numeric(unlist(features[r])),as.numeric(unlist(features[c])),binNum,binNum)
      res = 2*MI/(entropy[r]+entropy[c])
      result = c(result,res)
    }
  }
  result = matrix(result,len,len)
  return(1-result)
}

#' Compute the Root-Mean-Squared Error of the data given a linear model
#' @title Root-Mean-Squared Error (RMSE)
#' @param intercept The intercept of the linear model
#' @param coefficient The coefficient of the linear model
#' @param data The data for computing the RMSE
#' @param clusterNum The possible indices for the cluster having its RMSE computed
#' @return the RMSE of the data given the linear model
#' @export
RMSE <- function(intercept, coefficient, data, clusterNum){
  len = length(clusterNum)
  sum = 0
  for(i in 1:len){
    error = coefficient*clusterNum[i]+intercept-data[i]
    sum = sum + error*error
  }
  result = sqrt(sum/len)
  return(result)
}

#' Use the L-method to compute the knee (optimal number of clusters) for the hierarchical clustering
#' @title L-method (knee)
#' @param su The matrix of SU values, output of the 'SU' function
#' @return The knee (optimal number of clusters) of the hierarchical clustering
#' @export
knee <- function(su){
  hc = hclust(as.dist(su))
  len = length(hc$height)-1
  hc$height = rev(hc$height)
  RMSEc = numeric()
  for(i in 2:len){
    df1 = data.frame(height = hc$height[2:i], clusterNum = 2:i)
    l1 = lm(height~clusterNum, df1)
    rmse1 = RMSE(l1$coefficients[1],l1$coefficients[2],df1$height,df1$clusterNum)
    df2 = data.frame(height = hc$height[(i+1):len], clusterNum = (i+1):len)
    l2 = lm(height~clusterNum, df2)
    rmse2 = RMSE(l2$coefficients[1],l2$coefficients[2],df2$height,df2$clusterNum)
    RMSE = (i-1)*rmse1/len + (len-i+1)*rmse2/len
    RMSEc = c(RMSEc,RMSE)
  }
  return(which.min(RMSEc)+1)
}

#' Find the best cluster within the features using SOM
#' @title Choose Best Cluster (chooseCluster)
#' @param features A list of the computed features, output of the 'featureEx' function
#' @param su The SU values of the features, output of the 'SU' function
#' @param knee The knee of the clusters, output of the 'knee' function
#' @param xdim The parameter for the function som::som
#' @param ydim The parameter for the function som::som
#' @return The indices of the features of the chosen cluster
#' @export
chooseCluster <- function(features, su, knee, xdim, ydim){
  hc = hclust(as.dist(su))
  ct = cutree(hc,k=knee)
  qerrorScore = numeric()
  for(i in 1:knee){
    df = as.data.frame(features[which(ct==i)])
    df = unname(df)
    if(length(df)!=1){
      s = som(df,xdim,ydim)
      qerrorScore = c(qerrorScore,qerror(s))
    }
    else
      qerrorScore = c(qerrorScore, Inf)
  }
  return(which(ct==which.min(qerrorScore)))
}

#' Principal Component Analysis
#' @title Principal Component Analysis (pca)
#' @param features A list of the computed features, output of the 'featureEx' function
#' @param clustersChosen A vector containing the indices of the clusters chosen, output of the 'chooseCluster' function
#' @return The set of data from PCA with the lowest standard deviation
#' @export
pca <- function(features, clustersChosen){
  df = as.data.frame(features[clustersChosen])
  df = unname(df)
  pc = prcomp(df,scale. = TRUE)
  return(pc$x[,which.min(pc$sdev)])
}

#' Apply Empirical Mode Decomposition to the data
#' @title Empirical Mode Decomposition (EMDR)
#' @param data The data to be decomposed using EMD
#' @return The EMD object from the EMD R package
#' @export
emdR <- function(data){
  result = EMD::emd(data)
  return(result)
}

plotResidual <- function(emd){
  len = length(emd$residue)
  x = 1:len
  y = emd$residue
  df = data.frame(y,x)
  lo = loess(y~x,df)
  plot.ts(y)
  lines(predict(lo),col = "red",lwd = 2)
}

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

selfRUL.h2o.fftloess <- function(data,groupSize,degradStartPercent){
  fftloess = get.fftloess(data,groupSize)
  rul = c(rep(nrow(fftloess),as.integer(nrow(fftloess)*degradStartPercent)),seq(from=(nrow(fftloess)), to = 1, length.out = (nrow(fftloess))-as.integer(nrow(fftloess)*degradStartPercent)))
  fftloess.df = data.frame(rul = rul, fftloess)
  fftloess.hex = as.h2o(fftloess.df)
  fftloess.dl = h2o.deeplearning(x = 2:(ncol(fftloess)),y = 1,training_frame = fftloess.hex)
  predictions <- h2o.predict(fftloess.dl, fftloess.hex)
  predictions = as.data.frame(predictions)
  plot.ts(predictions$predict)
}

selfRUL.h2o.fftloess.logic <- function(data,groupSize,degradStartPercent){
  fftloess = get.fftloess(data,groupSize)
  rul = c(rep(TRUE,as.integer(nrow(fftloess)*degradStartPercent)),rep(FALSE,(nrow(fftloess))-as.integer(nrow(fftloess)*degradStartPercent)))
  fftloess.df = data.frame(rul = rul, fftloess)
  fftloess.hex = as.h2o(fftloess.df)
  fftloess.dl = h2o.deeplearning(x = 2:(ncol(fftloess)),y = 1,training_frame = fftloess.hex)
  predictions <- h2o.predict(fftloess.dl, fftloess.hex)
  predictions = as.data.frame(predictions)
  plot.ts(predictions$predict)
}

selfRUL.h2o.fftloess.logic.skip <- function(data,groupSize,degradStartPercent,smoothSampleNum){
  fftloess = get.fftloess.skip(data,groupSize,smoothSampleNum)
  #rul = c(rep(nrow(fftloess),as.integer(nrow(fftloess)*degradStartPercent)),seq(from=(nrow(fftloess)), to = 1, length.out = (nrow(fftloess))-as.integer(nrow(fftloess)*degradStartPercent)))
  rul = c(rep(TRUE,as.integer(nrow(fftloess)*degradStartPercent)),rep(FALSE,(nrow(fftloess))-as.integer(nrow(fftloess)*degradStartPercent)))
  fftloess.df = data.frame(rul = rul, fftloess)
  fftloess.hex = as.h2o(fftloess.df)
  fftloess.dl = h2o.deeplearning(x = 2:(ncol(fftloess)),y = 1,training_frame = fftloess.hex)
  predictions <- h2o.predict(fftloess.dl, fftloess.hex)
  predictions = as.data.frame(predictions)
  plot.ts(predictions$predict)
}

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
