meanPSD <- function(data){
  result = sum(abs(fft(data)))/length(data)
  result = 20*log10(result/(10e-5))
  return(result)
}

pca <- function(data){
  data.pca = prcomp(data,scale. = TRUE,center = TRUE)
  indexMax = which.max(data.pca$sdev)
  return(data.pca$x[,indexMax])
}

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

entropyR <- function(data, binNum){
  data.dis = entropy::discretize(data,binNum)
  return(entropy::entropy(data.dis))
}

mutInfo <- function(data1, data2, binNum1, binNum2){
  data.dis = entropy::discretize2d(data1,data2,binNum1,binNum2)
  return(entropy::mi.plugin(data.dis))
}

featureEx <- function(dataH, dataV, binNumH, binNumV, groupSize, max.order){
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
  done = c(0,0)
  for(i in 0:(max.order-1)){
    index = max.order - i
    lenH = length(which(is.na(as.data.frame(arMH)[,index])))

    if((lenH == len)&&(done[1] == 0))
      arMH = arMH[,-index]
    else
      done[1] = 1

    lenV = length(which(is.na(arMV[,index])))
    if((lenV == len)&&(done[2] == 0))
      arMV = arMV[,-index]
    else
      done[2] = 1
    if((done[1]==1)&&(done[2]==1))
      break
  }
  arMH[is.na(arMH)] = 0
  arMV[is.na(arMV)] = 0
  # p2p = data.frame(Horizontal = peak2peakH, Vertical = peak2peakV)
  # maxPeak = data.frame(Horizontal = maxPeakH, Vertical = maxPeakV)
  # RMS = data.frame(Horizontal = RMSH, Vertical = RMSV)
  # Kurtosis = data.frame(Horizontal = KurtosisH, Vertical = KurtosisV)
  # Skewness = data.frame(Horizontal = SkewnessH, Vertical = SkewnessV)
  # Entropy = data.frame(Horizontal = EntropyH, Vertical = EntropyV)
  # arithmeticMeanPSD = data.frame(Horizontal = amPSDH, Vertical = amPSDV)
  # LineIntegral = data.frame(Horizontal = LIH, Vertical = LIV)
  # arModels = data.frame(Horizontal = arMH, Vertical = arMV)
  # Energy = data.frame(Horizontal = EnergyH, Vertical = EnergyV)
  result = list(peak2peakH,peak2peakV,maxPeakH,maxPeakV,RMSH,RMSV,KurtosisH,KurtosisV,SkewnessH,SkewnessV,MI,EntropyH,EntropyV,amPSDH,amPSDV,LIH,LIV,EnergyH,EnergyV)
  arResultH = as.list(as.data.frame(arMH))
  arResultV = as.list(as.data.frame(arMV))
  result = append(result,arResultH)
  result = append(result,arResultV)
  return(result)
}
# [1]pk2pkH [2]pk2pkV [3]MaxPeakH [4]MaxPeakV [5]rmsH [6]rmsV [7]KurtosisH [8]KurtosisV [9]SkewnessH [10]SkewnessV [11]MutInfo [12]EntropyH [13]EntropyV [14]amPSDH [15]amPSDV [16]LineIntH [17]LineIntV [18]EnergyH [19]EnergyV [20...]arModelH [21...]arModelV

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
