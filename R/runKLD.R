#' Convert time series data to SAX, group Train and Test sequences, compute each distribution then graph the KL divergence between the two sequences
#' @title KL Divergence of a single time series (runKLDVec)
#' @param rawData A vector of the raw data/time series data
#' @param saxGroupSize The group size when converting data into SAX representation
#' @param saxAlphabetSize The alphabet size when converting data into SAX representation
#' @param dMarkovSize The size when grouping observed data to predict the next hidden state
#' @param disStep The step size of sampling when computing the distribution
#' @param disWindowSize The size of the sample when computing the distribution
#' @param k The maximum number of nearest neighbours when computing the KL divergence
#' @return This function does not return anything. However it does plot a graph of the KLD
#' @export
runKLDVec <- function(rawData, saxGroupSize, saxAlphabetSize, dMarkovSize, disStep, disWindowSize, k){
  dataSAX = rspatiotemp::runSAX(rawData, saxGroupSize, saxAlphabetSize)
  dataO = convert(dataSAX, dMarkovSize, saxAlphabetSize)
  dataH = dataSAX[dMarkovSize:length(dataSAX)]
  halfway = length(dataH)/2
  dataOTr = dataO[0:halfway]
  dataOTe = dataO[halfway:length(dataH)]
  dataHTr = dataH[0:halfway]
  dataHTe = dataH[halfway:length(dataH)]
  dataCount = .Call('rspatiotemp_count', PACKAGE = 'rspatiotemp', dataOTr,dataHTr,saxAlphabetSize^dMarkovSize,saxAlphabetSize)
  disTrTr = .Call('rspatiotemp_distribution', PACKAGE = 'rspatiotemp', dataCount, dataOTr, dataHTr, saxAlphabetSize^dMarkovSize, saxAlphabetSize, disStep, disWindowSize)
  disTrTe = .Call('rspatiotemp_distribution', PACKAGE = 'rspatiotemp', dataCount, dataOTe, dataHTe, saxAlphabetSize^dMarkovSize, saxAlphabetSize, disStep, disWindowSize)
  plot(FNN::KL.divergence(10^disTrTr,10^disTrTe,k=k))
}

#' Convert time series data to SAX, simulate test sequences using HMMs, compute each distribution then graph the KL divergence between the two sequences
#' @title KL Divergence of a single time series and its simulation (runKLDVecSim)
#' @param rawData A vector of the raw data/time series data
#' @param saxGroupSize The group size when converting data into SAX representation
#' @param saxAlphabetSize The alphabet size when converting data into SAX representation
#' @param delay The number of zeros to pad on the hidden sequence before converting to SAX in order to delay it
#' @param dMarkovSize The size when grouping observed data to predict the next hidden state
#' @param disStep The step size of sampling when computing the distribution
#' @param disWindowSize The size of the sample when computing the distribution
#' @param k The maximum number of nearest neighbours when computing the KL divergence
#' @return This function does not return anything. However it does plot a graph of the KLD
#' @export
runKLDVecSim <- function(rawData, saxGroupSize, saxAlphabetSize, delay, dMarkovSize, disStep, disWindowSize, k){
  dataSAX = rspatiotemp::convert(rspatiotemp::runSAX(rawData, saxGroupSize, saxAlphabetSize),dMarkovSize,saxAlphabetSize)
  dataSAXD = rspatiotemp::runSAX(c(rep(0,delay),rawData), saxGroupSize, saxAlphabetSize)
  obsLen = length(dataSAX)
  dataProbMatX = rspatiotemp::createProbMatX(dataSAX,dataSAXD[0:obsLen],dMarkovSize,saxAlphabetSize^dMarkovSize,saxAlphabetSize)
  dataSim = rspatiotemp::simulateHid(dataProbMatX,dataSAX,dMarkovSize,saxAlphabetSize^dMarkovSize,saxAlphabetSize)

  countData = rspatiotemp::count(dataSAX,dataSAXD,saxAlphabetSize^dMarkovSize,saxAlphabetSize)
  dataDis = rspatiotemp::distribution(countData,dataSAX,dataSAXD,saxAlphabetSize^dMarkovSize,saxAlphabetSize,disStep,disWindowSize)
  dataSimDis = rspatiotemp::distribution(countData,dataSAX,dataSim,saxAlphabetSize^dMarkovSize,saxAlphabetSize,disStep,disWindowSize)

  plot(FNN::KL.divergence(10^dataDis,10^dataSimDis,k=k))
}

#' Convert time series data to SAX, compute each distribution then graph the KL divergence between the two sequences
#' @title KL Divergence of a two time series (runKLDVecVec)
#' @param rawData1 A vector of the first set of raw data/time series data
#' @param rawData2 A vector of the first set of raw data/time series data
#' @param saxGroupSize The group size when converting data into SAX representation
#' @param saxAlphabetSize The alphabet size when converting data into SAX representation
#' @param delay The number of zeros to pad on the hidden sequence before converting to SAX in order to delay it
#' @param dMarkovSize The size when grouping observed data to predict the next hidden state
#' @param disStep The step size of sampling when computing the distribution
#' @param disWindowSize The size of the sample when computing the distribution
#' @param k The maximum number of nearest neighbours when computing the KL divergence
#' @return This function does not return anything. However it does plot a graph of the KLD
#' @export
runKLDVecVec <- function(rawData1, rawData2, saxGroupSize, saxAlphabetSize, delay, dMarkovSize, disStep, disWindowSize,k){
  dataSAX1 = rspatiotemp::convert(rspatiotemp::runSAX(rawData1, saxGroupSize, saxAlphabetSize),dMarkovSize,saxAlphabetSize)
  dataSAXD1 = rspatiotemp::runSAX(c(rep(0,delay),rawData1), saxGroupSize, saxAlphabetSize)
  dataSAX2 = rspatiotemp::convert(rspatiotemp::runSAX(rawData2, saxGroupSize, saxAlphabetSize),dMarkovSize,saxAlphabetSize)
  dataSAXD2 = rspatiotemp::runSAX(c(rep(0,delay),rawData2), saxGroupSize, saxAlphabetSize)

  countData = rspatiotemp::count(dataSAX1,dataSAXD1,saxAlphabetSize^dMarkovSize,saxAlphabetSize)
  dataDis1 = rspatiotemp::distribution(countData,dataSAX1,dataSAXD1,saxAlphabetSize^dMarkovSize,saxAlphabetSize,disStep,disWindowSize)
  dataDis2 = rspatiotemp::distribution(countData,dataSAX2,dataSAXD2,saxAlphabetSize^dMarkovSize,saxAlphabetSize,disStep,disWindowSize)

  plot(FNN::KL.divergence(10^dataDis1,10^dataDis2,k=k))
}

#' Convert time series data to SAX, randomize test sequence, compute each distribution then graph the KL divergence between the two sequences
#' @title KL Divergence of a single time series and a random data set (runKLDVecRand)
#' @param rawData A vector of the raw data/time series data
#' @param saxGroupSize The group size when converting data into SAX representation
#' @param saxAlphabetSize The alphabet size when converting data into SAX representation
#' @param delay The number of zeros to pad on the hidden sequence before converting to SAX in order to delay it
#' @param dMarkovSize The size when grouping observed data to predict the next hidden state
#' @param disStep The step size of sampling when computing the distribution
#' @param disWindowSize The size of the sample when computing the distribution
#' @param k The maximum number of nearest neighbours when computing the KL divergence
#' @return This function does not return anything. However it does plot a graph of the KLD
#' @export
runKLDVecRand <- function(rawData, saxGroupSize, saxAlphabetSize, delay, dMarkovSize, disStep, disWindowSize,k){
  dataSAX1 = rspatiotemp::convert(rspatiotemp::runSAX(rawData, saxGroupSize, saxAlphabetSize),dMarkovSize,saxAlphabetSize)
  dataSAXD1 = rspatiotemp::runSAX(c(rep(0,delay),rawData), saxGroupSize, saxAlphabetSize)

  countData = rspatiotemp::count(dataSAX1,dataSAXD1,saxAlphabetSize^dMarkovSize,saxAlphabetSize)
  dataDis1 = rspatiotemp::distribution(countData,dataSAX1,dataSAXD1,saxAlphabetSize^dMarkovSize,saxAlphabetSize,disStep,disWindowSize)
  max = saxAlphabetSize-1
  dataDis2 = rspatiotemp::distribution(countData,dataSAX1,sample(0:max,length(dataSAXD1),T),saxAlphabetSize^dMarkovSize,saxAlphabetSize,disStep,disWindowSize)

  plot(FNN::KL.divergence(10^dataDis1,10^dataDis2,k=k))
}

#' Convert time series data to SAX, group Train and Test sequences, compute each distribution then graph the KL divergence between every combination of two sequences
#' @title KL Divergence of multiple time series (runKLDVec)
#' @param rawData A vector of the raw data/time series data
#' @param saxGroupSize The group size when converting data into SAX representation
#' @param saxAlphabetSize The alphabet size when converting data into SAX representation
#' @param dMarkovSize The size when grouping observed data to predict the next hidden state
#' @param disStep The step size of sampling when computing the distribution
#' @param disWindowSize The size of the sample when computing the distribution
#' @param k The maximum number of nearest neighbours when computing the KL divergence
#' @return This function does not return anything. However it does create a png of all the KL divergences
#' @export
runKLDMtx <- function(rawData, saxGroupSize, saxAlphabetSize, dMarkovSize, disStep, disWindowSize, k){
  tempSAX = rspatiotemp::runSAX(rawData[,1], saxGroupSize, saxAlphabetSize)
  dataSAXOO = .Call('rspatiotemp_convert', PACKAGE = 'rspatiotemp', tempSAX, dMarkovSize, saxAlphabetSize)
  dataSAXHO = tempSAX[dMarkovSize:length(tempSAX)]
  dataSAXO = matrix(dataSAXOO,length(dataSAXOO),1)
  dataSAXH = matrix(dataSAXHO,length(dataSAXHO),1)
  for(i in 2:ncol(rawData)){
    tempSAX = rspatiotemp::runSAX(rawData[,i], saxGroupSize, saxAlphabetSize)
    tempSAXO = .Call('rspatiotemp_convert', PACKAGE = 'rspatiotemp', tempSAX, dMarkovSize, saxAlphabetSize)
    tempSAXH = tempSAX[dMarkovSize:length(tempSAX)]
    dataSAXO = cbind(dataSAXO, tempSAXO)
    dataSAXH = cbind(dataSAXH,tempSAXH)
  }
  halfway = nrow(dataSAXH)/2
  colNames = colnames(rawData)
  for(outer in 96:ncol(rawData)){
    dataOTr = dataSAXO[,outer]
    dataHTr = dataSAXH[,outer]
    dataCount = .Call('rspatiotemp_count', PACKAGE = 'rspatiotemp', dataOTr,dataHTr,saxAlphabetSize^dMarkovSize,saxAlphabetSize)
    disTrTr = .Call('rspatiotemp_distribution', PACKAGE = 'rspatiotemp', dataCount, dataOTr, dataHTr, saxAlphabetSize^dMarkovSize, saxAlphabetSize, disStep, disWindowSize)
    for(inner in 1:ncol(rawData)){
      if(outer!=inner){
        dataOTe = dataSAXO[,inner]
        dataHTe = dataSAXH[,inner]
      }
      else{
        dataOTe = dataOTr[halfway:length(dataHTr)]
        dataHTe = dataHTr[halfway:length(dataHTr)]
        dataOTr = dataOTr[0:halfway]
        dataHTr = dataHTr[0:halfway]
        dataCount = .Call('rspatiotemp_count', PACKAGE = 'rspatiotemp', dataOTr,dataHTr,saxAlphabetSize^dMarkovSize,saxAlphabetSize)
        disTrTr = .Call('rspatiotemp_distribution', PACKAGE = 'rspatiotemp', dataCount, dataOTr, dataHTr, saxAlphabetSize^dMarkovSize, saxAlphabetSize, disStep, disWindowSize)
      }
      disTrTe = .Call('rspatiotemp_distribution', PACKAGE = 'rspatiotemp', dataCount, dataOTe, dataHTe, saxAlphabetSize^dMarkovSize, saxAlphabetSize, disStep, disWindowSize)
      filename<- paste("Trian v",as.character(outer), " - ","Test v",as.character(inner), " unlog.png", sep="")
      png(file=filename)
      #plot(FNN::KL.divergence(10^disTrTr,10^disTrTe,k=k))
      len = length(disTrTe)-1
      cont = 0
      tryCatch({
        if(max(10^disTrTe[0:len])>=0.1){
          hist(10^disTrTe)
          cont = 2
        }
      },
      error = function(e){ print(filename) })
      dev.off()
      if(cont>1){
        filename<- paste("Dis v",as.character(outer), " - ","v",as.character(inner), ".png", sep="")
        png(file=filename)
        tryCatch({ plot(FNN::KL.divergence(10^disTrTr,10^disTrTe,k=k)) },
                 error = function(e){ print(filename) })
        dev.off()
      }
    }
  }
}

repeatKLDSim <- function(n){
  hbSAX = convert(rspatiotemp::runSAX(waveletdecomp::heartbeat, 10, 20),2,20)
  hbSAXD = runSAX(c(rep(0,7),waveletdecomp::heartbeat), 10, 20)

  obsLen = length(hbSAX)
  hbProbMatX = createProbMatX(hbSAX,hbSAXD[0:obsLen],2,400,20)


  countHB = count(hbSAX,hbSAXD,400,20)
  hbDis = distribution(countHB,hbSAX,hbSAXD,400,20,1,100)

  for(i in 0:n){
    hbSim = simulateHid(hbProbMatX,hbSAX,2,400,20)
    hbSimDis = distribution(countHB,hbSAX,hbSim,400,20,1,100)
    kld = FNN::KL.divergence(10^hbDis,10^hbSimDis,1000)
    filename<- paste(0,".png", sep="")
    png(file=filename)
    plot(kld)
    dev.off()
  }
}

