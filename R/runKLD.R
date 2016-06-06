#' Convert time series data to SAX, group Train and Test sequences, compute each distribution then graph the KL divergence between the two sequences
#' @title KL Divergence of a single time series (runKLDVec)
#' @param rawData A vector of the raw data/time series data
#' @param saxGroupSize The group size when converting data into SAX representation
#' @param saxAlphabetSize The alphabet size when converting data into SAX representation
#' @param dMarkovSize The size when grouping observed data to predict the next hidden state
#' @param disStep The step size of sampling when computing the distribution
#' @param disWindowSize The size of the sample when computing the distribution
#' @param k The maximum number of nearest neighbours when computing the KL divergence
#' @return This function does not return anything. However it does create a png of the KL divergence
#' @export
runKLDVec <- function(rawData, saxGroupSize, saxAlphabetSize, dMarkovSize, disStep, disWindowSize, k){
  dataSAX = .Call('rspatiotemp_runSAX', PACKAGE = 'rspatiotemp', rawData, saxGroupSize, saxAlphabetSize,F)
  dataO = .Call('rspatiotemp_convert', PACKAGE = 'rspatiotemp', dataSAX, dMarkovSize, saxAlphabetSize)
  dataH = dataSAX[dMarkovSize:length(dataSAX)]
  halfway = length(dataH)/2
  dataOTr = dataO[0:halfway]
  dataOTe = dataO[halfway:length(dataH)]
  dataHTr = dataH[0:halfway]
  dataHTe = dataH[halfway:length(dataH)]
  dataCount = .Call('rspatiotemp_count', PACKAGE = 'rspatiotemp', dataOTr,dataHTr,saxAlphabetSize^dMarkovSize,saxAlphabetSize)
  disTrTr = .Call('rspatiotemp_distribution', PACKAGE = 'rspatiotemp', dataCount, dataOTr, dataHTr, saxAlphabetSize^dMarkovSize, saxAlphabetSize, disStep, disWindowSize)
  disTrTe = .Call('rspatiotemp_distribution', PACKAGE = 'rspatiotemp', dataCount, dataOTe, dataHTe, saxAlphabetSize^dMarkovSize, saxAlphabetSize, disStep, disWindowSize)
  png('rplot.png')
  plot(FNN::KL.divergence(10^disTrTr,10^disTrTe,k=k))
  dev.off()
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
  for(outer in 1:ncol(rawData)){
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
      tryCatch({ hist(10^disTrTe) },
               error = function(e){ print(filename) })
      dev.off()
      filename<- paste("Trian v",as.character(outer), " - ","Test v",as.character(inner), " log.png", sep="")
      png(file=filename)
      tryCatch({ hist(disTrTe) },
               error = function(e){ print(filename) })
      dev.off()
    }
  }
}

