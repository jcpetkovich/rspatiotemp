detectOutlierPositionFft <- function(signal, thresholdFreq = 0.1, frequencyAmp = 0.01){
  fftSignal = stats::fft(signal)
  if(abs(max(signal)) > abs(min(signal)))
    outlier = max(signal)
  else
    outlier = min(signal)

  if(any(abs(fftSignal) > frequencyAmp)){
    outlierIndex = match(outlier,signal)
    return(outlierIndex)
  }
  else
    return(NaN)
}
#' Detecting outliers in a series using FFT(Fast Fourier Transform)
#' @title Detect outliers FFT (detectOutliersFft)
#' @param data A vector containing the time series data
#' @param stepSize The size of half a window that is sampled from the data
#' @param thresholdFreq Threshold frequency compared with frequency response of the signal to determine outliers
#' @param frequencyAmp Frequency Amplitude
#' @return A vector with the indices of the detected outliers
#' @export
detectOutliersFft <- function(data, stepSize, thresholdFreq = 0.1, frequencyAmp = 0.01){
  outliers = integer()
  for(i in seq(stepSize, length(data)-stepSize,stepSize)){
    lB = i-stepSize
    uB = i+stepSize
    outlierPos = detectOutlierPositionFft(data[lB:uB],thresholdFreq,frequencyAmp)

    if(!is.nan(outlierPos)){
      outlierPos = outlierPos + i - stepSize
      outliers = c(outliers,outlierPos)
    }
  }
  return(outliers)
}
