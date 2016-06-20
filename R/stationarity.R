#' Computes the stationarityy using forecast::ndiffs(test = "adf") of several windows of the original data
#' @title forecast ndiffs Window (ndiffsWindow)
#' @param data The time series data
#' @param stepSize The about of elements to jump after each window
#' @param windowSize The size of each window being sampled
#' @return A vector of the computed values for each window returned from the ndiffs function
#' @export
ndiffsWindow <- function(data, stepSize, windowSize){
  indexSeq = seq(from = 0, to = length(data)-windowSize, by = stepSize)
  result = integer()
  for(i in indexSeq){
    resultTemp = forecast::ndiffs(data[i:(i+windowSize)],test = "adf")
    result = c(result,resultTemp)
  }
  return(result)
}

#' Computes the stationarityy using tseries::adf.test of several windows of the original data
#' @title tseries adf.test Window (adf.testWindow)
#' @param data The time series data
#' @param stepSize The about of elements to jump after each window
#' @param windowSize The size of each window being sampled
#' @return A vector of the computed values for each window returned from the adf.test function
#' @export
adf.testWindow <- function(data, stepSize, windowSize){
  indexSeq = seq(from = 0, to = length(data)-windowSize, by = stepSize)
  result = numeric()
  for(i in indexSeq){
    resultTemp = tseries::adf.test(data[i:(i+windowSize)])[4]
    result = c(result,as.numeric(resultTemp))
  }
  return(result)
}
