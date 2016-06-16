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

adf.testWindow <- function(data, stepSize, windowSize){
  indexSeq = seq(from = 0, to = length(data)-windowSize, by = stepSize)
  result = numeric()
  for(i in indexSeq){
    resultTemp = tseries::adf.test(data[i:(i+windowSize)])[4]
    result = c(result,as.numeric(resultTemp))
  }
  return(result)
}
