#' @export
readToVec <- function(data.path){
  data = read.csv(data.path)
  data$X = NULL
  data = unlist(data)
  data = as.numeric(data)
  return(data)
}
