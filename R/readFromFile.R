#' @export
readToVec <- function(data.path){
  data = read.csv(data.path)
  data$X = NULL
  data = unlist(data)
  data = as.numeric(data)
  return(data)
}

#' @export
readToVecPart <- function(data.path, maxSize){
  data = read.csv(data.path)
  data$X = NULL
  data = unlist(data)
  data = as.numeric(data)
  return(data[1:maxSize])
}

convertData <- function(file.path){
  all.files = file.path(file.path,list.files(file.path))
  len = length(all.files)
  for(i in 1:len){
    data = read.csv(all.files[i])
  }
}

readFEMTO <- function(data.path,save1){
  all.files = file.path(data.path,list.files(data.path))
  data = read.csv(all.files[1])[5:6]
  horizontal = matrix(data[1],ncol=1)
  vertical = matrix(data[2],ncol=1)
  for(i in all.files[-1]){
    data = read.csv(i)[5:6]
    horizontal = cbind(horizontal,data[1])
    vertical = cbind(vertical,data[2])
  }
  save(horizontal,file = paste("bearingDataH",save1,".Rd",sep=""))
  save(vertical,file = paste("bearingDataV",save1,".Rd",sep=""))
}
