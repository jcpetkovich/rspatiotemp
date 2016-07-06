#' @export
lifeTab.allFile.vec <-function(data.path,timeInterval,colH,colV){
  all.files = file.path(data.path,list.files(data.path))
  fileGrpSize = length(read.table(all.files[1])[,1])
  fileNum = length(all.files)
  allData1 = lapply(all.files,function(file){read.table(file)[,colH]})
  allData2 = lapply(all.files,function(file){read.table(file)[,colV]})
  allData1 = unlist(allData1)
  allData2 = unlist(allData2)
  allData1 = matrix(allData1,fileGrpSize,fileNum)
  allData2 = matrix(allData2,fileGrpSize,fileNum)
  accDeg = rspatiotemp::accDegrad(allData1,allData2,1,1,1,0.5,"exp","poly")
  lifeTable = rspatiotemp::createLifeTab(accDeg,timeInterval)
  return(lifeTable)
}

#' @export
RUL.vec <- function(lifeTable, data.path){
  all.files = file.path(data.path,list.files(data.path))
  fileGrpSize = length(read.table(all.files[1])[,1])
  fileNum = length(all.files)
  allData1 = lapply(all.files,function(file){read.table(file)[,7]})
  allData2 = lapply(all.files,function(file){read.table(file)[,8]})
  allData1 = unlist(allData1)
  allData2 = unlist(allData2)
  allData1 = matrix(allData1,fileGrpSize,fileNum)
  allData2 = matrix(allData2,fileGrpSize,fileNum)
  accDeg = rspatiotemp::accDegrad(allData1,allData2,2,2,2,0.5,"exp","poly")
  time = numeric()
  for(i in 1:length(accDeg)){
    timeRUL = rspatiotemp::computeRUL(lifeTable,accDeg[i])
    time = c(time,timeRUL)
  }
  return(cbind(accDeg,time))
}
#' @export
RUL.mat <- function(lifeTable,accDeg){
  time = numeric()
  for(i in 1:length(accDeg)){
    timeRUL = rspatiotemp::computeRUL(lifeTable,accDeg[i])
    time = c(time,timeRUL)
  }
  return(cbind(accDeg,time))
}

#' @export
alldata <- function(data.path,col){
  all.files = file.path(data.path,list.files(data.path))
  fileGrpSize = length(read.table(all.files[1])[,1])
  fileNum = length(all.files)
  allData = lapply(all.files,function(file){read.table(file)[,col]})
  allData = unlist(allData)
  allData = matrix(allData,fileGrpSize,fileNum)
  return(allData)
}

#' @export
lifeTab <- function(allData1,allData2,timeInterval){
  accDeg = rspatiotemp::accDegrad(allData1,allData2,1,1,1,0.5,"exp","poly")
  lifeTable = rspatiotemp::createLifeTab(accDeg,timeInterval)
  return(lifeTable)
}

#' @export
graphRUL <- function(data.path){
  all.files = file.path(data.path,list.files(data.path))
  fileGrpSize = length(read.csv(all.files[1])[,1])
  fileNum = length(all.files)
  allRULdata = lapply(all.files,function(file){read.csv(file)[,3]})
  allRULdata = unlist(allRULdata)
  allRULdata = matrix(allRULdata, fileGrpSize,fileNum)
  return(allRULdata)
}
