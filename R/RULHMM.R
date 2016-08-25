#library(wavethresh)
#library(dplyr)

##' @import parallel

#' Extract features using Wavelet Packet Decomposition, compute nodal energies then convert to SAX.
#' @title To WPD, to SAX (toWPDSAX)
#' @param timeSeries The time series to be converted.
#' @param SAXalphabetSize The alphabet size used for the SAX conversion.
#' @param SAXgroupSize The size of each group to be converted to a single SAX symbol.
#' @param exp2 Two to the power of exp2 will be the size of each group fed into the WPD of each group fed into the WPD.
#' @return The SAX sequence of the given time series.
#' @export
toWPDSAX <- function(timeSeries, SAXalphabetSize, SAXgroupSize, exp2){
  by = 2^exp2
  splitIndex = seq(from=0,to=length(timeSeries),by=by)
  splitLen = length(splitIndex)-1
  wpAll = matrix(0,by,splitLen)
  for(i in 1:splitLen){
    lowBound = splitIndex[i]+1
    upBound = splitIndex[i+1]
    pieceWP = wavethresh::wp(timeSeries[lowBound:upBound])$wp
    wpAll[,i] = getEnergy(pieceWP)
  }
  wpAll = as.numeric(wpAll)
  wpSAX = rspatiotemp::runSAX(wpAll,SAXgroupSize,SAXalphabetSize)
  return(wpSAX)
}

toWPDSAX.path <- function(data.path, SAXalphabetSize, SAXgroupSize, exp2){
  timeSeries = rspatiotemp::readToVec(data.path)
  return(toWPDSAX(timeSeries,SAXalphabetSize,SAXgroupSize,exp2))
}

#' Compute nodal energies of the wavelet packet decomposition.
#' @title Nodal Energy of WPD (wpdEnergy)
#' @param timeSeries The time series data fed into the WPD.
#' @param exp2 Two to the power of exp2 will be the size of each group fed into the WPD.
#' @return A vector containing the nodal energy values.
#' @export
wpdEnergy <- function(timeSeries, exp2){
  by = 2^exp2
  splitIndex = seq(from=0,to=length(timeSeries),by=by)
  splitLen = length(splitIndex)-1
  wpAll = matrix(0,by,splitLen)
  for(i in 1:splitLen){
    lowBound = splitIndex[i]+1
    upBound = splitIndex[i+1]
    pieceWP = wavethresh::wp(timeSeries[lowBound:upBound])$wp
    wpAll[,i] = getEnergy(pieceWP)
  }
  wpAll = as.numeric(wpAll)
  return(wpAll)
}

#' Train Hidden Markov Model using the depmixS4 R package.
#' @title Train HMM depmixS4 (trainHMM.dm)
#' @param energySeries The nodal energies of the time series, output of 'wpdEnergy' function.
#' @param nstates The alphabet size for the hidden state.
#' @return The probability matrices of the hmm and the viterbi sequence.
#' @export
trainHMM.dm <- function(energySeries, nstates){
  mix = depmixS4::depmix(energy~1, data = energySeries, nstates = nstates)
  fit = depmixS4::fit(mix)
  trans = matrix(as.numeric(fit@trDens),nstates,nstates,byrow = T)
  init = fit@init
  emis = numeric()
  for(i in 1:nstates){
    param = fit@response[[i]]
    emis = c(emis,param[[1]]@parameters$coefficients,param[[1]]@parameters$sd)
  }
  emis = matrix(emis,ncol = nstates)
  vitSeq = fit@posterior$state
  rownames(emis) = c("mean","sd")
  hmm = list(Transition = trans,Initial = init, Emission = emis, vitSeq = vitSeq)
  return(hmm)
}

#' Train Hidden Markov Models using the Baum Welch Algorithm
#' @title Train Hidden Markov Model (trainHMM)
#' @param ObsSeq A discrete observed time series. Output from toWPDSAX.
#' @param obsAlphabetSize The alphabet size of the observed sequence.
#' @param hidAlphabetSize The alphabet size of the hidden sequence.
#' @return The probability matrices, mean, standard deviation vectors and the final failure state of the hidden markov model.
#' @export
trainHMM <- function(ObsSeq, obsAlphabetSize, hidAlphabetSize){
  return(rspatiotemp::train(matrix(1/hidAlphabetSize,hidAlphabetSize,hidAlphabetSize),matrix(1/obsAlphabetSize,hidAlphabetSize,obsAlphabetSize),rep(1/hidAlphabetSize,hidAlphabetSize),ObsSeq))
}

#' Compute the most probable hidden sequence given an observed sequence and an HMM.
#' @title Viterbi Algorithm (viterbiR)
#' @param ObsSeq A discrete observed time series. Output from toWPDSAX.
#' @param probMat The probability matrices of the hidden markov model. Output from trainHMM.
#' @return The most probably hidden sequence.
#' @export
viterbiR <- function(ObsSeq, probMat){
  return(rspatiotemp::viterbi(probMat$Transition,probMat$Emission,probMat$Initial,ObsSeq))
}

#' Compute mean and standard deviation of the time in each state
#' @title Compute mean and standard deviation (meanStdDevR)
#' @param hidSeq The hidden sequence. Output of viterbiR.
#' @param hidAlphabetSize The alphabet size of the hidden sequence
#' @return Vectors of the mean and median values.
#' @export
meanStdDevR <- function(hidSeq, hidAlphabetSize,tabB = FALSE){
  tab = rspatiotemp::formatRULHMM(hidSeq)
  return(meanStdDev(tab$redSeq,tab$repSeq,hidAlphabetSize,tabB))
}

#' Create a single hidden markov model from a given time series.
#' @title Create HMM (createModel)
#' @param timeSeries The time series to be converted.
#' @param SAXalphabetSize The alphabet size used for the SAX conversion.
#' @param SAXgroupSize The size of each group to be converted to a single SAX symbol.
#' @param exp2 Two to the power of exp2 will be the size of each group fed into the WPD.
#' @param hidAlphabetSize The alphabet size of the hidden sequence.
#' @return A HMM generated from the time series.
#' @export
createModel <- function(timeSeries, SAXalphabetSize, SAXgroupSize, exp2, hidAlphabetSize, tab){
  obsSAX = toWPDSAX(timeSeries,SAXalphabetSize, SAXgroupSize,exp2)
  probMat = trainHMM(obsSAX, SAXalphabetSize,hidAlphabetSize)
  vitSeq = viterbiR(obsSAX,probMat)
  MaSD = meanStdDevR(vitSeq,hidAlphabetSize,tab)
  if(tab)
    return(list(Transition = probMat$Transition, Emission = probMat$Emission, Initial = probMat$Initial, Mean  = MaSD$mean, StdDev = MaSD$stdDev,Tab = MaSD$tab, Final = vitSeq[length(vitSeq)]))
  return(list(Transition = probMat$Transition, Emission = probMat$Emission, Initial = probMat$Initial, Mean  = MaSD$mean, StdDev = MaSD$stdDev, Final = vitSeq[length(vitSeq)]))
}

#' Create a single hidden markov model from the given time series located at the file path.
#' @title Create HMM (createModel.path)
#' @param data.path The file path to the data.
#' @param SAXalphabetSize The alphabet size used for the SAX conversion.
#' @param SAXgroupSize The size of each group to be converted to a single SAX symbol.
#' @param exp2 Two to the power of exp2 will be the size of each group fed into the WPD.
#' @param hidAlphabetSize The alphabet size of the hidden sequence.
#' @return A HMM generated from the time series.
#' @export
createModel.path <- function(data.path,SAXalphabetSize,SAXgroupSize,exp2,hidAlphabetSize, tab){
  timeSeries = rspatiotemp::readToVec(data.path)
  return(createModel(timeSeries,SAXalphabetSize,SAXgroupSize,exp2,hidAlphabetSize,tab))
}

#' Create a several hidden markov model for the sets of data located in the given file path.
#' @title Create HMM (createModels.path)
#' @param data.path The file path to the sets data.
#' @param SAXalphabetSize The alphabet size used for the SAX conversion.
#' @param SAXgroupSize The size of each group to be converted to a single SAX symbol.
#' @param exp2 Two to the power of exp2 will be the size of each group fed into the WPD.
#' @param hidAlphabetSize The alphabet size of the hidden sequence.
#' @param tab True or False, using the modified tab version of the RUL HMM algorithm.
#' @return A list of HMM generated from the time series.
#' @export
createModels.path <- function(data.path,SAXalphabetSize,SAXgroupSize,exp2,hidAlphabetSize, tab){
  all.files = file.path(data.path,list.files(data.path))
  models = lapply(all.files, function(file){createModel.path(file,SAXalphabetSize,SAXgroupSize,exp2,hidAlphabetSize,tab)})
  return(models);
}

#' Create a single hidden markov model from a given time series using the depmixS4 package.
#' @title Create HMM (createModel.dm)
#' @param timeSeries The time series to be converted.
#' @param exp2 Two to the power of exp2 will be the size of each group fed into the WPD.
#' @param hidAlphabetSize The alphabet size of the hidden sequence.
#' @return A HMM generated from the time series.
#' @export
createModel.dm <- function(timeSeries, exp2, hidAlphabetSize){
  wpd = wpdEnergy(timeSeries,exp2)
  wpd = data.frame(energy = wpd)
  hmm = trainHMM.dm(wpd,hidAlphabetSize)
  MaSD = meanStdDevR(hmm$vitSeq,hidAlphabetSize)
  return(list(Transition = hmm$Transition,Emission = hmm$Emission,Initial = hmm$Initial, Mean  = MaSD$mean, StdDev = MaSD$stdDev, Final = hmm$vitSeq[length(hmm$vitSeq)]))
}

#' Create a single hidden markov model from the given time series located at the file path using the depmixS4 package.
#' @title Create HMM (createModel.path.dm)
#' @param data.path The file path to the data.
#' @param exp2 Two to the power of exp2 will be the size of each group fed into the WPD.
#' @param hidAlphabetSize The alphabet size of the hidden sequence.
#' @return A HMM generated from the time series.
#' @export
createModel.path.dm <- function(data.path, exp2, hidAlphabetSize){
  timeSeries = rspatiotemp::readToVec(data.path)
  return(createModel.dm(timeSeries,exp2,hidAlphabetSize))
}

#' Create a several hidden markov model for the sets of data located in the given file path using the depmixS4 package.
#' @title Create HMM (createModels.path.dm)
#' @param data.path The file path to the sets data.
#' @param exp2 Two to the power of exp2 will be the size of each group fed into the WPD.
#' @param hidAlphabetSize The alphabet size of the hidden sequence.
#' @return A list of HMM generated from the time series.
#' @export
createModels.path.dm <- function(data.path, exp2, hidAlphabetSize){
  all.files = file.path(data.path,list.files(data.path))
  models = mclapply(all.files, function(file){createModel.path.dm(file,exp2,hidAlphabetSize)}, mc.cores = detectCores())
  return(models)
}

#' Create a single hidden markov model from a given time series using the depmixS4 package with a tab of all states.
#' @title Create HMM (createModel.dm.tab)
#' @param timeSeries The time series to be converted.
#' @param exp2 Two to the power of exp2 will be the size of each group fed into the WPD.
#' @param hidAlphabetSize The alphabet size of the hidden sequence.
#' @return A HMM generated from the time series.
#' @export
createModel.dm.tab <- function(timeSeries, exp2, hidAlphabetSize){
  #wpd = wpdEnergy(timeSeries,exp2)
  wpd = getEnergies.wpd(timeSeries,exp2,3)
  wpd = data.frame(energy = wpd[,1])
  hmm = trainHMM.dm(wpd,hidAlphabetSize)
  MaSD = meanStdDevR(hmm$vitSeq,hidAlphabetSize,T)
  return(list(Transition = hmm$Transition,Emission = hmm$Emission,Initial = hmm$Initial, Mean  = MaSD$mean, StdDev = MaSD$stdDev,Tab = MaSD$tab, Final = hmm$vitSeq[length(hmm$vitSeq)]))
}

#' Create a single hidden markov model from the given time series located at the file path using the depmixS4 package with a tab of all states.
#' @title Create HMM (createModel.path.dm.tab)
#' @param data.path The file path to the data.
#' @param exp2 Two to the power of exp2 will be the size of each group fed into the WPD.
#' @param hidAlphabetSize The alphabet size of the hidden sequence.
#' @return A HMM generated from the time series.
#' @export
createModel.path.dm.tab <- function(data.path, exp2, hidAlphabetSize){
  #timeSeries = rspatiotemp::readToVec(data.path)
  load(data.path)
  data = unlist(data)
  data = as.numeric(data)
  return(createModel.dm.tab(data,exp2,hidAlphabetSize))
}

#' Create a single hidden markov model from the given time series located at the file path using the depmixS4 package with a tab of all states.
#' @title Create HMM (createModel.path.dm.tab)
#' @param data.path The file path to the data.
#' @param exp2 Two to the power of exp2 will be the size of each group fed into the WPD.
#' @param hidAlphabetSize The alphabet size of the hidden sequence.
#' @param maxSize Only generate the model using that amount of data from the file.
#' @return A HMM generated from the time series.
#' @export
createPartModel.path.dm.tab <- function(data.path, exp2, hidAlphabetSize,maxSize){
  #timeSeries = rspatiotemp::readToVecPart(data.path,maxSize)
  load(data.path)
  data = unlist(data)
  data = as.numeric(data)
  return(createModel.dm.tab(data[1:maxSize],exp2,hidAlphabetSize))
}

#' Create a several hidden markov model for the sets of data located in the given file path using the depmixS4 package with a tab of all states.
#' @title Create HMM (createModels.path.dm.tab)
#' @param data.path The file path to the sets data.
#' @param exp2 Two to the power of exp2 will be the size of each group fed into the WPD.
#' @param hidAlphabetSize The alphabet size of the hidden sequence.
#' @return A list of HMM generated from the time series.
#' @export
createModels.path.dm.tab <- function(data.path, exp2, hidAlphabetSize){
  all.files = file.path(data.path,list.files(data.path))
  models = lapply(all.files, function(file){createModel.path.dm.tab(file,exp2,hidAlphabetSize)})
  return(models)
}

#' Create a several hidden markov model for the sets of data located in the given file path using the depmixS4 package with a tab of all states.
#' @title Create HMM (createModels.path.dm.tab)
#' @param data.path The file path to the sets data.
#' @param exp2 Two to the power of exp2 will be the size of each group fed into the WPD.
#' @param hidAlphabetSize The alphabet size of the hidden sequence.
#' @param maxSize Only generate the model using that amount of data from the file.
#' @return A list of HMM generated from the time series.
#' @export
createPartModels.path.dm.tab <- function(data.path, exp2, hidAlphabetSize,maxSize){
  all.files = file.path(data.path,list.files(data.path))
  models = lapply(all.files, function(file){createModel.path.dm.tab(file,exp2,hidAlphabetSize)})
  return(models)
}
#' Compute pessimistic remaining useful life using HMMs
#' @title Remaining useful life using HMM (RUL.HMM)
#' @param timeSeries The time series to be converted.
#' @param models The generated hidden markov models. Output from createModel(s).*
#' @param SAXalphabetSize The alphabet size used for the SAX conversion.
#' @param SAXgroupSize The size of each group to be converted to a single SAX symbol.
#' @param exp2 Two to the power of exp2 will be the size of each group fed into the WPD.
#' @param lastStateScanNum The amount of states scanned when finding the last sequence state.
#' @param confidenceCoef A constant used when computing the lower and upper bounds of the RUL.
#' @param tab True or False, using the modified tab version of the RUL HMM algorithm.
#' @return The estimated pessimistic remaining useful life with an upper and lower bound.
#' @export
RUL.HMM <- function(timeSeries, models, SAXalphabetSize, SAXgroupSize, exp2, lastStateScanNum, confidenceCoef, tab){
  #Most probable model
  testSAX = toWPDSAX(timeSeries,SAXalphabetSize,SAXgroupSize,exp2)
  modLen = length(models)
  probScore = numeric()
  for(i in 1:modLen){
    score = rspatiotemp::probObsLog(models[[i]]$Transition,models[[i]]$Emission,testSAX)
    probScore = c(probScore,score)
  }
  chosenModelIndex = which.max(probScore)
  #last state
  vitSeq = viterbiR(testSAX,models[[chosenModelIndex]])
  if(tab){
    testModel = createModel(timeSeries,SAXalphabetSize,SAXgroupSize,exp2,nrow(models[[1]]$Transition),T)
    tab = models[[chosenModelIndex]]$Tab - testModel$Tab
    tab[tab < 0] = 0
    mid = sum(tab*models[[chosenModelIndex]]$Mean)
    error = sum(tab*confidenceCoef*models[[chosenModelIndex]]$StdDev)
    return(list(Lower = mid-error, Middle = mid, Upper = mid+error))
  }
  lowerBound = length(vitSeq)-lastStateScanNum
  upperBound = length(vitSeq)
  lastState = mode(vitSeq[lowerBound:upperBound])
  #Critical Path
  possibleStates = nrow(models[[chosenModelIndex]]$Transition)-1
  path = rspatiotemp::criticalPath(models[[chosenModelIndex]]$Transition,lastState, models[[chosenModelIndex]]$Final,0:possibleStates)
  return(computeRULBounds(path,models[[chosenModelIndex]]$Mean,models[[chosenModelIndex]]$StdDev,confidenceCoef))
}

#' Compute pessimistic remaining useful life using HMMs using the package depmixS4
#' @title Remaining useful life using HMM (RUL.HMM.dm)
#' @param timeSeries The time series to be converted.
#' @param models The generated hidden markov models. Output from createModel(s).dm.*
#' @param exp2 Two to the power of exp2 will be the size of each group fed into the WPD.
#' @param lastStateScanNum The amount of states scanned when finding the last sequence state.
#' @param confidenceCoef A constant used when computing the lower and upper bounds of the RUL.
#' @return The estimated pessimistic remaining useful life with an upper and lower bound.
#' @export
RUL.HMM.dm <- function(timeSeries, models, exp2, lastStateScanNum, confidenceCoef){
  hidAlphabetSize = nrow(models[[1]]$Transition)
  wpdNum = wpdEnergy(timeSeries,exp2)
  wpd = data.frame(energy = wpdNum)
  hmm = trainHMM.dm(wpd,hidAlphabetSize)
  len = length(models)
  probScore = numeric()
  for(i in 1:len){
    probScore = c(probScore,probObsLog(models[[i]]$Transition,models[[i]]$Emission,wpdNum))
  }
  chosenModelIndex = which.max(probScore)
  len = length(hmm$vitSeq)
  lower = len - lastStateScanNum
  lastState = mode(hmm$vitSeq[lower:len])
  possibleStates = hidAlphabetSize-1
  path = rspatiotemp::criticalPath(models[[chosenModelIndex]]$Transition,lastState, models[[chosenModelIndex]]$Final,0:possibleStates)
  return(computeRULBounds(path,models[[chosenModelIndex]]$Mean,models[[chosenModelIndex]]$StdDev,confidenceCoef))
}

#' Compute the remaining useful life using HMMs using the package depmixS4 and tabs of states
#' @title Remaining useful life using HMM (RUL.HMM)
#' @param timeSeries The time series to be converted.
#' @param models The generated hidden markov models. Output from createModel(s).dm.tab*
#' @param exp2 Two to the power of exp2 will be the size of each group fed into the WPD.
#' @param confidenceCoef A constant used when computing the lower and upper bounds of the RUL.
#' @return The estimated remaining useful life with an upper and lower bound.
#' @export
RUL.dm.tab <- function(timeSeries, models, exp2, confidenceCoef){
  hidAlphabetSize = nrow(models[[1]]$Transition)
  #wpdNum = wpdEnergy(timeSeries,exp2)
  #wpd = data.frame(energy = wpdNum)
  wpdNum = getEnergies.wpd(timeSeries,exp2,3)
  wpd = data.frame(energy = wpdNum[,1])
  hmm = trainHMM.dm(wpd,hidAlphabetSize)
  MaSD = meanStdDevR(hmm$vitSeq,hidAlphabetSize,T)
  len = length(models)
  probScore = numeric()
  for(i in 1:len){
    probScore = c(probScore,probObsLog(models[[i]]$Transition,models[[i]]$Emission,wpdNum))
  }
  chosenModelIndex = which.max(probScore)
  tab = models[[chosenModelIndex]]$Tab - MaSD$tab
  tab[tab < 0] = 0
  mid = sum(tab*models[[chosenModelIndex]]$Mean)
  error = sum(tab*confidenceCoef*models[[chosenModelIndex]]$StdDev)
  return(list(Model = chosenModelIndex, Lower = mid-error, Middle = mid, Upper = mid+error))
}

# boundSeq = as.integer(seq(from = 0, to = 7175680, by = 2560))
# registerDoMC()
# RUL = foreach(i = 1:10,.combine = 'rbind',.inorder = TRUE) %dopar% {
#   unlist(RUL.dm.tab(h1_1[1:(717568*i)],list(modelH1_1_HMM,modelH1_1_HMM),6,0.1))
# }

#' Compute the remaining useful life using HMMs using the package depmixS4 and tabs of states
#' @title Remaining useful life using HMM (RUL.HMM)
#' @param timeSeries The time series to be converted.
#' @param trainData.path Path to the training data for partial model generation.
#' @param models The generated hidden markov models. Output from createModel(s).dm.tab*
#' @param exp2 Two to the power of exp2 will be the size of each group fed into the WPD.
#' @param confidenceCoef A constant used when computing the lower and upper bounds of the RUL.
#' @return The estimated remaining useful life with an upper and lower bound.
#' @export
RULpartial.dm.tab <- function(timeSeries, trainData.path, models, exp2, confidenceCoef){
  hidAlphabetSize = nrow(models[[1]]$Transition)
  wpdNum = wpdEnergy(timeSeries,exp2)
  wpd = data.frame(energy = wpdNum)
  hmm = trainHMM.dm(wpd,hidAlphabetSize)
  MaSD = meanStdDevR(hmm$vitSeq,hidAlphabetSize,T)
  partialModels = createPartModels.path.dm.tab(trainData.path,exp2,hidAlphabetSize,length(timeSeries))
  chosenModelIndex = pickPartModel.dm.tab(timeSeries, partialModels, exp2, hidAlphabetSize,wpdNum)
  tab = models[[chosenModelIndex]]$Tab - MaSD$tab
  tab[tab < 0] = 0
  mid = sum(tab*models[[chosenModelIndex]]$Mean)
  error = sum(tab*confidenceCoef*models[[chosenModelIndex]]$StdDev)
  return(list(Model = chosenModelIndex, Lower = mid-error, Middle = mid, Upper = mid+error))
}

#' Choosing best partial model using the forward algorithm.
#' @title Pick Partial Model Depmix Tab (pickPartModel.dm.tab)
#' @param timeSeries The test time series data to be matched with a partial model.
#' @param partialModels Generated partial models to be paired.
#' @param wpdNum The observed sequence converted to nodal energies.
#' @param The index of the paired model.
#' @export
pickPartModel.dm.tab <- function(timeSeries, partialModels,exp2, hidAlphabetSize, wpdNum){
  len = length(partialModels)
  probScore = numeric()
  for(i in 1:len){
    probScore = c(probScore,probObsLog(partialModels[[i]]$Transition,partialModels[[i]]$Emission,wpdNum))
  }
  chosenModelIndex = which.max(probScore)
  return(chosenModelIndex)
}
mode <- function(v){
  uniqv = unique(v)
  return(uniqv[which.max(tabulate(match(v, uniqv)))])
}

#' Compute Nodal Energies using Wavelet Packet Decomposition
#' @title Compute Nodal Energy WPD (getEnergies.wpd)
#' @param timeSeries The time series data.
#' @param exp2 Two to the power of exp2 will be the size of each group fed into the WPD of each group fed into the WPD.
#' @param levels The level of the WPD to compute the energy for.
#' @return A matrix of energies for each frequency band on the level chosen
#' @export
getEnergies.wpd <- function(timeSeries, exp2, levels){
  by = 2^exp2
  splitIndex = seq(from=0,to=length(timeSeries),by=by)
  splitLen = length(splitIndex)-1
  levelLen = (2^levels)-1
  registerDoMC()
  energies = foreach(i = 1:splitLen, .combine = 'c', .inorder = TRUE) %dopar% {
    lowBound = splitIndex[i]+1
      upBound = splitIndex[i+1]
      pieceWP = wavethresh::wp(timeSeries[lowBound:upBound])
      lapply(0:levelLen, function(pkt) getpacket(pieceWP, level = (exp2-levels), index = pkt)) %>% lapply(rms) %>% {total <- sum(unlist(.)); unlist(.) / total}
  }
  # energies = numeric()
  # for(i in 1:splitLen){
  #   lowBound = splitIndex[i]+1
  #   upBound = splitIndex[i+1]
  #   pieceWP = wavethresh::wp(timeSeries[lowBound:upBound])
  #   energy = lapply(0:levelLen, function(pkt) getpacket(pieceWP, level = (exp2-levels), index = pkt)) %>% lapply(rms) %>% {total <- sum(unlist(.)); unlist(.) / total}
  #   energies = c(energies, energy)
  # }
  energies <- matrix(energies, ncol = (levelLen+1), byrow = T)
  return(energies)
}
