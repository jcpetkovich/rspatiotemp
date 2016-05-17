#include <iostream>
#include <math.h>
#include <omp.h>
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

//Possibility for additional training after initial training
//Baye's Theorem? (Probably won't)
//Instead return the vector of the size for each row of probabilities.
//Calculate new probability from there


//Convert Original data into viterbi friendly form, by taking the groupings and giving it the new alphabetSize
// [[Rcpp::export]]
RObject convert(std::vector<int> data, int groupSize, int alphabetSize){
  std::vector<int> convertData;

  int tile = 4096; //tbd

  convertData.reserve(data.size() - groupSize+1);
  convertData.resize(data.size() - groupSize+1);

#pragma omp parallel for
  for(int ii = 0; ii < (data.size()-groupSize+1); ii+= tile){
    for(int i = ii; i < (ii+ tile) && i < (data.size() - groupSize+1); i++){
      int index = 0;
      for(int k = 0; k < groupSize; k++)
        index += data.at(i+k) * std::pow(alphabetSize, groupSize - k - 1);
      convertData.at(i) = index;
    }
  }
  return wrap(convertData);
}

//This only output transition probabilities for itself
//It takes in one sequence of data and produces and matrix of probabilities to predict the next values in the same set of data
List createProbMat(std::vector<int> data, int groupSize, int alphabetSize){
  //create empty Matrix for prob storing
  int rowSize = std::pow(alphabetSize,groupSize);
  arma::mat probMat(rowSize,alphabetSize);
  arma::vec counter(rowSize);
  probMat.zeros();
  counter.zeros();
  int tile = 4096; //tbd
#pragma omp parallel for
  for(int ii = 0; ii < (data.size() - groupSize); ii+=tile){
    for(int i = ii; i < (ii+tile) && i < (data.size()-groupSize); i++){
      int indexRow = 0;
      for(int k = 0; k < groupSize; k++)
        indexRow += data.at(k+i) * std::pow(alphabetSize,groupSize - k - 1);
      probMat(indexRow, data.at(i+groupSize))++;
      counter(indexRow)++;
    }
  }
#pragma omp parallel for
  for(int ii = 0; ii < rowSize; ii+=tile){
    for(int i = ii; i < (ii+tile) && i < rowSize; i++){
      if(counter(i) != 0)
        probMat.row(i) /= counter(i);
    }
  }
  return List::create(Named("Probability") = probMat, Named("Counter") = counter, Named("AlphabetSize") = alphabetSize, Named("GroupSize") = groupSize);
}

//createProbMat(std::vector<int> Observed Data, std::vector<int> Hidden Data, int Group Size, int Observed Alphabet Size, int Hidden Alphabet Size)
// [[Rcpp::export]]
List createProbMatX(std::vector<int> dataO, std::vector<int> dataH, int groupSize, int alphabetSizeO, int alphabetSizeH){
  int rowSize = std::pow(alphabetSizeO,groupSize);
  List transData = createProbMat(dataH,groupSize,alphabetSizeH);
  arma::mat emisProb(rowSize, alphabetSizeH);
  arma::vec counter(rowSize);
  counter.zeros();
  emisProb.zeros();

  int tile = 4096; //tbd
#pragma omp parallel for
  for(int ii = 0; ii < (dataO.size() - groupSize); ii+=tile){
    for(int i = ii; i < (ii+tile) && i < (dataO.size() - groupSize + 1); i++){
      int indexRow = 0;
      for(int k = 0; k < groupSize; k++)
        indexRow += dataO.at(k+i) * std::pow(alphabetSizeO,groupSize - k - 1);
      emisProb(indexRow,dataH.at(i+1))++;
      counter(indexRow)++;
    }
  }
#pragma omp parallel for
  for(int ii = 0; ii < rowSize; ii+=tile){
    for(int i = ii; i < (ii+tile) && i < rowSize; i++){
      if(counter(i) != 0)
        emisProb.row(i) /= counter(i);
    }
  }
  return List::create(Named("Transition") = transData, Named("Emission") = List::create(Named("Probability") = emisProb, Named("Counter") = counter, Named("AlphabetSize") = alphabetSizeO, Named("GroupSize") = groupSize));
}

//update probabilities. input previously calculated probabilities and new parameters to update the previous probabilities
// [[Rcpp::export]]
List updateProb(List prevData, std::vector<int> dataO, std::vector<int> dataH, int groupSize, int alphabetSizeO, int alphabetSizeH){
  List prevTrans = prevData[0];
  List prevEmis = prevData[1];
  int prevASizeH = prevTrans[2];
  int prevASizeO = prevEmis[2];
  int prevGSize = prevTrans[3];
  if((prevASizeH != alphabetSizeH)||(prevASizeO != alphabetSizeO)||(prevGSize != groupSize)){
    Rcerr<<"Something went wrong, one or more of the sizes don't match\n";
    return prevData;
  }
  else{
    List newData = createProbMatX(dataO,dataH,groupSize,alphabetSizeO,alphabetSizeH);
    arma::mat prevTransProb = prevTrans[0];
    arma::vec prevTransCount = prevTrans[1];
    arma::mat prevEmisProb = prevEmis[0];
    arma::vec prevEmisCount = prevEmis[1];

    List newTrans = newData[0];
    List newEmis = newData[1];
    arma::mat newTransProb = newTrans[0];
    arma::vec newTransCount = newTrans[1];
    arma::mat newEmisProb = newEmis[0];
    arma::vec newEmisCount = newEmis[1];

    int tile = 4096; //tbd

#pragma omp parallel for
    for(int ii = 0; ii < prevTransProb.n_rows; ii+=tile){
      for(int i = ii; i < (ii+tile) && i < prevTransProb.n_rows; i++){
        prevTransProb.row(i) *= prevTransCount.at(i);
        newTransProb.row(i) *= newTransCount.at(i);
      }
    }
#pragma omp parallel for
    for(int ii = 0; ii < prevEmisProb.n_rows; ii+=tile){
      for(int i = 0; i < (ii+tile) && i < prevEmisProb.n_rows;i++){
        prevEmisProb.row(i) *= prevEmisCount.at(i);
        newEmisProb.row(i) *= newEmisCount.at(i);
      }
    }
    newTransProb += prevTransProb;
    newEmisProb += prevEmisProb;
    newTransCount += prevTransCount;
    newEmisCount += prevEmisCount;

#pragma omp parallel for
    for(int ii = 0; ii < prevTransProb.n_rows; ii+=tile){
      for(int i = ii; i < (ii+tile) && i < prevTransProb.n_rows; i++){
        if(newTransCount.at(i) != 0)
          newTransProb.row(i) /= newTransCount.at(i);
      }
    }

#pragma omp parallel for
    for(int ii = 0; ii < prevEmisProb.n_rows; ii+=tile){
      for(int i = 0; i < (ii+tile) && i < prevEmisProb.n_rows;i++){
        if(newEmisCount.at(i) != 0)
          newEmisProb.row(i) /= newEmisCount.at(i);
      }
    }

    List transData = List::create(Named("Probability") = newTransProb, Named("Counter") = newTransCount, Named("AlphabetSize") = prevASizeH, Named("GroupSize") = prevGSize);
    List emisData = List::create(Named("Probability") = newEmisProb, Named("Counter") = newEmisCount, Named("AlphabetSize") = prevASizeO, Named("GroupSize") = prevGSize);
    return List::create(Named("Transition") = transData, Named("Emission") = emisData);
  }
}
