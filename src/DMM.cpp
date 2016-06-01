#include <iostream>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

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
// [[Rcpp::export]]
List createProbMatX(std::vector<int> dataO, std::vector<int> dataH, int groupSize, int alphabetSizeO, int alphabetSizeH){
  arma::mat transProb(alphabetSizeH,alphabetSizeH);
  arma::mat emisProb(pow(alphabetSizeO,groupSize),alphabetSizeH);
  arma::mat revEmisProb(alphabetSizeH,pow(alphabetSizeO,groupSize));

  transProb.zeros();
  emisProb.zeros();
  revEmisProb.zeros();

  std::vector<int> transCounter;
  transCounter.reserve(alphabetSizeH);
  transCounter.resize(alphabetSizeH);

  std::vector<int> emisCounter;
  emisCounter.reserve(alphabetSizeH);
  emisCounter.resize(alphabetSizeH);

  std::vector<int> revEmisCounter;
  revEmisCounter.reserve((int)pow(alphabetSizeO,groupSize));
  revEmisCounter.resize((int)pow(alphabetSizeO,groupSize));

  int tile = 4096;
#pragma omp parallel for
  for(int ii = 1; ii < dataH.size();ii+=tile){
    for(int i = ii; i < (ii+tile) && i < dataH.size(); i++){
      //transition (x given y), transProb(x,y)
      transProb.at(dataH.at(i),dataH.at(i-1))++;
      transCounter[dataH.at(i-1)]++;
      //emission (z given x), emisProb(z,x)
      if(i>= groupSize-1){
        int index = 0;
        for(int k = 0; k < groupSize;k++)
          index += dataO.at(i-k)*pow(alphabetSizeO,k);
        emisProb(index,dataH.at(i))++;
        emisCounter[dataH.at(i)]++;
        revEmisProb(dataH.at(i),index)++;
        revEmisCounter[index]++;
      }
    }
  }

  for(int i = 0; i < pow(alphabetSizeO,groupSize); i++){
    if(revEmisCounter[i]!=0)
      revEmisProb.col(i) /= revEmisCounter[i];
  }
  for(int i = 0; i < alphabetSizeH; i++){
    if(transCounter[i]!=0)
      transProb.col(i) /= transCounter[i];
    if(emisCounter[i]!=0)
      emisProb.col(i) /= emisCounter[i];
  }
  return List::create(Named("TransProb") = wrap(transProb), Named("TransCount") = wrap(transCounter), Named("EmisProb") = wrap(emisProb), Named("EmisCount") = wrap(emisCounter), Named("RevEmisProb") = wrap(revEmisProb), Named("RevEmisCount") = wrap(revEmisCounter));
}

/*
 //update probabilities. input previously calculated probabilities and new parameters to update the previous probabilities
 // [[Rcpp::export]]
 List updateProb(List prevData, std::vector<int> dataO, std::vector<int> dataH, int groupSize, int alphabetSizeO, int alphabetSizeH){
 List prevTrans = prevData[0];
 List prevEmis = prevData[1];
 int prevASizeH = prevTrans[2];
 int prevASizeO = prevEmis[2];
 int prevGSize = prevTrans[3];
 if((prevASizeH != alphabetSizeH)||(prevASizeO != alphabetSizeO)||(prevGSize != groupSize)){
 //Rcerr<<"Something went wrong, one or more of the sizes don't match\n";
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
 prevTransProb += newTransProb;
 prevEmisProb += newEmisProb;
 prevTransCount += newTransCount;
 prevEmisCount += newEmisCount;

#pragma omp parallel for
 for(int ii = 0; ii < prevTransProb.n_rows; ii+=tile){
 for(int i = ii; i < (ii+tile) && i < prevTransProb.n_rows; i++){
 if(prevTransCount.at(i) != 0)
 prevTransProb.row(i) /= prevTransCount.at(i);
 }
 }

#pragma omp parallel for
 for(int ii = 0; ii < prevEmisProb.n_rows; ii+=tile){
 for(int i = 0; i < (ii+tile) && i < prevEmisProb.n_rows;i++){
 if(prevEmisCount.at(i) != 0)
 prevEmisProb.row(i) /= prevEmisCount.at(i);
 }
 }

 List transData = List::create(Named("Probability") = prevTransProb, Named("Counter") = prevTransCount, Named("AlphabetSize") = prevASizeH, Named("GroupSize") = prevGSize);
 List emisData = List::create(Named("Probability") = prevEmisProb, Named("Counter") = prevEmisCount, Named("AlphabetSize") = prevASizeO, Named("GroupSize") = prevGSize);
 return List::create(Named("Transition") = transData, Named("Emission") = emisData);
 }
 }*/

// [[Rcpp::export]]
NumericVector simulateHid(List probMatX, std::vector<int> dataObs, int groupSize, int alphabetSizeO, int alphabetSizeH){
  std::vector<int> simHid;
  simHid.reserve(dataObs.size());
  simHid.resize(dataObs.size());

  NumericMatrix revEmisProb = probMatX.at(4);

  srand (time(NULL));
  //initial values
  for(int i = 0; i < groupSize-1; i++)
    simHid.at(i) = rand() % alphabetSizeH + 0;

  int tile = 4096; //tbd

  //simulate remaining based on revEmisProb and rng
#pragma omp parallel for
  for(int ii = 0; ii < dataObs.size()-groupSize;ii+=tile){
    for(int i = ii; i < (ii+tile) && i < dataObs.size()-groupSize; i++){
      int index = 0;

      for(int g = 0; g < groupSize; g++)
        index += dataObs.at(i+g)*pow(alphabetSizeO,groupSize-g-1);

      float randNum = (float)(rand() % 101 + 0)/100;
      float culmProb = 0;

      for(int k = 0; k < alphabetSizeH; k++){
        culmProb += revEmisProb.at(k,index);

        if(culmProb >= randNum){
          simHid.at(i+groupSize-1) = k;
          break;
        }
      }
    }
  }
  return wrap(simHid);
}

// [[Rcpp::export]]
void tt(){
  srand (time(NULL));
  for(int i = 0; i < 10; i++){
    std::cout<<rand() % 3 + 0<<std::endl;
  }
}
