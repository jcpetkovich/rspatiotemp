#include <iostream>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

//' Convert a SAX series to another SAX series with chosen group size.
//' @title DMarkov Convert (convert)
//' @param data A vector containing the SAX series to be converted
//' @param groupSize The chosen group size for sampling from original series
//' @param alphabetSize The original alphabet size of the inputted SAX data
//' @return A vector containing the new SAX series of alphabet size 'alphabetSize'^'groupSize'
//' @export
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

//' Create new ProbMatX object from a set of observed and hidden data. This ProbMatX object contains transition, emission and reverse emission probabilities
//' @title Create Probability Matrix X (createProbMatX)
//' @param dataO A vector containing the observed data in SAX
//' @param dataH A vector containing the hidden data in SAX
//' @param groupSize The chosen group size for DMarkov
//' @param alphabetSizeO The original alphabet size of the observed data
//' @param alphabetSizeH The original alphabet size of the hidden data
//' @return A ProbMatX object. A list containing the transition, emission and reverse emission probabilities trained from the input data
//' @export
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


//' Update ProbMatX object with more data
//' @title Update Probability Matrix X (updateProbMatX)
//' @param prevData The ProbMatX object created from 'createProbMatX' function
//' @param dataO the new observed data to be updated into the ProbMatX object
//' @param dataH the new hidden data to be updated into the ProbMatX object
//' @param groupSize The chosen group size for DMarkov. The groupSize must be the same as the one used when creating the ProbMatX object
//' @param alphabetSizeO The original alphabet size of the observed data. The alphabetSize must be the same as the one used when creating the ProbMatX object
//' @param alphabetSizeH The original alphabet size of the hidden data. The alphabetSize must be the same as the one used when creating the ProbMatX object
//' @return A ProbMatX object. A list containing the updated transition, emission and reverse emission probabilities from the input data
//' @export
// [[Rcpp::export]]
List updateProbMatX(List prevData, std::vector<int> dataO, std::vector<int> dataH, int groupSize, int alphabetSizeO, int alphabetSizeH){
  List newData = createProbMatX(dataO,dataH,groupSize,alphabetSizeO,alphabetSizeH);
  arma::mat prevTransProb = as<arma::mat>(prevData.at(0));
  std::vector<int> prevTransCount = as<std::vector<int>>(prevData.at(1));
  arma::mat prevEmisProb = as<arma::mat>(prevData.at(2));
  std::vector<int> prevEmisCount = as<std::vector<int>>(prevData.at(3));
  arma::mat prevRevEmisProb = as<arma::mat>(prevData.at(4));
  std::vector<int> prevRevEmisCount = as<std::vector<int>>(prevData.at(5));

  arma::mat newTransProb = as<arma::mat>(newData.at(0));
  std::vector<int> newTransCount = as<std::vector<int>>(newData.at(1));
  arma::mat newEmisProb = as<arma::mat>(newData.at(2));
  std::vector<int> newEmisCount = as<std::vector<int>>(newData.at(3));
  arma::mat newRevEmisProb = as<arma::mat>(newData.at(4));
  std::vector<int> newRevEmisCount = as<std::vector<int>>(newData.at(5));

  for(int i = 0; i < pow(alphabetSizeO,groupSize); i++){
    prevRevEmisProb.col(i) *= prevRevEmisCount[i];
    newRevEmisProb.col(i) *= newRevEmisCount[i];
  }
  for(int i = 0; i < alphabetSizeH; i++){
    prevTransProb.col(i) *= prevTransCount[i];
    newTransProb.col(i) *= newTransCount[i];
    prevEmisProb.col(i) *= prevEmisCount[i];
    newEmisProb.col(i) *= newEmisCount[i];
  }
  prevTransProb+=newTransProb;
  prevEmisProb+=newEmisProb;
  prevRevEmisProb+=newRevEmisProb;

  for(int i = 0; i < pow(alphabetSizeO,groupSize); i++){
    prevRevEmisCount[i]+=newRevEmisCount[i];
    if(prevRevEmisCount[i]!=0)
      prevRevEmisProb.col(i) /= prevRevEmisCount[i];
  }
  for(int i = 0; i < alphabetSizeH; i++){
    prevTransCount[i]+=newTransCount[i];
    prevEmisCount[i]+=newEmisCount[i];
    if(prevTransCount[i]!=0)
      prevTransProb.col(i) /= prevTransCount[i];
    if(prevEmisCount[i]!=0)
      prevEmisProb.col(i) /= prevEmisCount[i];
  }
  return List::create(Named("TransProb") = wrap(prevTransProb), Named("TransCount") = wrap(prevTransCount), Named("EmisProb") = wrap(prevEmisProb), Named("EmisCount") = wrap(prevEmisCount), Named("RevEmisProb") = wrap(prevRevEmisProb), Named("RevEmisCount") = wrap(prevRevEmisCount));
}

//' Simulate hidden series given an observed series and using the reverse emission probability from ProbMatX object
//' @title Simulate Hidden Series (simulateHid)
//' @param probMatX The ProbMatX object created from 'createProbMatX' function
//' @param dataObs The observed data used to simulate the hidden data
//' @param groupSize The groupSize of the DMarkov model. Should be the same as the one used to create ProbMatX object
//' @param alphabetSizeO The original alphabet size of the observed data. The alphabetSize must be the same as the one used when creating the ProbMatX object
//' @param alphabetSizeH The original alphabet size of the hidden data. The alphabetSize must be the same as the one used when creating the ProbMatX object
//' @return A vector containing the simulated hidden data
//' @export
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
