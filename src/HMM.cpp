#include <iostream>
#include <math.h>
#include <vector>
#include <omp.h>
#include <Rcpp.h>

#include <opencv2/core/core.hpp>
#include "CvHMM.h"

using namespace Rcpp;

//Rcpp::Matrix to a double[]
double* mtxToDbl(NumericMatrix *mtxData){
  int row = mtxData -> nrow();
  int col = mtxData -> ncol();
  double *dblData = new double[row*col];
  int counter = 0;
  for(int r = 0; r < row; r++){
    for(int c = 0; c < col; c++){
      dblData[counter] = mtxData -> at(r,c);
      counter++;
    }
  }
  return dblData;
}

//std::vector to double[]
double* vecToDbl(std::vector<double> &vecData){
  double *dblData = new double[vecData.size()];

  for(int i = 0; i < vecData.size(); i++)
    dblData[i] = vecData.at(i);

  return dblData;
}

//std::vector to int[]
int* vecToInt(std::vector<int> &vecData){
  int *dblData = new int[vecData.size()];

  int tile = 4096; //tbd
#pragma omp parallel for
  for(int ii = 0; ii < vecData.size(); ii+=tile){
    for(int i = ii; i < ii + tile && i < vecData.size(); i++)
      dblData[i] = vecData.at(i);
  }
  return dblData;
}

//cv::Mat to std::vector<int>
std::vector<int> matToVecInt(cv::Mat &Data){
  std::vector<int> data;
  data.reserve(Data.cols);
  data.resize(Data.cols);
  int tile = 4096; //tbd

#pragma omp parallel for
  for(int ii = 0; ii < Data.cols; ii+=tile){
    for(int i = ii; i < ii + tile && i < Data.cols; i++)
      data.at(i) = Data.at<int>(0,i);
  }
  return data;
}

std::vector<double> matToVecDbl(cv::Mat &Data){
  std::vector<double> data;
  data.reserve(Data.cols);
  data.resize(Data.cols);
  int tile = 4096; //tbd

#pragma omp parallel for
  for(int ii = 0; ii < Data.cols; ii+=tile){
    for(int i = ii; i < ii + tile && i < Data.cols; i++)
      data.at(i) = Data.at<double>(0,i);
  }
  return data;
}

//cv::Mat to std::vector<double> used for matrix
std::vector<double> matToMtxDbl(cv::Mat &Data){
  std::vector<double> data;
  data.reserve(Data.cols * Data.rows);
  data.resize(Data.cols * Data.rows);

#pragma omp parallel for
  for(int col = 0; col < Data.cols; col++)
    for(int row = 0; row < Data.rows; row++)
      data.at(row + col*Data.rows) = Data.at<double>(row,col);
  return data;
}

//' Forward algorithm for hidden markov models.
//' @title Forward Algorithm (forward)
//' @param transProb A matrix containing the transition probabilites from the hidden markov model
//' @param emisProb A matrix containing the emission probabilities from the hidden markov model
//' @param initProb A matrix containing the initial probabilities from the hidden markov model
//' @param dataO A vector containing the observed sequence of data
//' @param dataH A single symbol of the hidden data
//' @param indexH The index at which the single symbol of hidden data above is located in it's original sequence
//' @return A double(real) value. The probability that the hidden symbol is at 'indexH' given the observed sequence and probabilitiy matrices
//' @export
// [[Rcpp::export]]
double forward(NumericMatrix transProb, NumericMatrix emisProb, std:: vector<double> initProb, std::vector<int> dataO, int dataH, int indexH){
  //Error Checking
  if(indexH > dataO.size()){
    Rf_error("The index of the hidden data is larger than the size of the observed data\n");
    return 0;
  }
  if(dataH >= transProb.ncol()){
    Rf_error("Your hidden data value is larger than the alphabet Size\n");
    return 0;
  }
  //End of Error Checking

  int m = transProb.ncol();
  int maxLength = 0;
  maxLength = indexH-1;

  NumericMatrix probForw(maxLength,m);
  int tile = 10; //tbd
  //initialize alpha Z0
#pragma omp parallel for
  for(int zz = 0; zz < m; zz+=tile){
    for(int z = zz; z < (zz+tile) && z < m; z++)
      probForw.at(0,z) = initProb.at(z) * emisProb.at(dataO.at(0),z);
  }
  //recursion alpha Z_k-1
  //a_t
  for(int k = 1; k < maxLength; k++){
    //z_t
#pragma omp parallel for
    for(int zz = 0; zz < m; zz+=tile){
      for(int z = zz; z < (zz+tile) && z < m; z++){
        double sum = 0;
        //z_(t-1)
        for(int i = 0; i < m; i++)
          sum += transProb(z,i) * (probForw.at(k-1,i));
        probForw.at(k,z) = sum * emisProb.at(dataO.at(k), z);
      }
    }
  }
  //final
  double sum = 0;
  for(int i = 0; i < m; i++)
    sum += transProb(dataH,i) * probForw(maxLength-1,i);
  double result = sum*emisProb.at(dataO.at(maxLength), dataH);
  if(result == 0)
    Rf_warning("The calculated probability is too small for the double variable type to hold and your result has been rounded to zero");
  return result;
}

//' Compute the probability of that the given observed and hidden sequences occur in a trained system
//' @title Viterbi Probability Value (viterbiProbVal)
//' @param transProb A matrix containing the transition probabilites from the hidden markov model
//' @param emisProb A matrix containing the emission probabilities from the hidden markov model
//' @param initProb A matrix containing the initial probabilities from the hidden markov model
//' @param dataO A vector containing the observed sequence of data
//' @param dataH A vector containing the hidden sequence of data
//' @return A vector containing the probability to get to each index of the hidden state given the observed state
//' @export
// [[Rcpp::export]]
RObject viterbiProbVal(NumericMatrix transProb, NumericMatrix emisProb, std:: vector<double> initProb, std::vector<int> dataO, std::vector<int> dataH, bool logLand = false){
  int maxLength = std::min(dataO.size(),dataH.size());
  std::vector<double> forwProb;
  forwProb.reserve(maxLength);
  forwProb.resize(maxLength);
  if(logLand){
    forwProb.at(0) = log(initProb.at(dataH.at(0))) + log(emisProb(dataH.at(0),dataO.at(0)));
    int tile = 10; //tbd
#pragma omp parallel for
    for(int tt = 1; tt < maxLength; tt+= tile){
      for(int t = tt; t < (tt+tile) && t < maxLength; t++)
        forwProb.at(t) = log(emisProb.at(dataH.at(t),dataO.at(t))) + log(transProb.at(dataH.at(t-1),dataH.at(t))) + forwProb.at(t-1);
    }
    return wrap(forwProb);
  }
  else{
    //initial
    forwProb.at(0) = initProb.at(dataH.at(0)) * emisProb(dataH.at(0),dataO.at(0));
    int tile = 10; //tbd
#pragma omp parallel for
    for(int tt = 1; tt < maxLength; tt+= tile){
      for(int t = tt; t < (tt+tile) && t < maxLength; t++)
        forwProb.at(t) = emisProb.at(dataH.at(t),dataO.at(t)) * transProb.at(dataH.at(t-1),dataH.at(t)) * forwProb.at(t-1);
    }
    return wrap(forwProb);
  }
}
// [[Rcpp::export]]
double maxVec(std::vector<double> data){
  double max = data.at(0);
  for(int i = 1; i < data.size(); i++){
    if(data.at(i) > max)
      max = data.at(i);
  }
  return max;
}

//' @export
// [[Rcpp::export]]
double probObsLog(NumericMatrix transProb, NumericMatrix emisProb, std::vector<double> dataO){
  int m = transProb.ncol();

  NumericMatrix probForw(dataO.size(),m);
  int tile = 10; //tbd
  //initialize alpha Z0
#pragma omp parallel for
  for(int zz = 0; zz < m; zz+=tile){
    for(int z = zz; z < (zz+tile) && z < m; z++){
      double obsVal = fabs(emisProb.at(0,z)-dataO.at(0)) + emisProb.at(0,z);
      probForw.at(0,z) = log(R::pnorm5(obsVal,emisProb.at(0,z),emisProb.at(1,z),1,0));
    }
  }
  //recursion alpha Z_k-1
  //a_t
  for(int k = 1; k < dataO.size(); k++){
    //z_t
#pragma omp parallel for
    for(int zz = 0; zz < m; zz+=tile){
      for(int z = zz; z < (zz+tile) && z < m; z++){
        std::vector<double> sum;
        //z_(t-1)
        for(int i = 0; i < m; i++){
          double obsVal = fabs(emisProb.at(0,z)-dataO.at(i)) + emisProb.at(0,z);
          sum.push_back(log(transProb(z,i)) + probForw.at(k-1,i) + log(R::pnorm5(obsVal,emisProb.at(0,i),emisProb.at(1,i),1,0)));
        }
        double max = maxVec(sum);
        double logProb = 0;
        for(int i = 0; i < m; i++)
          logProb += exp(sum.at(i)-max);
        logProb = log(logProb);
        probForw.at(k,z) = max + logProb;
      }
    }
  }

  double result = 0;
  for(int i = 0; i < m; i++)
    result += probForw(dataO.size()-1,i);
  if(result == 0)
    Rf_warning("The calculated probability is too small for the double variable type to hold and your result has been rounded to zero");
  return result;
}

//' Determine the most probable hidden sequence given an observed sequence
//' @param transProb A matrix containing the transition probabilites from the hidden markov model
//' @param emisProb A matrix containing the emission probabilities from the hidden markov model
//' @param initProb A matrix containing the initial probabilities from the hidden markov model
//' @param dataV A vector containing the observed sequence of data
//' @return A vector containing the most probable sequence of hidden states
//' @export
// [[Rcpp::export]]
RObject viterbi(NumericMatrix transProb, NumericMatrix emisProb, std::vector<double> initProb, std::vector<int> dataV){
  cv::Mat trans = cv::Mat(transProb.ncol(),transProb.nrow(),CV_64F, mtxToDbl(&transProb)).clone();
  cv::Mat emis = cv::Mat(emisProb.ncol(),emisProb.nrow(),CV_64F, mtxToDbl(&emisProb)).clone();
  cv::Mat init = cv::Mat(1,initProb.size(),CV_64F, vecToDbl(initProb)).clone();
  int *data = vecToInt(dataV);
  cv::Mat seq = cv::Mat(1,dataV.size(),CV_32S,data);
  cv::Mat vitData;
  CvHMM hmm;
  hmm.viterbi(seq.row(0),trans,emis,init,vitData);
  return wrap(matToVecInt(vitData));
}

//' Train Hidden Markov model
//' @title Train (train)
//' @param transProb A matrix containing the transition probabilites from the hidden markov model
//' @param emisProb A matrix containing the emission probabilities from the hidden markov model
//' @param initProb A matrix containing the initial probabilities from the hidden markov model
//' @param dataV A vector containing the observed sequence of data
//' @return A list of new probability matrices
//' @export
// [[Rcpp::export]]
List train(NumericMatrix transProb, NumericMatrix emisProb, std::vector<double> initProb, std::vector<int> dataV){
  cv::Mat trans = cv::Mat(transProb.nrow(),transProb.ncol(),CV_64F, mtxToDbl(&transProb)).clone();
  cv::Mat emis = cv::Mat(emisProb.nrow(),emisProb.ncol(),CV_64F, mtxToDbl(&emisProb)).clone();
  cv::Mat init = cv::Mat(1,initProb.size(),CV_64F, vecToDbl(initProb)).clone();
  CvHMM hmm;

  int *data = vecToInt(dataV);
  cv::Mat seq = cv::Mat(1,dataV.size(),CV_32S,data);

  int maxIterations = 1000;
  hmm.train(seq,maxIterations,trans,emis,init);
  NumericVector mtxTrans = wrap(matToMtxDbl(trans));
  mtxTrans.attr("dim") = Dimension(transProb.nrow(),transProb.ncol());

  NumericVector mtxEmis = wrap(matToMtxDbl(emis));
  mtxEmis.attr("dim") = Dimension(emisProb.nrow(),emisProb.ncol());

  NumericVector vecInit = wrap(matToVecDbl(init));
  return List::create(Named("Transition") = mtxTrans, Named("Emission") = mtxEmis, Named("Initial") = vecInit);
}
