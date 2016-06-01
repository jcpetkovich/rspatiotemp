#include <iostream>
#include <math.h>
#include <vector>
#include <omp.h>
#include <Rcpp.h>

#include <opencv2/core/core.hpp>
#include "CvHMM.h"

using namespace Rcpp;

//Rcpp::Matrix to a double[]
double* mtxToDbl(NumericMatrix &mtxData){
  int row = mtxData.nrow();
  int col = mtxData.ncol();
  double *dblData = new double[row*col];

#pragma omp parallel for
  for(int c = 0; c < col; c++){
    for(int r = 0; r < row; r++)
      dblData[c*row + r] = mtxData.at(r,c);
  }
  return dblData;
}

//std::vector to double[]
double* vecToDbl(std::vector<double> &vecData){
  double *dblData = new double[vecData.size()];
  int tile = 4096; //tbd

#pragma omp parallel for
  for(int ii = 0; ii < vecData.size(); ii+= tile){
    for(int i = ii; i < ii + tile && i < vecData.size(); i++)
      dblData[i] = vecData.at(i);
  }

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
      data.at(row*Data.cols + col) = Data.at<double>(row,col);
  return data;
}

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

// [[Rcpp::export]]
RObject viterbiProbVal(NumericMatrix transProb, NumericMatrix emisProb, std:: vector<double> initProb, std::vector<int> dataO, std::vector<int> dataH){
  int maxLength = std::min(dataO.size(),dataH.size());
  NumericVector forwProb(maxLength);
  //initial
  forwProb.at(0) = initProb.at(dataH.at(0) * emisProb(dataO.at(0),dataH.at(0)));
  int tile = 10; //tbd
#pragma omp parallel for
  for(int tt = 1; tt < maxLength; tt+= tile){
    for(int t = tt; t < (tt+tile) && t < maxLength; t++)
      forwProb.at(t) = emisProb.at(dataO.at(t),dataH.at(t)) * transProb.at(dataH.at(t),dataH.at(t-1)) * forwProb.at(t-1);
  }
  return forwProb;
}

// [[Rcpp::export]]
RObject viterbi(NumericMatrix transProb, NumericMatrix emisProb, std::vector<double> initProb, std::vector<int> dataV){

  //both transition and emission probablities will probably need to be matrices
  cv::Mat trans = cv::Mat(transProb.ncol(),transProb.nrow(),CV_64F, mtxToDbl(transProb)).clone();

  cv::Mat emis = cv::Mat(emisProb.ncol(),emisProb.nrow(),CV_64F, mtxToDbl(emisProb)).clone();

  //It is probably okay if the initial Probablities
  cv::Mat init = cv::Mat(1,initProb.size(),CV_64F, vecToDbl(initProb)).clone();

  //also okay if input data is a vector
  int *data = vecToInt(dataV);
  cv::Mat seq = cv::Mat(1,dataV.size(),CV_32S,data);

  //print input data
  //std::cout<<"input data: ";
  //for (int j=0;j<seq.cols;j++)
  //  std::cout << seq.at<int>(0,j);
  //std::cout << "\n";

  cv::Mat vitData;
  CvHMM hmm;
  //hmm.printModel(trans,emis,init);
  hmm.viterbi(seq.row(0),trans,emis,init,vitData);

  //print viterbi outputs
  //for (int i=0;i<vitData.cols;i++)
  //    std::cout << vitData.at<int>(0,i);
  //std::cout << "\n";
  //std::cout<<"Output Data: \n";
  return wrap(matToVecInt(vitData));
}

// [[Rcpp::export]]
RObject train(NumericMatrix transProb, NumericMatrix emisProb, std::vector<double> initProb, std::vector<int> dataV){
  //both transition and emission probablities will probably need to be matrices
  cv::Mat trans = cv::Mat(transProb.ncol(),transProb.nrow(),CV_64F, mtxToDbl(transProb)).clone();

  cv::Mat emis = cv::Mat(emisProb.ncol(),emisProb.nrow(),CV_64F, mtxToDbl(emisProb)).clone();

  //It is probably okay if the initial Probablities
  cv::Mat init = cv::Mat(1,initProb.size(),CV_64F, vecToDbl(initProb)).clone();

  //also okay if input data is a vector
  int *data = vecToInt(dataV);
  cv::Mat seq = cv::Mat(1,dataV.size(),CV_32S,data);

  CvHMM hmm;
  int maxIterations = 1000;
  hmm.train(seq,maxIterations,trans,emis,init);
  //hmm.printModel(trans,emis,init);

  NumericVector mtxTrans = wrap(matToMtxDbl(trans));
  mtxTrans.attr("dim") = Dimension(transProb.nrow(),transProb.ncol());

  NumericVector mtxEmis = wrap(matToMtxDbl(emis));
  mtxEmis.attr("dim") = Dimension(emisProb.nrow(),emisProb.ncol());

  NumericVector vecInit = wrap(matToVecDbl(init));
  return List::create(Named("Transition") = mtxTrans, Named("Emission") = mtxEmis, Named("Initial") = vecInit);
}
