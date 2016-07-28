#include <math.h>
#include <iostream>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
double rValue(int alpha, int i, int t, int rType){
  double rResult = 0;
  double divBy = 0;
  //Compute rValue
  switch(rType){
  //Exponential
  case 213:
    rResult = exp(-alpha*(t-i));
    break;
    //Normalized Exponential
  case 190:
    for(int k = 1; k <= t; k++)
      divBy += exp(-alpha*(t-k));
    rResult = t*exp(-alpha*(t-i))/divBy;
    break;
    //Polynomial Type 1
  case 161:
    rResult = pow((i/t),alpha);
    break;
    //Polynomial Type 2
  case 162:
    rResult = pow((t-i+1),alpha);
    break;
    //Normalized Polynomial Type 1
  case 127:
    for(int k = 1; k <= t; k++)
      divBy += pow((k/t),alpha);
    rResult = t*pow((i/t),alpha)/divBy;
    break;
    //Normalized Polynomial Type 2
  case 128:
    for(int k = 1; k <= t; k++)
      divBy += pow((t-k+1),alpha);
    rResult = t*pow((t-i+1),alpha)/divBy;
    break;
  }
  return rResult;
}
// [[Rcpp::export]]
double vValue(double vH, double vV, double weight, int lambda){
  return pow((pow(weight*vH,lambda) + pow((1-weight)*vV,lambda)),(1/lambda));
}

double fValue(int j, int beta, double weight, int lambda, std::vector<double> *dataH, std::vector<double> *dataV, int fType){
  double fResult = 0;
  double sum = 0;
  switch(fType){
  //Polynomial
  case 233:
    fResult = pow(vValue(dataH -> at(j),dataV -> at(j), weight, lambda),beta);
    break;
    //Normalized polynomial
  case 199:
    for(int k = 1; k < dataH -> size(); k++)
      sum += vValue(dataH -> at(k),dataV -> at(k), weight, lambda);
    fResult = pow(sum/dataH->size(),1/beta);
    return fResult;
    break;
  }
  return fResult;
}
double degradation(std::vector<double> dataH, std::vector<double> dataV, int alpha, int beta, int lambda, double weight, int i, int t, int rType, int fType){
  double R = rValue(alpha,i,t,rType);
  double degrad = 0;
  for(int j = 1; j < dataH.size(); j++)
    degrad += fValue(j,beta,weight,lambda,&dataH,&dataV,fType);
  degrad *= R;///dataH.size();
  return degrad;
}

//' Compute the accumulated degradation at each time interval.
//' Accumulated Degradation (accDegrad)
//' @param dataH A matrix containing the horizontal vibration data.
//' @param dataV A matrix containing the vertical vibration data.
//' @param alpha A choosen constant value used when computing the R value. The R value takes into account the influence of the time.
//' @param beta A choosen constant value used when computing the F value. The F value takes into account the influence of the acceleration.
//' @param weight The weight must be less or equal to 1. It is the percent weight of the horizontal data.
//' @param rType The method used to compute the R value. See details for more details.
//' @param fType The method used to compute the F value. See details for more details.
//' @param onlyFinal Default to false. If true returns only the final accumulated degradtion value.
//' @return A vector containing all the accumulated degradation values at each time interval.
//' @details rType: \cr
//' "exp" = Exponential [R = e^(alpha*(t-i))] \cr
//' "Nexp" = Normalized Exponential [R = t*e^(alpha*(t-i)/[sum k(1:t)]e^(alpha(t-k)))] \cr
//' "poly1" = Polynomial 1 [R = (i/t)^alpha] \cr
//' "poly2" = Polynomial 2 [R = (t-i+1)^alpha] \cr
//' "Npoly1" = Normalized Polynomial 1 [R = t*(i/t)^alpha/[sum k(1:t)](i/t)^alpha] \cr
//' "Npoly2" = Normalized Polynomial 2 [R = t*(t-i+1)^alpha/[sum k(1:t)](t-i+1)^alpha]\cr \cr
//' fType: \cr
//' "poly" = Polynomial [F = V(i,j)^beta]
//' "Npoly" = Normalized Polynomial [F = ([sum j(1:data length)]V(i,j)^beta/datalength)^(1/beta)]
//' @export
// [[Rcpp::export]]
NumericVector accDegrad(NumericMatrix dataH, NumericMatrix dataV, int alpha, int beta, int lambda, double weight, std::string rType, std::string fType,bool onlyFinal = false){
  std::vector<double> degrad;
  degrad.reserve(dataH.ncol());
  degrad.resize(dataH.ncol());
  int rTypeN = (int)rType[0]+(int)rType[rType.length()-1];
  int fTypeN = (int)fType[0]+(int)fType[fType.length()-1];
  int tile = 4096; //tbd
#pragma omp parallel for
  for(int ii = 1; ii <= dataH.ncol(); ii+=tile){
    for(int i = ii; i <= dataH.ncol() && i < ii+tile; i++){
      NumericVector tempH = dataH(_,i-1);
      NumericVector tempV = dataV(_,i-1);
      degrad.at(i-1) = degradation(as<std::vector<double>>(tempH),as<std::vector<double>>(tempV), alpha, beta, lambda, weight,i, dataH.ncol(),rTypeN,fTypeN);
    }
  }
#pragma omp parallel for
  for(int ii = 2; ii <= dataH.ncol(); ii+=tile){
    for(int i = ii; i <= dataH.ncol() && i < ii+tile; i++){
      degrad.at(i-1)+=degrad.at(i-2);
    }
  }
  if(!onlyFinal)
    return wrap(degrad);
  return wrap(degrad.at(dataH.ncol()-1));
}

//' Create a Life Table comparing the accumulated Degradation to the time remaining
//' @title Create Life Table (createLifeTab)
//' @param accDeg The output of the 'accDegrad' function. A vector containing the Accumulated Degredation values of the time series data.
//' @param timeInterval The interval of time inbetween each sample of data.
//' @return A Life Table used for RUL computation
//' @export
// [[Rcpp::export]]
List createLifeTab(std::vector<double> accDeg, double timeInterval){
  std::vector<double> intervalsRemain;
  intervalsRemain.reserve(accDeg.size());
  intervalsRemain.resize(accDeg.size());
  double maxTime = timeInterval*(accDeg.size()-1);
  int tile = 4096; //tbd
  for(int ii = 0; ii < accDeg.size(); ii+=tile){
    for(int i = ii; i < accDeg.size() && i < ii+tile; i++){
      intervalsRemain.at(i) = maxTime - i*timeInterval;
    }
  }
  NumericMatrix lifeTab(accDeg.size(),2);
  NumericVector tempAccDeg = wrap(accDeg);
  NumericVector tempInterval = wrap(intervalsRemain);
  lifeTab(_,0) = tempAccDeg;
  lifeTab(_,1) = tempInterval;
  colnames(lifeTab) = CharacterVector::create("accDeg","Life");
  return List::create(Named("LifeTab") = lifeTab, Named("SetNum") = 1);
}

//' Updates Life Table produced from 'createLifeTab' function
//' @title Update Life Table (updateLifeTab)
//' @param lifeTab The life table created from createLifeTab or updateLifeTab.
//' @param accDeg The output of the 'accDegrad' function. A vector containing the Accumulated Degredation values of the time series data.
//' @return A Life Table used for RUL computation.
//' @export
// [[Rcpp::export]]
List updateLifeTab(List lifeTab, NumericVector accDeg){
  NumericMatrix lifeTabMat = lifeTab.at(0);
  int setNum = lifeTab[1];
  arma::mat lifeTabMatA = as<arma::mat>(lifeTabMat);
  lifeTabMatA.col(0)*=setNum;
  arma::vec accDegVecA = as<arma::vec>(accDeg);
  lifeTabMatA.col(0)+=accDegVecA;
  setNum ++;
  lifeTabMatA.col(0)/=setNum;
  return List::create(Named("LifeTab") = wrap(lifeTabMatA), Named("SetNum") = setNum);
}

//' Determines the Remaining useful life of a system based on its accumulated degradation value
//' @title Compute Remaining Useful Life (computeRUL)
//' @param lifeTab The life table created from createLifeTab or updateLifeTab.
//' @param accDegVal The final accumulated degredation value of your data.
//' @return The estimated remaining useful life of the system.
//' @export
// [[Rcpp::export]]
double computeRUL(List lifeTab, double accDegVal){
  NumericMatrix lifeTabMat = lifeTab[0];
  for(int i = 0; i < lifeTabMat.nrow(); i++){
    if(lifeTabMat.at(i,0)<accDegVal)
      return lifeTabMat.at(i,1);
  }
  return 0;
}

//' @export
// [[Rcpp::export]]
List meanAndStd(List probMatX){
  std::vector<double> mean;
  std::vector<double> stdev;
  arma::mat revEmisProb = as<arma::mat>(probMatX[4]);
  std::vector<int> revEmisCount = as<std::vector<int>>(probMatX[5]);
  mean.reserve(revEmisCount.size());
  mean.resize(revEmisCount.size());
  stdev.reserve(revEmisCount.size());
  stdev.resize(revEmisCount.size());
  for(int i = 0; i < revEmisCount.size();i++)
    revEmisProb.col(i)*=revEmisCount.at(i);
  for(int c = 0; c < revEmisCount.size(); c++){
    double meanTemp = 0;
    for(int r = 0; r < revEmisProb.n_rows; r++)
      meanTemp += (r+1)*revEmisProb.at(r,c);
    mean.at(c) = meanTemp/revEmisCount.at(c);
  }
  for(int c = 0; c < revEmisCount.size(); c++){
    double stdTemp = 0;
    for(int r = 0; r < revEmisProb.n_rows; r++){
      if(revEmisProb.at(r,c)!=0)
        stdTemp += pow(revEmisProb.at(r,c)-mean.at(c),2)*revEmisProb.at(r,c);
    }
    stdev.at(c) = sqrt(stdTemp/revEmisCount.at(c));
  }
  return List::create(Named("TransProb") = probMatX[0], Named("TransCount") = probMatX[1], Named("EmisProb") = probMatX[2], Named("EmisCount") = probMatX[3], Named("RevEmisProb") = probMatX[4], Named("RevEmisCount") = probMatX[5], Named("Mean") = wrap(mean), Named("StdDev") = wrap(stdev));
}


