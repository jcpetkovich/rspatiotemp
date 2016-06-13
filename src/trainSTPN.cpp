#include <math.h>
#include <omp.h>
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;

std::map<int,int> factTable;

void fact(int num, bool add = true){
  int diff = 1;
  if(add == false)
    diff = -1;

  for(int i = 2; i <=num; i++){
    if(factTable.find(i) == factTable.end())
      factTable.insert (std::pair<int,int> (i,diff));

    else{
      factTable.find(i) -> second +=diff;
      if(factTable.find(i) -> second == 0)
        factTable.erase(i);
    }
  }
}

double compute(){
  double val = 0;
  for(std::map<int,int>::iterator it = factTable.begin(); it != factTable.end(); it++)
    val += log10(it -> first)*(it -> second);
  return val;
}

int rowSum(NumericVector rowVec){
  int sum = 0;
#pragma omp parallel for reduction (+:sum)
  for(int i = 0; i < rowVec.length(); i++)
    sum += rowVec.at(i);
  return sum;
}

//' Count the number of times an Observed state 'm' is followed by a Hidden Symbol 'n'
//' @title Count Occurrences (count)
//' @param dataO A vector containing observed data
//' @param dataH A vector containing hidden data
//' @param oSize The alphabet size/cardinality of the observed data
//' @param hSize The alphabet size/cardinality of the hidden data
//' @return A NumericMatrix of times each combination of Observed followed by Hidden has occured
//' @export
// [[Rcpp::export]]
NumericMatrix count(std::vector<int> dataO, std::vector<int> dataH, int oSize, int hSize){
  arma::mat tab(oSize,hSize);
  tab.zeros();
  int maxIt = std::min(dataO.size(), dataH.size()-1);
  int tile = 4096; //tbd
#pragma omp parallel for
  for(int ii = 0; ii < maxIt; ii+=tile)
    for(int i = ii; i < (ii+tile) && i < maxIt; i++)
      tab(dataO.at(i),dataH.at(i+1))++;
  return wrap(tab);
}

//' A causality measurement between two sets of data
//' @title Causality Measurement (cosMeasure)
//' @param tabTrained A NumericMatrix returned from the 'count' function for the "training" data
//' @param testObs A vector containing the test observed data
//' @param testHid A vector containing the test hidden data
//' @param testOSize The alphabet size/cardinality of the observed data
//' @param testHSize The alphabet size/cardinality of the hidden data
//' @export
// [[Rcpp::export]]
double cosMeasure(NumericMatrix tabTrained, std::vector<int> testObs, std::vector<int> testHid, int testOSize, int testHSize){

  NumericMatrix tabTest = count(testObs,testHid,testOSize,testHSize);
  /*for(int r = 0; r < testOSize; r++){
   for(int c = 0; c < testHSize; c++){
   Rcout<<"["<<r<<","<<c<<"]"<<fact(tabTrained.at(r,c) + tabTest.at(r,c))/(fact(tabTrained.at(r,c))*fact(tabTest.at(r,c)))<<std::endl;
   }
  }*/
  //#pragma omp parallel for reduction(*:val)
  for(int r = 0; r < testOSize; r++){
    fact(rowSum(tabTest.row(r)));
    fact(rowSum(tabTrained.row(r))+testHSize-1);
    fact(rowSum(tabTest.row(r)) + rowSum(tabTrained.row(r))+testHSize-1,false);
    for(int c = 0; c < testHSize; c++){
      fact(tabTrained.at(r,c) + tabTest.at(r,c));
      fact(tabTrained.at(r,c),false);
      fact(tabTest.at(r,c),false);
    }
  }
  //print();
  double val = compute();
  //Rcout<<compute()<<std::endl;
  factTable.clear();
  return val;
}

//' Compute the causality measurement between several subsequences and a trained set of data
//' @title Compute Distribution (distribution)
//' @param tabTrained Trained count matrix created from 'count' function
//' @param testObs A vector containing the test observed data
//' @param testHid A vector containing the test hidden data
//' @param testOSize The alphabet size of the observed data
//' @param testHSize The alphabet size of the hidden data
//' @param step The step size for subsequence sampling
//' @param subSeqSize The size of each subsequence extracted from test sequence
//' @return A vector containing the causality measurement (Lambda) of the several small subsequences sampled
//' @export
// [[Rcpp::export]]
NumericVector distribution(NumericMatrix tabTrained, std::vector<int> testObs, std::vector<int> testHid, int testOSize, int testHSize, int step, int subSeqSize){
  std::vector<double> cosMeasDis;
  cosMeasDis.reserve((int)((testObs.size()-subSeqSize)/step)+1);
  cosMeasDis.resize((int)((testObs.size()-subSeqSize)/step)+1);
  for(int i = 0; i < testObs.size()-subSeqSize; i+=step){
    //Rcout<<"i: "<<i<<" ";

    std::vector<int> subSeqObs;
    std::vector<int> subSeqHid;
    //Rcout<<"A";
    subSeqObs.assign(testObs.begin()+i,testObs.begin()+i+subSeqSize);
    //Rcout<<"B";
    subSeqHid.assign(testHid.begin()+i,testHid.begin()+i+subSeqSize);
    //Rcout<<"C"<<std::endl;
    cosMeasDis.at(i/step) = cosMeasure(tabTrained, subSeqObs,subSeqHid, testOSize, testHSize);
  }
  //Rcout<<cosMeasDis.size()<<std::endl;
  return wrap(cosMeasDis);
}


//' Compute energy of a binary vector of visible data and a binary vector of hidden data
//' @title Energy (E)
//' @param v A vector containing the binary visible data
//' @param h A vector containing the binary hidden data
//' @param a A vector containing the bias values for the visible data
//' @param b A vector containing the bias values for the hidden data
//' @param w A matrix containing the weights between the visible and hidden data
//' @return A double value. The energy of the binary vector of visible data and binary vector of hidden data
//' @export
// [[Rcpp::export]]
double E(std::vector<int> v, std::vector<int> h, std::vector<double> a, std::vector<double> b, NumericMatrix w){
  double val = 0;
  int maxSide = 0; // v > h
  if(h.size()>v.size())
    maxSide = 1;

  int maxSize = std::max(v.size(),h.size());
  int minSize = std::min(v.size(),h.size());
#pragma omp parallel for reduction (-:val)
  for(int i = 0; i < maxSize;i++){
    if((i<v.size())&&v.at(i))
      val -= a.at(i);
    if(i<h.size()&&h.at(i))
      val -= b.at(i);
    if(maxSide){
      for(int j = 0; j < minSize;j++)
        if(v.at(j)&&h.at(i))
          val-= w.at(j,i);
    }
    else{
      for(int j = 0; j < minSize;j++)
        if(v.at(i)&&h.at(j))
          val-= w.at(i,j);
    }
  }
  return val;
}

//' Compute free energy of a binary vector of visible data for all possible hidden data of chosen size
//' @title Free Energy (F)
//' @param v A vector containing the binary visible data
//' @param hSize The size of vector containing the hidden data
//' @param a A vector containing the bias values for the visible data
//' @param b A vector containing the bias values for the hidden data
//' @param w A matrix containing the weights between the visible and hidden data
//' @return A double value. The free energy for that vector of visible data
//' @export
// [[Rcpp::export]]
double F(std::vector<int> v, int hSize, std::vector<double> a, std::vector<double> b, NumericMatrix w){
  double val = 0;
  int maxIt = pow(2,hSize);

  std::vector<int> h(hSize,0);
  h.at(hSize-1) = -1;
  for(int i = 0; i < maxIt;i++){
    int index = hSize-1;
    h.at(index)++;
    while(h.at(index) > 1){
      index--;
      h.at(index)++;
      h.at(index+1)-=2;
    }
    val += exp(E(v,h,a,b,w));
  }

  return log(val)*-1;
}

