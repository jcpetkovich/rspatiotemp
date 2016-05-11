#include <iostream>
#include <time.h>
#include <vector>
//#include <omp.h>
#include <Rcpp.h>

#include <opencv2/core/core.hpp>
#include "CvHMM.h"

using namespace Rcpp;

//Note: Parallelize these later

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
    int tile = 10; //tbd

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

    int tile = 10; //tbd
#pragma omp parallel for
    for(int ii = 0; ii < vecData.size(); ii+=tile){
        for(int i = ii; i < ii + tile && i < vecData.size(); i++)
            dblData[i] = vecData.at(i);
    }

    return dblData;
}

//cv::Mat to std::vector
std::vector<int> matToVec(cv::Mat &vitData){
    std::vector<int> data;

    data.reserve(vitData.cols);
    data.resize(vitData.cols);
    int tile = 10; //tbd

#pragma omp parallel for
    for(int ii = 0; ii < vitData.cols; ii+=tile){
        for(int i = ii; i < ii + tile && i < vitData.cols; i++)
            data.at(i) = vitData.at<int>(0,i);
    }

    return data;
}

// [[Rcpp::export]]
RObject viterbi(NumericMatrix transProb, NumericMatrix emisProb, std::vector<double> initProb, std::vector<int> dataV){

    //both transition and emission probablities will probably need to be matrices
    cv::Mat trans = cv::Mat(transProb.ncol(),transProb.nrow(),CV_64F, mtxToDbl(transProb)).clone();

    cv::Mat emis = cv::Mat(emisProb.ncol(),emisProb.nrow(),CV_64F, mtxToDbl(emisProb)).clone();

    //It is probably okay if the initial Probablities
    cv::Mat init = cv::Mat(1/*This should always be 1*/,2/*columns from matrix*/,CV_64F, vecToDbl(initProb)).clone();

    //also okay if input data is a vector
    int *data = vecToInt(dataV);
    cv::Mat seq = cv::Mat(1,25,CV_32S,data);

    //print input data
    std::cout<<"input data: ";
    for (int j=0;j<seq.cols;j++)
        std::cout << seq.at<int>(0,j);
    std::cout << "\n";

    cv::Mat vitData;
    CvHMM hmm;
    hmm.viterbi(seq.row(0),trans,emis,init,vitData);

    //print viterbi outputs
    for (int i=0;i<vitData.cols;i++)
        std::cout << vitData.at<int>(0,i);
    std::cout << "\n";

    return wrap(matToVec(vitData));
}
