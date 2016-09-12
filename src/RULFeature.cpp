#include <Rcpp.h>
#include <omp.h>
using namespace Rcpp;

//' @export
// [[Rcpp::export]]
double mean(std::vector<double> data){
  double sum = 0;
  if(data.size() != 0){
#pragma omp parallel for reduction(+:sum)
    for(int i = 0; i < data.size(); i++)
      sum += data.at(i);
    return sum/data.size();
  }
  return 0;
}

//' @export
// [[Rcpp::export]]
double stdDev(std::vector<double> data, double mean){
  double sum = 0;
  if(data.size() != 0){
#pragma omp parallel for reduction(+:sum)
    for(int i = 0; i < data.size(); i++)
      sum += (data.at(i)-mean)*(data.at(i)-mean);
    return sqrt(sum/data.size());
  }
  return 0;
}

//' @export
// [[Rcpp::export]]
double rms(std::vector<double> data){
  double sum = 0;
  if(data.size() != 0){
#pragma omp parallel for reduction(+:sum)
    for(int i = 0; i < data.size(); i++)
      sum += data.at(i)*data.at(i);
    return sqrt(sum/data.size());
  }
  return 0;
}

//' @export
// [[Rcpp::export]]
double kurtosis(std::vector<double> data, double mean, double stdDev){
  double sum = 0;
  if((data.size() != 0)&&(stdDev != 0)){
#pragma omp parallel for reduction(+:sum)
    for(int i = 0; i < data.size(); i++)
      sum += (data.at(i)-mean)*(data.at(i)-mean)*(data.at(i)-mean);
    return sum/(stdDev*stdDev*stdDev*stdDev);
  }
  return 0;
}

//' @export
// [[Rcpp::export]]
double skewness(std::vector<double> data, double mean, double stdDev){
  double sum = 0;
  if((data.size() != 0)&&(stdDev != 0)){
#pragma omp parallel for reduction(+:sum)
    for(int i = 0; i < data.size(); i++)
      sum += (data.at(i)-mean)*(data.at(i)-mean)*(data.at(i)-mean);
    return sum/((data.size()-1)*stdDev*stdDev*stdDev);
  }
  return 0;
}

//' @export
// [[Rcpp::export]]
double lineIntegral(std::vector<double> data){
  double sum = 0;
  if(data.size() != 0){
#pragma omp parallel for reduction(+:sum)
    for(int i = 0; i < data.size()-1; i++)
      sum += fabs(data.at(i+1)-data.at(i));
    return sum;
  }
  return 0;
}

//' @export
// [[Rcpp::export]]
double energy(std::vector<double> data){
  double sum = 0;
  if(data.size() != 0){
#pragma omp parallel for reduction(+:sum)
    for(int i = 0; i < data.size(); i++)
      sum += data.at(i)*data.at(i);
    return sum;
  }
  return 0;
}

//' @export
// [[Rcpp::export]]
List allFeatures(std::vector<double> dataH, std::vector<double> dataV){
  int len = dataH.size();
  //mean
  double sumMeanH = 0;
  double sumMeanV = 0;
  //pk2pk
  double peakMaxH [2]  = {0,0};
  double peakMinH [2]  = {0,0};
  int isPositiveH = -1;
  double peakMaxV [2] = {0,0};
  double peakMinV [2] = {0,0};
  int isPositiveV = -1;

  //rms
  double sumRMSH = 0;
  double sumRMSV = 0;

  //line integral
  double sumLIH = 0;
  double sumLIV = 0;

  for(int i = 0; i < len; i++){
    //mean
    sumMeanH += dataH.at(i);
    sumMeanV += dataV.at(i);

    //pk2pk
    if(i != 0){
      if(isPositiveH == 1){
        if(dataH[i] < dataH[(i-1)]){
          peakMaxH[0] = peakMaxH[0] + dataH[(i-1)];
          peakMaxH[1] = peakMaxH[1] + 1;
          isPositiveH = 0;
        }
      }
      else if(isPositiveH == 0){
        if(dataH[i] > dataH[(i-1)]){
          peakMinH[0] = peakMinH[0] + dataH[(i-1)];
          peakMinH[1] = peakMinH[1] + 1;
          isPositiveH = 1;
        }
      }
      else{
        if(dataH[i] < dataH[(i-1)])
          isPositiveH = 0;
        else if(dataH[i] > dataH[(i-1)])
          isPositiveH = 1;

      }
      if(isPositiveV == 1){
        if(dataV[i] < dataV[(i-1)]){
          peakMaxV[0] = peakMaxV[0] + dataV[(i-1)];
          peakMaxV[1] = peakMaxV[1] + 1;
          isPositiveV = 0;
        }
      }
      else if(isPositiveV == 0){
        if(dataV[i] > dataV[(i-1)]){
          peakMinV[0] = peakMinV[0] + dataV[(i-1)];
          peakMinV[1] = peakMinV[1] + 1;
          isPositiveV = 1;
        }
      }
      else{
        if(dataV[i] < dataV[(i-1)])
          isPositiveV = 0;
        else if(dataV[i] > dataV[(i-1)])
          isPositiveV = 1;
      }
    }

    //rms
    sumRMSH += dataH.at(i) * dataH.at(i);
    sumRMSV += dataV.at(i) * dataV.at(i);

    //line integral
    if(i != (dataH.size()-1)){
      sumLIH += fabs(dataH.at(i+1)-dataH.at(i));
      sumLIV += fabs(dataV.at(i+1)-dataV.at(i));
    }

    //energy
  }
  //mean
  double meanH = sumMeanH/dataH.size();
  double meanV = sumMeanV/dataV.size();

  //pk2pk
  double meanPeakMaxH, meanPeakMinH,meanPeakMaxV,meanPeakMinV;
  if(peakMaxH[1]!=0)
    meanPeakMaxH = peakMaxH[0]/peakMaxH[1];
  if(peakMinH[1]!=0)
    meanPeakMinH = peakMinH[0]/peakMinH[1];
  if(peakMaxV[1]!=0)
    meanPeakMaxV = peakMaxV[0]/peakMaxV[1];
  if(peakMinV[1]!=0)
    meanPeakMinV = peakMinV[0]/peakMinV[1];
  double pk2pkH = meanPeakMaxH - meanPeakMinH;
  double pk2pkV = meanPeakMaxV - meanPeakMinV;

  //rms
  double rmsH = sqrtf(sumRMSH/dataH.size());
  double rmsV = sqrtf(sumRMSV/dataV.size());

  //line integral
  //return sumLIH, sumLIV

  //energy
  double energyH = sumRMSH;
  double energyV = sumRMSV;

  //stdDev
  double sdH = 0;
  double sdV = 0;

  //kurtosis
  double kurtosisH = 0;
  double kurtosisV = 0;

  //skewness
  double skewH = 0;
  double skewV = 0;

//#pragma omp parallel for reduction(+:sdH, +:sdV, +:kurtosisH, +:kurtosisV, +:skewH, +:skewV)
  for(int i = 0; i < len; i++){
    //stdDev
    sdH += (dataH.at(i)-meanH)*(dataH.at(i)-meanH);
    sdV += (dataV.at(i)-meanV)*(dataV.at(i)-meanV);
    kurtosisH += (dataH.at(i)-meanH)*(dataH.at(i)-meanH)*(dataH.at(i)-meanH);
    kurtosisV += (dataV.at(i)-meanV)*(dataV.at(i)-meanV)*(dataV.at(i)-meanV);
    skewH += (dataH.at(i)-meanH)*(dataH.at(i)-meanH)*(dataH.at(i)-meanH);
    skewV += (dataV.at(i)-meanV)*(dataV.at(i)-meanV)*(dataV.at(i)-meanV);
  }
  //stdDev
  sdH = sqrtf(sdH/dataH.size());
  sdV = sqrtf(sdV/dataV.size());

  //kurtosis
  kurtosisH = kurtosisH/(sdH*sdH*sdH*sdH);
  kurtosisV = kurtosisV/(sdV*sdV*sdV*sdV);

  //skewness
  skewH = skewH/((dataH.size()-1)*sdH*sdH*sdH);
  skewV = skewV/((dataV.size()-1)*sdV*sdV*sdV);

  return(List::create(Named("pk2pkH") = pk2pkH, Named("pk2pkV") = pk2pkV, Named("rmsH") = rmsH, Named("rmsV") = rmsV, Named("lineIntH") = sumLIH, Named("lineIntV") = sumLIV, Named("energyH") = energyH, Named("energyV") = energyV, Named("kurtosisH") = kurtosisH, Named("kurtosisV") = kurtosisV, Named("skewH") = skewH, Named("skewV") = skewV));
}
