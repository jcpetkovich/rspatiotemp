#include <math.h>
#include <Rcpp.h>
using namespace Rcpp;

//note: 4096 is the largest power of 2 that divides into 20480

//' Computes nodal energy of a single column of a WPD. Called by the 'getEnergy' function.
//' @title Compute single nodal energy (nodeEnergy)
//' @param node A NumericVector containing a single column of the WPD values.
//' @param size The length of each column.
//' @return The nodal energy value of that column of the WPD.
double nodeEnergy(NumericVector node,int size){
  double result = 0;
  for(int i = 0; i < size;i++)
    result += node.at(i)*node.at(i);
  result/=size;
  return sqrt(result);
}

//' Compute nodal energy of a WPD.
//' @title Compute nodal energy (getEnergy)
//' @param wpdData A NumericMatrix containing the WPD.
//' @return A vector of the nodal energy valyes of the WPD.
//' @export
// [[Rcpp::export]]
NumericVector getEnergy(NumericMatrix wpdData){
  std::vector<double> energy;
  energy.reserve(wpdData.ncol());
  energy.resize(wpdData.ncol());
  int tile = 4096;
  int size = wpdData.nrow();
#pragma omp parallel for
  for(int ii = 0; ii < wpdData.ncol(); ii+=tile){
    for(int i = ii ; i < wpdData.ncol() && i < ii+tile;i++){
      NumericVector node = wpdData.column(i);
      energy.at(i) = nodeEnergy(node,size);
    }
  }
  return wrap(energy);
}

//' Group viterbi sequence values accordingly for mean and stdDev calculation
//' @title Format Viterbi Sequence (formatRULHMM)
//' @param seq The viterbi sequence to be formatted
//' @return A list of the reduced sequence and the repeat sequence. The formatted sequence.
//' @export
// [[Rcpp::export]]
List formatRULHMM(std::vector<int> seq){
  std::vector<int> redSeq;
  std::vector<int> repSeq;
  int wpLength =  seq.size();
  int currSAX = seq.at(0);
  int rep = 1;
  for(int i = 1; i <wpLength; i++){
    if(currSAX == seq.at(i))
      rep++;
    else{
      redSeq.push_back(currSAX);
      repSeq.push_back(rep);
      rep = 1;
      currSAX = seq.at(i);
    }
  }
  redSeq.push_back(currSAX);
  repSeq.push_back(rep);
  return List::create(Named("redSeq") = redSeq, Named("repSeq") = repSeq);
}

//' Compute the mean and standard deviation of the formatted viterbi sequence.
//' @title Mean and StdDev of viterbi (meanStdDev)
//' @param redSeq The reduced viterbi sequence from the 'formatRULHMM' function.
//' @param repSeq The repeat viterbi sequence from the 'formatRULHMM' function.
//' @param hidAlphabetSize The alphabet size of the hidden sequence.
//' @param tab Whether or not you are using *.tab functions.
//' @return The mean and standard deviation vectors for each state.
//' @export
// [[Rcpp::export]]
List meanStdDev(std::vector<int> redSeq, std::vector<int> repSeq, int hidAlphabetSize, bool tab = false){
  std::vector<double> mean;
  std::vector<double> stdDev;
  std::vector<double> counter;
  mean.reserve(hidAlphabetSize);
  mean.resize(hidAlphabetSize);
  stdDev.reserve(hidAlphabetSize);
  stdDev.resize(hidAlphabetSize);
  counter.reserve(hidAlphabetSize);
  counter.resize(hidAlphabetSize);

  for(int i = 0; i < hidAlphabetSize; i++)
    counter[i] = 0;
  for(int i = 0; i < redSeq.size(); i++){
    mean[redSeq.at(i)-1] += repSeq.at(i);
    stdDev[redSeq.at(i)-1] += repSeq.at(i)*repSeq.at(i);
    counter[redSeq.at(i)-1] ++;
  }
  for(int i = 0; i < hidAlphabetSize; i++){
    if(counter[i]!=0){
      mean[i] /= counter[i];
      stdDev[i] = sqrt(stdDev[i]/counter[i] - mean[i]*mean[i]);
    }
  }
  if(tab)
    return List::create(Named("mean") = wrap(mean), Named("stdDev") = wrap(stdDev), Named("tab") = wrap(counter));
  return List::create(Named("mean") = wrap(mean), Named("stdDev") = wrap(stdDev));
}

struct node{
  int state;
  double probScore;
  std::vector<int> possibleStates;
  node *prev = NULL;
};

//' Finding the critical path from a start state to an end state given the transition probabilities.
//' @title Find Critical Path (criticalPath)
//' @param transProb A NumericMatrix of the transition probabilities.
//' @param startState The beginning state in the path.
//' @param endState The ending state in the path.
//' @param possibleStates A vector containing the possible states in the path.
//' @return A vector containing the critical path. The most probably path from the start state to the end state where each state is visited at most once.
//' @export
// [[Rcpp::export]]
NumericVector criticalPath(NumericMatrix transProb, int startState, int endState, std::vector<int> possibleStates){
  if(startState == endState)
    return startState;
  double initProb = transProb.at(startState,endState);
  possibleStates.erase(possibleStates.begin()+startState);
  std::vector<node*> paths;
  std::vector<node*> finishedPaths;
  node* head = new node;
  head -> state = startState;
  for(int i = 0; i < possibleStates.size();i++){
    if(transProb.at(startState,possibleStates.at(i))>initProb){
      node *path = new node;
      path -> state = possibleStates.at(i);
      path -> possibleStates = possibleStates;
      path -> possibleStates.erase(path -> possibleStates.begin() + i);
      path -> prev = head;
      path -> probScore = transProb.at(startState,possibleStates.at(i));
      paths.push_back(path);
    }
  }
  if(paths.size()!=0){
    while(true){
      std::vector <node*> nextPath;
      for(int path = 0; path < paths.size();path++){
        for(int state = 0; state < paths.at(path) -> possibleStates.size(); state++){
          if((paths.at(path) -> probScore * transProb.at(paths.at(path) -> state, paths.at(path) -> possibleStates.at(state)))>initProb){
            node *newpath = new node;
            newpath -> state = paths.at(path) -> possibleStates.at(state);
            newpath -> possibleStates = paths.at(path) -> possibleStates;
            newpath -> possibleStates.erase(newpath -> possibleStates.begin() + state);
            newpath -> prev = paths.at(path);
            newpath -> probScore = paths.at(path) -> probScore * transProb.at(paths.at(path) -> state, paths.at(path) -> possibleStates.at(state));
            if(newpath -> state == endState)
              finishedPaths.push_back(newpath);
            else
              nextPath.push_back(newpath);
          }
        }
      }
      paths = nextPath;
      if(paths.size() == 0)
        break;
      else if(paths.at(0) -> possibleStates.size() == 0)
        break;
    }
  }
  std::vector<int> seq;
  int indexMax = -1;
  if(finishedPaths.size()!=0){
    double probMax = finishedPaths.at(0) -> probScore;
    indexMax = 0;
    for(int i = 1; i < finishedPaths.size();i++){
      if(finishedPaths.at(i) -> probScore > probMax){
        probMax = finishedPaths.at(i) -> probScore;
        indexMax = i;
      }
    }
    node* tail = finishedPaths.at(indexMax);
    while(true){
      seq.push_back(tail -> state);
      if(tail -> prev == NULL)
        break;
      tail = tail -> prev;
    }
  }
  else if(indexMax == -1)
    seq = {startState,endState};
  return wrap(seq);
}

//' Compute the RUL and it's upper and lower bounds.
//' @title Compute RUL and bounds (computeRULBounds)
//' @param criticalSeq The critical path from the 'criticalPath' function.
//' @param meanVec The mean vector from the 'meanStdDev' function.
//' @param stdDev The standard deviation vector from the 'meanStdDev' function.
//' @param confCoef A confidence constant used for the computation of the upper and lower bounds.
//' @return The estimated RUL and its lower and upper bounds.
//' @export
// [[Rcpp::export]]
List computeRULBounds(std::vector<int> criticalSeq, std::vector<double> meanVec, std::vector<double> stdDev, double confCoef){
  double lower = 0;
  double mean = 0;
  double upper = 0;
  for(int i = 0; i < criticalSeq.size();i++){
    upper += meanVec.at(criticalSeq.at(i)) + confCoef*stdDev.at(criticalSeq.at(i));
    lower += meanVec.at(criticalSeq.at(i)) - confCoef*stdDev.at(criticalSeq.at(i));
    mean += meanVec.at(criticalSeq.at(i));
  }
  return List::create(Named("Lower") = lower, Named("Mean") = mean, Named("Upper") = upper);
}

//' Compute the probability of that the given continuous observed and discrete hidden sequences occuring in a trained system.
//' @title Viterbi Probability Value (viterbiProbVal)
//' @param transProb A matrix containing the transition probabilites from the hidden markov model.
//' @param emisProb A matrix containing the emission probabilities from the hidden markov model.
//' @param obsSeq A vector containing the observed sequence of data.
//' @param hidSeq A vector containing the hidden sequence of data.
//' @return The probability of the given continous observed and discrete hidden sequences occuring in a trained system.
//' @export
// [[Rcpp::export]]
NumericVector viterbiProbDepmix(NumericMatrix transProb, NumericMatrix emisProb, std::vector<double> obsSeq, std::vector<int> hidSeq){
  std::vector<double> forwProb;
  forwProb.reserve(obsSeq.size());
  forwProb.resize(obsSeq.size());
  double obsVal = fabs(emisProb.at(0,hidSeq.at(0)-1)-obsSeq.at(0)) + emisProb.at(0,hidSeq.at(0)-1);
  forwProb.at(0) = log(R::pnorm5(obsVal,emisProb.at(0,hidSeq.at(0)-1),emisProb.at(1,hidSeq.at(0)-1),1,0));
  for(int i = 1; i < obsSeq.size(); i++){
    double obsVal = fabs(emisProb.at(0,hidSeq.at(i)-1)-obsSeq.at(i)) + emisProb.at(0,hidSeq.at(i)-1);
    forwProb.at(i) = log(R::pnorm5(obsVal,emisProb.at(0,hidSeq.at(i)-1),emisProb.at(1,hidSeq.at(i)-1),1,0)) + log(transProb.at(hidSeq.at(i-1)-1,hidSeq.at(i)-1)) + forwProb.at(i-1);
  }
  return wrap(forwProb.at(forwProb.size()-1));
}

struct nodeVit{
  nodeVit* prev;
  int state = 0;
  double logProb;
};

//' @export
// [[Rcpp::export]]
NumericVector viterbiCont(NumericMatrix transProb, NumericMatrix emisProb, std::vector<double> obsSeq){
  nodeVit** prevNodes = new nodeVit*[transProb.nrow()];
  nodeVit** tailNodes1 = new nodeVit*[transProb.nrow()];
  for(int i = 0; i < transProb.nrow(); i++){
    nodeVit *tempNode = new nodeVit;
    tempNode -> prev = NULL;
    tempNode -> state = i;
    double obsVal = 1 - fabs(emisProb.at(0,i)-obsSeq.at(0)) + emisProb.at(0,i);
    tempNode -> logProb = log(R::pnorm5(obsVal,emisProb.at(0,i),emisProb.at(1,i),1,0));
    tailNodes1[i] = tempNode;
  }
  prevNodes = tailNodes1;
  for(int ob = 1; ob < obsSeq.size(); ob++){
    nodeVit** tailNodes = new nodeVit*[transProb.nrow()];
    for(int row = 0; row < transProb.nrow(); row++){
      double maxLogProb = (prevNodes[0] -> logProb) + log(transProb.at(0,row));
      int prevState = 0;
      for(int prevRow = 1; prevRow < transProb.nrow(); prevRow++){
        if(((prevNodes[prevRow] -> logProb) + log(transProb.at(prevRow,row))) > maxLogProb){
          maxLogProb = (prevNodes[prevRow] -> logProb) + log(transProb.at(prevRow,row));
          prevState = prevRow;
        }
      }
      nodeVit *tempNode = new nodeVit;
      tempNode -> prev = prevNodes[prevState];
      tempNode -> state = row;
      double obsVal = 1 - fabs(emisProb.at(0,row)-obsSeq.at(ob)) + emisProb.at(0,row);
      tempNode -> logProb = maxLogProb + log(R::pnorm5(obsVal,emisProb.at(0,row),emisProb.at(1,row),1,0));
      tailNodes[row] = tempNode;
    }
    prevNodes = tailNodes;
  }
  double maxProbVal = prevNodes[0] -> logProb;
  int finalState = 0;
  for(int i = 1; i < transProb.nrow(); i++){
    if((prevNodes[i] -> logProb) > maxProbVal){
      maxProbVal = prevNodes[i] -> logProb;
      finalState = i;
    }
  }
  std::vector<int> vitSeq;
  nodeVit *tailNode = prevNodes[finalState];
  vitSeq.push_back(tailNode -> state);
  while(true){
    if(tailNode -> prev == NULL)
      break;
    tailNode = tailNode -> prev;
    vitSeq.push_back(tailNode -> state);
  }
  std::reverse(vitSeq.begin(),vitSeq.end());
  return(wrap(vitSeq));
}
