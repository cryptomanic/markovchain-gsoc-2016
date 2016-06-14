// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "sampler.h"

using namespace Rcpp;
using namespace RcppArmadillo;
using namespace RcppParallel;
using namespace std;

struct Sum : public Worker
{   
  arma::vec output;
  int n;
  
  Sum() : n(0), output(100000) {}
  Sum(const Sum& sum, Split) :n(0), output(100000) {}
  
  void operator()(std::size_t begin, std::size_t end) {
    arma::vec states(5);
    for(int i=0;i<5;i++) states[i] = i;
    
    arma::vec probs(5);
    probs[0]  = 0.1;
    probs[1]  = 0.2;
    probs[2]  = 0.3;
    probs[3]  = 0.3;
    probs[4]  = 0.1;
    
    arma::vec rstat = rsample(states, 1, false, probs);
    output[n++] = rstat[0];
    
  }
  
  void join(const Sum& rhs) { 
    for(int i=0;i<rhs.n;i++) {
      output[n++] = rhs.output[i];
    }
  }
};

// [[Rcpp::export]]
NumericVector parallelVectorSum(int n) {
  
  Sum sum;
  
  parallelReduce(0, n, sum);
  
  NumericVector data = wrap(sum.output);
  int sz = sum.n;
  
  NumericVector ans;
  for(int i=0;i<sz;i++) ans.push_back(data[i]);
  
  return ans;
}
