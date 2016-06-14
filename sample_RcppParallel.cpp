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
  const arma::cube mat;
  const vector<vector<string> > names;
  const int num_mat;
  const vector<int> size_emat;
  
  vector<vector<string> > output;
  
  Sum(const arma::cube &pmat, const int &pnum_mat, const vector<vector<string> > &pnames, 
      const vector<int> psize_emat) : mat(pmat), num_mat(pnum_mat), names(pnames), size_emat(psize_emat) {}
  Sum(const Sum& sum, Split) : mat(sum.mat), num_mat(sum.num_mat), names(sum.names), size_emat(sum.size_emat) {}
  
  void operator()(std::size_t begin, std::size_t end) {
    
    vector<string> temp(num_mat);
    
    arma::vec probs(size_emat[0]);
    arma::vec pout(size_emat[0]);
    
    for(int i=0;i<probs.size();i++) {
      probs[i] = 1.0/size_emat[0];
      pout[i] = i;
    }
    
    arma::vec elm = rsample(pout, 1, false, probs);
    string t0;
    
    for(int p=begin;p<end;p++) {
      
      t0 = names[0][elm[0]];
      
      
      for(int i=0;i<num_mat;i++) {
        
        int j = 0;
        for(j=0;j<size_emat[i];j++) {
          if(names[i][j] == t0) break;
        }
        
        arma::vec prob(size_emat[i]);
        arma::vec pou(size_emat[i]);
        for(int k=0;k<probs.size();k++) {
          prob[k] = mat(i, j, k);
          pou[k] = k;
        }
        
        
        arma::vec elmt = rsample(pou, 1, false, prob);
        string t0 = names[i][elmt[0]];
        
        temp[i] = t0;
        
      }
      output.push_back(temp);  
    }
    
    
  }
  
  void join(const Sum& rhs) { 
    for(int i=0;i<rhs.output.size();i++) {
      output.push_back(rhs.output[i]);
    }
  }
};

// [[Rcpp::export]]
CharacterMatrix parallelVectorSum(List object, int n) {
  
  int num_matrix = object.size();
  int max_dim_mat = 0;
  vector<int> size_emat(num_matrix);
  
  for(int i=0;i<num_matrix;i++) {
    S4  ob = object[i];
    CharacterVector stat = ob.slot("states");
    if(stat.size() > max_dim_mat) max_dim_mat = stat.size();
    size_emat[i] = stat.size();
  }
  
  vector<vector<string> > names(num_matrix, vector<string>(max_dim_mat));
  arma::cube mat(max_dim_mat, max_dim_mat, num_matrix);
  mat.fill(0);
  
  for(int i=0;i<num_matrix;i++) {
    S4  ob = object[i];
    NumericMatrix mtrx = ob.slot("transitionMatrix"); 
    CharacterVector stat = ob.slot("states");
    
    for(int j=0;j<mtrx.nrow();j++) {
      for(int k=0;k<mtrx.ncol();k++) {
        mat(i,j,k) = mtrx(j,k);
      }
      names[i][j] = stat[j];
    }
    
  }
  
   Sum sum(mat, num_matrix, names, size_emat);
  
  parallelReduce(0, n, sum);
  
  long long counter=1;
  for(int i=0;i<sum.output.size();i++) {
    for(int j=0;j<sum.output[i].size();j++) {
      cout<<sum.output[i][j]<<" ";
    } cout<<counter++<<endl;
  }
  
  
  return CharacterMatrix();
  
}
