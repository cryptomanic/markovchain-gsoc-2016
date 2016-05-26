allTransMat <- function(seqMat) {
  # number of sequence
  s <- nrow(seqMat)
  # size of each sequence
  n <- ncol(seqMat)
  # possible states
  sname <- unique(as.vector(seqMat))
  m <- length(sname)
  
  lstTransMat <- list()
  for (i in 1:s) {
    for (j in 1:s) {
      # temporary matrix
      tmat <- matrix(0, nrow = m, ncol = m)
      colnames(tmat) <- sname
      rownames(tmat) <- sname
      
      for (k in 1:(n - 1)) {
        dat <- tmat[seqMat[i, k + 1], seqMat[j, k]]
        dat <- dat + 1
        tmat[seqMat[i, k + 1], seqMat[j, k]] <- dat
      }
      
      csum <- colSums(tmat)
      for (l in 1:length(csum)) {
        if (csum[l] == 0) {
          tmat[, l] <- 1 / m
          csum[l] <- 1
        }
      }
      
      tmat <- t(t(tmat) / csum)
      
      lstTransMat[[s * (i - 1) + j]] <- tmat
    }
  }
  return(lstTransMat)
}

# smat <- matrix(c('a','a','b','b','a','c','a',
#          'b','a','b','c','c','b','a',
#          'a','b','b','a','c','a','b',
#          'c','c','c','c','a','a','a'), nrow = 4, byrow = T)

allSeq2freqProb(seqMat) {
  # number of sequence
  s <- nrow(seqMat)
  # size of each sequence
  n <- ncol(seqMat)
  # possible states
  sname <- unique(as.vector(seqMat))
  m <- length(sname)
  
  
}

library(parallel)
# Calculate the number of cores
no_cores <- detectCores() - 1
# Initiate cluster
cl <- makeCluster(no_cores)

statesNames=c("a","b","c")
mcA<-new("markovchain", states=statesNames, transitionMatrix=matrix(c(0.2,0.5,0.3,0,0.2,0.8,0.1,0.8,0.1),nrow=3, byrow=TRUE, dimnames=list(statesNames,statesNames)))
mcB<-new("markovchain", states=statesNames, transitionMatrix=matrix(c(0.2,0.5,0.3,0,0.2,0.8,0.1,0.8,0.1),nrow=3, byrow=TRUE, dimnames=list(statesNames,statesNames)))
mcC<-new("markovchain", states=statesNames, transitionMatrix=matrix(c(0.2,0.5,0.3,0,0.2,0.8,0.1,0.8,0.1),nrow=3, byrow=TRUE, dimnames=list(statesNames,statesNames))) 
mclist <- new("markovchainList", markovchains = list(mcA, mcB, mcC))   
mc <- mclist

clusterExport(cl, "mclist")

f <- function(x) {
  n <- length(mclist@markovchains)
  seq <- character(length = n)
  
  t0  <- (mclist@markovchains[[1]]@states)[1]
  
  
  for(i in 1:n) {
    stateName <- mclist@markovchains[[i]]@states
    t0 <- sample(x = stateName, size = 1, 
           prob = mclist@markovchains[[i]]@transitionMatrix[which(stateName == t0 ), ]) 
    
    seq[i] <- t0
  } 
  return(seq)
}

data <- parSapply(cl, 1:10,f)

stopCluster(cl)


struct SequenceMC : public Worker
{   
  // source
  List tmatList;
  
  // store results
  List output;
  List temp
  CharacterVector help
  
  // constructors
  sequenceMC(List tmat, int n) : tmatList(tmat), help(n) {}
  
  // process just the elements of the range I've been asked to
  void operator()(std::size_t begin, std::size_t end) {
    int n = end-begin;
    int m = help.size()
    temp = markovchainListRcpp(n, tmatList, false); 
    
    for(int i = 0;i < n;i++) {
      for(int j = 0;j < m;j++) {
        help[j] = temp[1][m*i+j];
      }  
      output.push_back(help);
    }
  }
};

// [[Rcpp::export]]
List markovchainSequenceParallelRcpp(S4 object, int n) {
  
  List tmatList = object.slot("markovchains")
  
  SequenceMC sequenceMC(tmatList, tmatList.size());
  
  // call it with parallelFor
  parallelReduce(0, n, sequenceMC);
  
  return sequenceMC.output;
}
