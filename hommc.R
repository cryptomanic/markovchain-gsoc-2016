# all transition matrices
# n*s*s n = order s = number of categorical sequences
# verified using two examples from research paper
allTransMat <- function(data, order = 2) {
  n <- order # order
  uelement <- sort(unique(as.character(data))) # unique element
  m <- length(uelement) # dim of trans-matrix
  s <- nrow(data) # number of categorical sequence
  lseq <- ncol(data) # length of each categorical sequence
  
  # store all transition matrices
  allTmat <- array(dim = c(length(uelement), length(uelement), n*s*s))
  
  t <- 1 # help
  
  for(i in 1:s) {
    for(j in 1:s) {
      x <- data[j, ] # jth sequence
      y <- data[i, ] # ith sequence
      
      # jumps
      for(h in 1:n) {
        # column wise
        allTmat[ , , t] <- t(createSequenceMatrix(matrix(c(x[1:(lseq-h)], y[-(1:h)]), ncol = 2, byrow = FALSE), 
                                                  toRowProbs = TRUE, possibleStates = uelement, sanitize = TRUE))
        t <- t + 1
      }
    }
  }
  return(allTmat)
}

# distribution of each categorical sequence based on the frequency
# verified using two examples from research paper
allFreqProbMat <- function(data) {
 
  uelement <- sort(unique(as.character(data))) # unique element
  m <- length(uelement) # dim of trans-matrix
  s <- nrow(data) # number of categorical sequence
  
  # frequency based probability for all sequences
  freqMat <- array(0, dim = c(m, 1, s), dimnames = list(uelement)) 
  
  for(i in 1:s) {
    idata <- data[i, ] # ith categorical sequence
    
    # populate frequency matrix
    for(j in idata) {
      freqMat[j, 1, i] <- freqMat[j, 1, i] + 1
    }
    
    # normalization
    freqMat[, , i] <- freqMat[, , i] / sum(freqMat[, , i])
  }
  
  return(freqMat)
}

# objective function to pass to solnp
fn2 <- function(params, ...) {
  
  hdata <- list(...)
  
  # calculate error
  error <- 0
  
  # number of categorical sequence
  s <- hdata$s
  
  # order
  n <- hdata$n
  
  # number of uniq states || dimension of t-matrix
  m <- hdata$m
  
  # array of transition matrices
  allTmat <- hdata$allTmat
  
  # all frequency matrix
  freqMat <- hdata$freqMat
  
  # norm
  Norm <- hdata$Norm
  
  for(i in 1:s) {
    helper <- matrix(0, nrow = m*n, ncol = 1)
    for(j in 1:s) {
      helper2 <- matrix(0, nrow = m, ncol = 1)
      y <- n * (j - 1 + s * (i - 1))
      
      for(k in 1:n) {
        helper2 <- helper2 + params[y + k] * (allTmat[ , , y + k] %*% matrix(freqMat[ , , j]))
      }
      
      helper[1:m, ] <- helper[1:m, ] + helper2
      
      if(i == j && n>= 2) {
        for(k in 2:n) {
          p <- (k - 1) * m
          helper[(p + 1):(p + m)] <- freqMat[ , , j]
        }
      }
    }
    error <- error + sum(abs((helper - freqMat[ , , i]) ^ Norm))
  }
  
  return(error ^ (1 / Norm))
}

# equality constraint function to pass to solnp
eqn2 <- function(params, ...) {
  
  hdata <- list(...)
  
  # number of categorical sequence
  s <- hdata$s
  
  # order
  n <- hdata$n
  
  toReturn <- numeric()
  
  for(i in 1:s) {
    toReturn[i] <- sum(params[((i - 1) * n * s + 1):(i * n * s)])
  }
  return(toReturn)
}

fitHighOrderMultivarMC <- function(seqMat, order = 2, Norm = 2) {
  
  # array of transition matrices
  allTmat <- allTransMat(seqMat, order = order)
  
  # array of freq probability
  freqMat <- allFreqProbMat(seqMat)
  
  n <- order # order
  uelement <- sort(unique(as.character(seqMat))) # unique element
  m <- length(uelement) # dim of trans-matrix
  s <- nrow(seqMat) # number of categorical sequence
  
  lmbda <- rep(1 / (n * s), n * s * s)
  
  fit <- solnp(pars = lmbda, fun =  fn2, eqfun = eqn2, eqB = rep(1, s), LB = rep(0, n * s * s), 
               control = list(trace = 0), allTmat = allTmat, freqMat = freqMat, n = n, m = m,
               s = s, Norm = Norm)
  
  return(fit$pars)
}

modelAccuracy <- function(seqMat, order = 2, Norm = 2) {
  n <- order # order
  uelement <- sort(unique(as.character(seqMat))) # unique element
  m <- length(uelement) # dim of trans-matrix
  s <- nrow(seqMat) # number of categorical sequence
  
  # array of transition matrices
  allTmat <- allTransMat(seqMat, order)
  params <- fitHighOrderMultivarMC(seqMat, order = order, Norm = 2)
  
  # initial probability distribution
  iprob <- array(0, dim = c(m, n, s), dimnames = list(uelement))
  iprob2 <- array(0, dim = c(m, n, s), dimnames = list(uelement))
  
  # populate initial distribution
  for(i in 1:s) {
    ini <- seqMat[i, 1] # first state
    iprob[ini, , i] <- 1
  }
  
  err <- matrix(0, s, 1) # to be returned
  
  
  len <- length(seqMat[1, ])
  
  for(t in 2:len) {
  
  for(i in 1:s) {
    helper <- matrix(0, nrow = m*n, ncol = 1)
    for(j in 1:s) {
      helper2 <- matrix(0, nrow = m, ncol = 1)
      y <- n * (j - 1 + s * (i - 1))
      
      for(k in 1:n) {
        helper2 <- helper2 + params[y + k] * (allTmat[ , , y + k] %*% matrix(iprob[ , k, j]))
      }
      
      helper[1:m, ] <- helper[1:m, ] + helper2
      
      if(i == j && n>= 2) {
        for(k in 2:n) {
          p <- (k - 1) * m
          helper[(p + 1):(p + m)] <- iprob[ , k-1, j]
        }
      }
      
    }
    
    for(j in 1:n)
      iprob2[ , j, i] <- helper[((j-1)*m+1):(j*m), ]
    
  }
  iprob <- iprob2
  
  
  for(i in 1:s) {
    d <- matrix(iprob[,,i])
    
    d  <- matrix(d, ncol = n, byrow = F)
    
    mval <- 0
    rowNo <- 0
    
    for(j in 1:m) {
      y <- d[j, ]
      
      if(max(y) > mval)  {
        mval <- max(y)
        rowNo <- j
      }
    }
    
    z <- uelement[rowNo]
    if(z == seqMat[i, t])
      err[i , 1] <- err[i , 1] + 1
  }
  
  
  }
  return(err)
  
}
