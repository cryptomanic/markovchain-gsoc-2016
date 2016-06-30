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


# A <- "66662626226266262444566122666262662622621226662126266226222626222226226666122622226222233232666626266262662662234331312161661662626222661626121626222266166226222344464616616666162226666266226262226222666632262222226262226226626662223334166166161666616662122222236666626"
# 
# B <- "16616111111666121661116621661116121622222616612166611166111161121616116262666366166222322666116266262661366111223226222161611621112216111126111161612161661612222332226666211611161616161166211661126266612616111161611661661616611662222222226666166616616611613335166666666"
# 
# C <- "66666662666666626666266622666666616266666666266126166162666666626662661666666633632122166161666666166616666666666626666666622662612666266266261626212662262622626662226626622612126622661221626221156361661226162661626266616166222123616161616661166666166616116666666616616"
# 
# D <- "62222334445433626663443333326634444342622622663454463666262662264454343446266226266266266262635554443626626262262662644444463662626262662222222223335545333626622622226232236322344445544662626222222255445526266262622334454443436262222222222234444544432226222626262222232"
# 
# E <- "62222334445433626623443443322634444342322633663454533266262662264444445446266226266266266262634444444626626266662622644444463362226262222222222223645555246626622622226232236322344445543362622263222255444436266262622334454444436262226222222234444544432226662626262222222"
# 
# data <- matrix(nrow = 5, ncol = 269)
# data[1, ] <- strsplit(A, "")[[1]]
# data[2, ] <- strsplit(B, "")[[1]]
# data[3, ] <- strsplit(C, "")[[1]]
# data[4, ] <- strsplit(D, "")[[1]]
# data[5, ] <- strsplit(E, "")[[1]]
