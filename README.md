# markovchain-gsoc-2016

## Comparing R vs C++ version of the functions

### markovchainSequence vs markovchainSequenceDeepak

```{r}
# > myMc <- as(matrix(c(.2, .8, .5, .5), byrow = TRUE, nrow = 2), "markovchain")
# > microbenchmark(markovchainSequenceDeepak(10000, myMc),markovchainSequence(10000, myMc))
# Unit: milliseconds
#                                   expr        min         lq       mean     median         uq       max neval
#  markovchainSequenceDeepak(10000, myMc)   7.409771   7.527679   7.922717   7.645513   7.848832  15.21632   100
#        markovchainSequence(10000, myMc) 180.709934 182.649039 188.409250 184.998183 187.537142 340.19028   100
```

### rmarkovchain vs rmarkovchainDeepak

```{r}
#  > myMc <- as(matrix(c(.2,.8,.5,.5), byrow = TRUE, nrow=2), "markovchain")
#  > myMc2 <- as(matrix(c(.3,.7,.1,.9), byrow = TRUE, nrow=2), "markovchain")
#  > myMc3 <- as(matrix(c(.4,.6,.3,.7), byrow = TRUE, nrow=2), "markovchain")
#  > myMcList <- new("markovchainList", markovchains = list(myMc, myMc2, myMc3), 
#  +               name = "Non - homogeneous Markov Chain")
#  > microbenchmark(rmarkovchainDeepak(n=100,object=myMcList,include.t0=TRUE), rmarkovchain(n=100,object=myMcList,include.t0=TRUE))
#  Unit: milliseconds
#                                                                expr       min        lq     mean    median        uq       max neval
#   rmarkovchainDeepak(n = 100, object = myMcList, include.t0 = TRUE)  2.689414  2.807666  3.27313  2.843361  2.902687  7.599461   100
#         rmarkovchain(n = 100, object = myMcList, include.t0 = TRUE) 12.938390 13.056885 13.73651 13.169860 13.628384 18.170999   100
```

### For both functions we can see the improvement in computational time.
