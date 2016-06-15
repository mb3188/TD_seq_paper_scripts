### R function to compute sere statistics
### from Schulze et al, BMC Genomics 2012, 13:524

#within and between condition variability
sigfun_Pearson <- function(observed) {
  #calculate lambda and expected values
  laneTotals<- colSums(observed);
  total <- sum(laneTotals)
  fullObserved <- observed[rowSums(observed)>0,];
  fullLambda <- rowSums(fullObserved)/total;
  fullLhat <- fullLambda > 0;
  fullExpected<- outer(fullLambda, laneTotals);

  #keep values
  fullKeep <- which(fullExpected > 0);
  
  #calculate degrees of freedom (nrow*(ncol -1) >> number of parameters - calculated (just lamda is calculated >> thats why minus 1)
  #calculate pearson and deviance for all values
  oeFull <- (fullObserved[fullKeep] - fullExpected[fullKeep])^2/ fullExpected[fullKeep] # pearson chisq test
  dfFull <- length(fullKeep) - sum(fullLhat!=0);
  
  return(c(sqrt(sum(oeFull)/dfFull)));
}
