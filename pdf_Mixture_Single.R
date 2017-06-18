library(signal) # Provides interp1 function
# interp1 <- function(x, v, xq)
# {
#   cat("x ", x, " v ", v, " xq ", xq, "\n")
#   x2 <- x[order(x)]
#   v2 <- v[order(x)]
#   
#   vq <- rep(0, length(xq))
#   for (idx in 1:length(xq)) {
#     i <- 1
#     while (xq[idx] >= x2[i]) i = i + 1
# 
#     # Now x[i - 1] < xq and x[i] >= xq
#     m <- (v2[i] - v2[i - 1]) / (x2[i] - x2[i - 1])
#     vq[idx] <- v2[i - 1] + m * (xq[idx] - x2[i - 1])
#   }
# 
#   return(vq)
# }

pdf_Mixture_Single <- function(x,p,mu,sigma){

    pseudo_normal <- dnorm(xDomain,mu,sigma)*pseudo_uniform
    cat("pseudo_normal ", pseudo_normal, "\n")  
    
    normFactor_uniform <- sum(pseudo_uniform)
    cat("normFactor_uniform", normFactor_uniform, "\n")  
    normFactor_normal <- sum(pseudo_normal)
    cat("normFactor_normal", normFactor_normal, "\n")  
    
    if (normFactor_uniform  == 0 || is.nan(normFactor_uniform)){
        normFactor_uniform <- 10^-8
    }

    if (normFactor_normal == 0 || is.nan(normFactor_normal)){
        normFactor_normal <- 10^-8
    }


    uniResultTemp <- interp1(xDomain, pseudo_uniform, x)
    #cat("uniResultTemp ", uniResultTemp, "\n")
    normResultTemp <- dnorm(x,mu,sigma)*uniResultTemp
    #cat("normResultTemp ", normResultTemp, "\n")
    
    uniResultTemp <- uniResultTemp/normFactor_uniform
    #cat("uniResultTemp ", uniResultTemp, "\n")
    normResultTemp <- normResultTemp/normFactor_normal
    #cat("normResultTemp ", normResultTemp, "\n")
    
    propNorm <- p
    propUniform <- 1-p

    normResult <- propNorm*normResultTemp
    uniResult <- propUniform*uniResultTemp

    if (sum(length(normResult)==length(uniResult))==2){
        tempResult <- normResult+uniResult
    } else {
        tempResult <- normResult+t(uniResult)
    }

    #xIndex = x-min(xDomain)+1;
    #results = tempResult(xIndex);
    #cat("tempResult ", tempResult, "\n")
    #tempResult <- tempResult[1]
    return(tempResult)
}
