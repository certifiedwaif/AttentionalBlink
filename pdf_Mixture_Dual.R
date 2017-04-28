pdf_Mixture_Dual <- function(x,p1,mu1,sigma1,p2,mu2,sigma2){
    pseudo_normal1 <- dnorm(xDomain,mu1,sigma1)*pseudo_uniform
    pseudo_normal2 <- dnorm(xDomain,mu2,sigma2)*pseudo_uniform

    normFactor_uniform <- sum(pseudo_uniform)
    normFactor_normal1 <- sum(pseudo_normal1)
    normFactor_normal2 <- sum(pseudo_normal2)

    if (normFactor_uniform  ==  0) {
        normFactor_uniform = 10^-8
    }

    if (normFactor_normal1 == 0){
        normFactor_normal1 <- 10^-8
    }

    if (normFactor_normal2 == 0){
        normFactor_normal2 <- 10^-8
    }

    uniResultTemp <- interp1(xDomain, pseudo_uniform, x)
    normResult1Temp <- dnorm(x,mu1,sigma1)*uniResultTemp
    normResult2Temp <- dnorm(x,mu2,sigma2)*uniResultTemp

    uniResultTemp <- uniResultTemp/normFactor_uniform
    normResultTemp1 <- normResult1Temp/normFactor_normal1
    normResultTemp2 <- normResult2Temp/normFactor_normal2

    propNorm1 <- p1*(1-p2)
    propNorm2 <- p1*p2
    propUniform <- 1-p1

    normResult1 <- propNorm1*normResultTemp1
    normResult2 <- propNorm2*normResultTemp2
    uniResult <- propUniform*uniResultTemp

    if (sum(length(normResult1)==length(uniResult))==2){
        tempResult <- normResult1+normResult2+uniResult
    } else {
        tempResult <- normResult1+normResult2+t(uniResult)
    }

    return(tempResult)
}
