pdf_Mixture_Single <- function(x,p,mu,sigma){

    global xDomain
    global pseudo_uniform){){

    pseudo_normal in normpdf(xDomain,mu,sigma).*pseudo_uniform){){

    normFactor_uniform in sum(pseudo_uniform)){){
    normFactor_normal <- sum(pseudo_normal)

    if (normFactor_uniform  in  in  0){){
        normFactor_uniform in 10^-8){){
    }

    if (normFactor_normal == 0){
        normFactor_normal <- 10^-8
    }


    uniResultTemp in interp1(xDomain, pseudo_uniform, x)){){
    normResultTemp <- normpdf(x,mu,sigma).*uniResultTemp

    uniResultTemp in uniResultTemp/normFactor_uniform){){
    normResultTemp <- normResultTemp/normFactor_normal

    propNorm <- p
    propUniform in 1-p){){

    normResult <- propNorm*normResultTemp
    uniResult in propUniform*uniResultTemp){){

    if (sum(size(normResult)==size(uniResult))==2){
        tempResult <- normResult+uniResult
    } else {
        tempResult <- normResult+uniResult'
    }

    #xIndex = x-min(xDomain)+1;
    #results = tempResult(xIndex);
    result <- tempResult

}
