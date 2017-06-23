library(R.matlab)
library(stringr)

# This script fits a single-episode model (M1) to error data from
# an attentional blink task.
# There are two targets and the participant clicks on which letters she thinks were presented.
# The responses are aggregated together into a single serial position error histogram,
# which will be fit with the mixture of the pseudouniform and the Gaussian.
#
# You should have a folder named 'ModelOutput', where the model output will
# be saved. It will overwrite any saved output with the same name
# ('ModelOutput_ThisSample_Single.mat'). You should also have a folder
# named 'Data', where the compiled data from each sample should be (called
# something like 'CompiledData_ThisSample.mat'). See the documentation for
# information regarding the format of this data file.
#
# The objective function pdf_Mixture_Single.m should be on the MATLAB path.
#
# After also running the function AB_Model_Dual, you should run
# AB_Compare_Models to put everything together and get your final parameter
# estimates.
#
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# This section contains some things that will need to be set locally,
# depending on your directory structure and your data.

# Provide the path.
thisPath <- 'data/'

# Provide a name for each sample,
# so files can be read and written with corresponding filenames.
sampleNames <- c("Warwick","MIT","Western","Berkeley","SydneyObject","SydneyWord")

debugDoJustOneFit<-TRUE #For debugging purposes, only fit a single set of observations

# Provide some properties of the data for each sample, in order
allNParticipants <- c(20, 11, 12, 12, 32, 31)# Number of participants
allNLetters <- c(26, 26, 26, 14, 20, 20)# Number of items in a stream on each trial
allNTrials <- c(100, 100, 100, 210, 280, 280)# Total number of trials per block across all lags
allNBlocks <- c(4, 4, 4, 2, 1, 1)# Total number of blocks

# Set some model-fitting parameters.
nReplicates <- 100# Number of times to repeat each fit with different starting values
smallNonZeroNumber <- 10^-5# Useful number for when limits can't be exactly zero but can be anything larger
fitMaxIter <- 10^4# Maximum number of fit iterations
fitMaxFunEvals <- 10^4# Maximum number of model evaluations

# Set some parameter bounds. Pat apparently found these were needed to 
# prevent over-fitting to blips in the distributions. These
# values are about right in most cases, but might need some tweaking if
# e.g. you were analysing data with an unusually high or low item rate.
muBound <- 4   #will only consider -4 to +4 for mu 
sigmaBound <- 4 #will only consider 0 to 4 for sigma

# Ordinarily you wouldn't want to change these, but you might want to
# provide a different function with a different number of parameters.
nFreeParameters <- 3
source("pdf_Mixture_Single.R")
pdf_normmixture_single <- pdf_Mixture_Single

# Just for diagnostics. Setting this to 1 will show the fitted
# distributions for every participant at every lag, so it's not practical
# to run it on large datasets.
plotFits <- 0

# -------------------------------------------------------------------------
# Declare global variables that need to be accessed by the objective
# function.

# Determine number of samples
nSamples <- length(sampleNames)

# Cycle through each sample
samplesToFit<- ifelse(debugDoJustOneFit,1,nSamples)
for (thisSample in 1:samplesToFit) {

    # Load the compiled data file for this sample.
    data <- readMat(str_c('data/CompiledData_', sampleNames[thisSample], '.mat'))
    allLags <- data$allLags
    allT1Error <- data$allT1Error
    allT1Pos <- data$allT1Pos
    allT1Resp <- data$allT1Resp
    allT2Error <- data$allT2Error
    allT2Pos <- data$allT2Pos
    allT2Resp <- data$allT2Resp

    # Extract the relevant task parameters specified above.
    nLetters <- allNLetters[thisSample]
    nTrials <- allNTrials[thisSample]
    nBlocks <- allNBlocks[thisSample]
    nParticipants <- allNParticipants[thisSample]

    # Work out the number of lags.
    listLags <- unique(allLags)
    listLags[is.nan(listLags)] <- c()
    nLags <- length(listLags)

    # Work out possible positions in the stream for T1.
    listT1Pos <- unique(allT1Pos)
    nT1Pos <- length(listT1Pos)

    # Get the list of T1 errors and extract some properties.
    listT1Errors <- unique(allT1Error)# List of unique observed T1 errors
    listT1Errors[is.nan(listT1Errors)] <- c()# Get rid of NaN values
    nT1Errors <- length(listT1Errors)# Number of unique observed T1 errors
    minT1Error <- min(listT1Errors)# Lowest (most negative) T1 error
    maxT1Error <- max(listT1Errors)# Highest (most positive) T1 error

    nTrialsPerLag <- (nTrials*nBlocks)/nLags# Calculate number of trials per lag

    # The following output to the command window is just to keep track of
    # what the program is doing.
    cat(sprintf('\n\n%s\n\n', str_to_upper(sampleNames[thisSample])))

    # Build empty matrices for storing parameter estimates for each
    # participant, at each lag. Also build matrices to store upper and
    # lower bounds, and minimum negative log likelihoods.

    # Note that this program will also fit the single-episode model (M1) to
    # T2, even though that won't be used at any point in the analysis
    # proper. However, it can be worthwhile running it as a sanity check.

    allT1Estimates_byParticipant <- array(data=NA, dim=c(nParticipants,nLags,nFreeParameters))
    allT1LowerBounds_byParticipant <- array(data=NA, dim=c(nParticipants,nLags,nFreeParameters))
    allT1UpperBounds_byParticipant <- array(data=NA, dim=c(nParticipants,nLags,nFreeParameters))
    allT1MinNegLogLikelihoods_byParticipant <- array(data=NA, dim=c(nParticipants,nLags))

    allT2Estimates_byParticipant <- array(data=NA, dim=c(nParticipants,nLags,nFreeParameters))
    allT2LowerBounds_byParticipant <- array(data=NA, dim=c(nParticipants,nLags,nFreeParameters))
    allT2UpperBounds_byParticipant <- array(data=NA, dim=c(nParticipants,nLags,nFreeParameters))
    allT2MinNegLogLikelihoods_byParticipant <- array(data=NA, dim=c(nParticipants,nLags))

    allT1T2Estimates_byParticipant <- array(data=NA, dim=c(nParticipants,nLags,nFreeParameters))
    allT1T2LowerBounds_byParticipant <- array(data=NA, dim=c(nParticipants,nLags,nFreeParameters))
    allT1T2UpperBounds_byParticipant <- array(data=NA, dim=c(nParticipants,nLags,nFreeParameters))
    allT1T2MinNegLogLikelihoods_byParticipant <- array(data=NA, dim=c(nParticipants,nLags))

    # Cycle through each participant.
    nParticipantsToFit<- ifelse(debugDoJustOneFit,1,nParticipants)
    for (thisParticipant in 1:nParticipantsToFit) {

        # Keep track of progress in the command window.
        print(sprintf('\nParticipant %d... ', thisParticipant))

        # Extract the relevant lists of T1 and T2 errors, T1 and T2 stream
        # positions, and corresponding lags.
        T1Error <- allT1Error[thisParticipant,,]
        T2Error <- allT2Error[thisParticipant,,]
        T1Pos <- allT1Pos[thisParticipant,,]
        T2Pos <- allT2Pos[thisParticipant,,]
        Lags <- allLags[thisParticipant,,]

        # Cycle through each lag.
        nLagsToFit<- ifelse(debugDoJustOneFit,1,nLags)
        for (thisLag in 1:nLagsToFit){

            # Keep track of progress in the command window.
            print(sprintf('L%d ', thisLag))

            # Find the lag value for this numbered lag. These are usually
            # the same, but this allows for the possibility that some lag
            # values aren't tested.
            thisLagVal <- listLags[thisLag]

            # Go through and replace NaN values. Here, we assume that a NaN
            # value is a guess, and replace it with a random sample from
            # the range of possible errors. You might want to think about
            # this if you have data with lots of missing values, but it
            # makes virtually no difference in these data.

            # Identify T1 trials at this lag.
            hasThisLag <- Lags==thisLagVal
            theseT1Error <- t(T1Error[hasThisLag])
            theseT1Pos <- t(T1Pos[hasThisLag])

            if (sum(is.nan(theseT1Error)) > 0){# If there is at least one NaN

                # Find NaNs.
                replaceCells <- which(is.nan(theseT1Error))

                # Find the T1 positions of the trials to be replaced.
                replacePositions <- theseT1Pos[replaceCells]

                # Cycle through each point to be replaced.
                for (thisReplace in 1:length(replaceCells)){

                    # Replace with a random possible value.
                    theseT1Error[replaceCells[thisReplace]] <- sample(1:nLetters, 1)-replacePositions[thisReplace]

                }

            }

            # Identify T2 trials at this lag.
            theseT2Error <- t(T2Error[hasThisLag])
            theseT2Pos <- t(T2Pos[hasThisLag])

            if (sum(is.nan(theseT2Error)) > 0){# If there is at least one NaN

                # Find NaNs.
                replaceCells <- which(is.nan(theseT2Error))

                # Find the T2 positions of the trials to be replaced.
                replacePositions <- theseT2Pos[replaceCells]

                # Cycle through each point to be replaced.
                for (thisReplace in 1:length(replaceCells)){

                    # Replace with a random possible value.
                    theseT2Error[replaceCells[thisReplace]] <- sample(nLetters, 1)-replacePositions[thisReplace]

                }

            }

            # Combine T1 and T2 distributions on a common scale (i.e. alter
            # T1 errors to reflect position relative to T2.
            theseT1T2Error <- c(theseT1Error-listLags[thisLag], theseT2Error)

            # Get minimum and maximum error values.
            minT2Error <- min(theseT2Error)
            maxT2Error <- max(theseT2Error)
            minT1T2Error <- min(theseT1T2Error)
            maxT1T2Error <- max(theseT1T2Error)

            # Unpack the parameter bounds to feed into the model fitting.
            # These are set near the top of the script, but need to be
            # unpacked for each scenario.

            # Unpack mean (latency) bounds for T1.
            mu_lb_T1 <- -muBound
            mu_ub_T1 <- muBound

            # Unpack SD (precision) bounds for T1.
            sigma_lb_T1 <- smallNonZeroNumber
            sigma_ub_T1 <- sigmaBound

            # Unpack mean (latency) bounds for T2.
            mu_lb_T2 <- -muBound
            mu_ub_T2 <- sigmaBound

            # Unpack SD (precision) bounds for T2.
            sigma_lb_T2 <- smallNonZeroNumber
            sigma_ub_T2 <- sigmaBound

            # Unpack mean (latency) bounds for compined T1 & T2.
            mu_lb_T1T2 <- -muBound-listLags[thisLag]
            mu_ub_T1T2 <- muBound

            # Unpack SD (precision) bounds for combined T1 & T2.
            sigma_lb_T1T2 <- smallNonZeroNumber
            sigma_ub_T1T2 <- sigmaBound

            # Fit the model to the T1 distribution.

            # Keep track of the minimum negative log likelihood on each
            # replicate. Start at infinity so the first replicate
            # automatically qualifies as the best candidate up to that
            # point.
            minNegLogLikelihood <- Inf

            # Calculate the domain of possible errors (xDomain).
            xPosition <- unique(theseT1Pos)
            minX_T1 <- min(xPosition)
            maxX_T1 <- max(xPosition)
            minErr <- 1-maxX_T1-1
            maxErr <- nLetters-minX_T1+1
            xDomain <- minErr:maxErr

            # Generate the 'pseudo-uniform' distribution, which is the
            # expected distribution of errors if a random guess was
            # provided on every trial. This isn't an actual uniform
            # distribution because the most extreme errors are only
            # possible on trials in which targets appear at their most
            # extreme positions.
            pseudo_uniform <- rep(0, length(xDomain))

            # Cycle through each possible T1 position.
            for (thisPosNo in 1:length(theseT1Pos)){

                # Identify the actual T1 position corresponding to the
                # position number. For example, the first position number
                # might be the 7th position in the stream.
                thisPos <- theseT1Pos[thisPosNo]

                # Add to the pseudo-uniform distribution one unit for every
                # possible error given that T1 position.
                pseudo_uniform[(1-thisPos-minErr+1):(nLetters-thisPos-minErr+1)] = pseudo_uniform[(1-thisPos-minErr+1):(nLetters-thisPos-minErr+1)]+rep(1,nLetters)
            }

            # Cycle through a number of replicates of the fitting
            # procedure with different randomised starting values across
            # the range dictated by the bounds.
            nReplicatesToFit<- ifelse(debugDoJustOneFit,1,nReplicates)
            for (thisReplicate in 1:nReplicatesToFit) {

                # Randomise starting values for each parameter.
                pGuess <- max(c(smallNonZeroNumber, runif(1)))
                muGuess <- (2*muBound*runif(1))-muBound
                sigmaGuess <- sigmaBound*runif(1)+smallNonZeroNumber

                # Compile to feed into the MLE function.
                parameterGuess <- c(pGuess, muGuess, sigmaGuess)
                parameterLowerBound <- c(smallNonZeroNumber, mu_lb_T1, sigma_lb_T1)
                parameterUpperBound <- c(1, mu_ub_T1, sigma_ub_T1)

                # Ensure guesses satisfy bounds, and round them marginally
                # up or down if necessary.
                for (i in 1:length(parameterGuess)) {
                    if (parameterGuess[i] < parameterLowerBound[i])
                        parameterGuess[i] <- parameterLowerBound[i]

                    if (parameterGuess[i] > parameterUpperBound[i])
                        parameterGuess[i] <- parameterUpperBound[i]
                }

                # Run the MLE function.
                # [currentEstimates, currentCIs] <- mle(theseT1Error, 'pdf', pdf_normmixture_single, 'start', parameterGuess, 'lower', parameterLowerBound, 'upper', parameterUpperBound, 'options', options)
                pdf_normmixture_single_par <- function(par)
                {
                    p <- par[1]
                    mu <- par[2]
                    sigma <- par[3]
                    result <- pdf_normmixture_single(theseT1Error, p, mu, sigma)
                    #cat("p ", p, " mu ", mu, " sigma ", sigma, " result ", result, "\n")
                    return(-exp(sum(log(result))))
                }                
                cat("parameterGuess", parameterGuess, "\n")
                fit <- optim(parameterGuess, pdf_normmixture_single_par, lower=parameterLowerBound, upper=parameterUpperBound, control=list(trace=6), method="L-BFGS-B")
                currentEstimates <- fit$par
                cat("currentEstimates=", currentEstimates, "\n")
                # Compute the negative log likelihood of the fitted model.
                thisNegLogLikelihood <- -sum(log(pdf_normmixture_single(theseT1Error,currentEstimates[1],currentEstimates[2],currentEstimates[3])))

                # Check whether this is lower than the lowest so far.
                if (minNegLogLikelihood > thisNegLogLikelihood){

                    # If so, store this as the current best estimate.
                    minNegLogLikelihood <- thisNegLogLikelihood
                    bestEstimates <- currentEstimates
                    # bestEstimateCIs <- currentCIs

                }

            }

            # Enter the best estimates into the parameter matrices.
            allT1Estimates_byParticipant[thisParticipant,thisLag,] <- bestEstimates
            # allT1LowerBounds_byParticipant[thisParticipant,thisLag,] <- bestEstimateCIs[1,]
            # allT1UpperBounds_byParticipant[thisParticipant,thisLag,] <- bestEstimateCIs[2,]
            allT1MinNegLogLikelihoods_byParticipant[thisParticipant,thisLag] <- minNegLogLikelihood
            STOP
            # Plot the distributions if required.
            # if (plotFits){
            #     pFig <- figure('Color','white','Name',['Participant ' num2str(thisParticipant) ', Lag ' num2str(thisLag)])##ok<UNRCH>
            #     subplot(1,3,1)
            #     tBars <- hist(theseT1Error, xDomain)
            #     tBars <- tBars/sum(tBars)
            #     bar(xDomain,tBars)
            #     hold on
            #     yModelPlot <- pdf_normmixture_single(xDomain, bestEstimates(1), bestEstimates(2), bestEstimates(3))
            #     plot(xDomain,yModelPlot,'r--')
            #     axis square
            #     axis([minErr-1 maxErr+1 0 0.75])
            #     title('T1')
            #     drawnow
            # }

            # Fit the model to the T2 distribution.

            # Keep track of the minimum negative log likelihood on each
            # replicate.
            minNegLogLikelihood <- Inf

            # Calculate the domain of possible errors (xDomain).
            xPosition <- unique(theseT2Pos)
            minX_T2 <- min(xPosition)
            maxX_T2 <- max(xPosition)
            minErr <- 1-maxX_T2-1
            maxErr <- nLetters-minX_T2+1
            xDomain <- minErr:maxErr

            # Generate the 'pseudo-uniform' distribution, which is the
            # expected distribution of errors if a random guess was
            # provided on every trial.
            pseudo_uniform <- rep(0, length(xDomain))

            # Cycle through each possible T2 position.
            for (thisPosNo in 1:length(theseT2Pos)){

                # Identify the actual T2 position corresponding to the
                # position number. For example, the first position number
                # might be the 7th position in the stream.
                thisPos <- theseT2Pos[thisPosNo]

                # Add to the pseudo-uniform distribution one unit for every
                # possible error given that T2 position.
                pseudo_uniform[(1-thisPos-minErr+1):(nLetters-thisPos-minErr+1)] = pseudo_uniform[(1-thisPos-minErr+1):(nLetters-thisPos-minErr+1)]+rep(1,nLetters)

            }

            # Cycle through a number of replicates of the fitting
            # procedure with different randomised starting values across
            # the range dictated by the bounds.

            for (thisReplicate in 1:nReplicates){

                # Randomise starting values for each parameter.
                pGuess <- max(c(smallNonZeroNumber, runif(1)))
                muGuess <- (2*muBound*runif(1))-muBound
                sigmaGuess <- sigmaBound*runif(1)+smallNonZeroNumber

                # Compile to feed into the MLE function.
                parameterGuess <- c(pGuess, muGuess, sigmaGuess)
                parameterLowerBound <- c(smallNonZeroNumber, mu_lb_T2, sigma_lb_T2)
                parameterUpperBound <- c(1, mu_ub_T2, sigma_ub_T2)

                # Ensure guesses satisfy bounds.
                parameterGuess <- max(parameterGuess,parameterLowerBound)
                parameterGuess <- min(parameterGuess,parameterUpperBound)

                # Run the MLE function.
                # [currentEstimates, currentCIs] <- mle(theseT2Error, 'pdf', pdf_normmixture_single, 'start', parameterGuess, 'lower', parameterLowerBound, 'upper', parameterUpperBound, 'options', options)
                pdf_normmixture_single_par <- function(par)
                {
                    p <- par[1]
                    mu <- par[2]
                    sigma <- par[3]
                    result <- pdf_normmixture_single(theseT1Error, p, mu, sigma)
                    cat("p ", p, " mu ", mu, " sigma ", sigma, " result ", result, "\n")
                    return(-exp(sum(log(result))))
                }                
                fit <- optim(parameterGuess, pdf_normmixture_single_par, lower=parameterLowerBound, upper=parameterUpperBound, control=list(trace=6), method="L-BFGS-B")
                currentEstimates <- fit$par

                # Compute the negative log likelihood of the fitted model.
                thisNegLogLikelihood <- sum(log(pdf_normmixture_single(theseT2Error,currentEstimates[1],currentEstimates[2],currentEstimates[3])))

                # Check whether this is lower than the lowest so far.
                if (minNegLogLikelihood > thisNegLogLikelihood){

                    # If so, store this as the current best estimate.
                    minNegLogLikelihood <- thisNegLogLikelihood
                    bestEstimates <- currentEstimates
                    # bestEstimateCIs <- currentCIs

                }

            }

            # Enter the best estimates into the parameter matrices.
            allT2Estimates_byParticipant[thisParticipant,thisLag,] <- bestEstimates
            # allT2LowerBounds_byParticipant[thisParticipant,thisLag,] <- bestEstimateCIs[1,]
            # allT2UpperBounds_byParticipant[thisParticipant,thisLag,] <- bestEstimateCIs[2,]
            allT2MinNegLogLikelihoods_byParticipant[thisParticipant,thisLag] <- minNegLogLikelihood

            # Plot the distributions if required.
            # if (plotFits){
            #     subplot(1,3,2)##ok<UNRCH>
            #     tBars <- hist(theseT2Error, xDomain)
            #     tBars <- tBars/sum(tBars)
            #     bar(xDomain,tBars)
            #     hold on
            #     yModelPlot <- pdf_normmixture_single(xDomain, bestEstimates(1), bestEstimates(2), bestEstimates(3))
            #     plot(xDomain,yModelPlot,'r--')
            #     axis square
            #     axis([minErr-1 maxErr+1 0 0.75])
            #     title('T2')
            #     drawnow
            # }


            # Fit the model to the combined T1 + T2 distribution. First,
            # combine the distributions on a common scale (i.e. alter T1
            # errors to reflect position relative to T2.

            # Keep track of the minimum negative log likelihood on each
            # replicate.
            minNegLogLikelihood <- Inf

            # Calculate the domain of possible errors (xDomain).
            theseT1T2Pos <- theseT2Pos# Errors are relative to T2
            xPosition <- unique(theseT1T2Pos)
            minX <- min(xPosition)
            maxX <- max(xPosition)
            minErr <- 1-maxX-1
            maxErr <- nLetters-minX+1
            xDomain <- minErr:maxErr

            # Generate the 'pseudo-uniform' distribution, which is the
            # expected distribution of errors if a random guess was
            # provided on every trial.
            pseudo_uniform = rep(0, length(xDomain))

            # Cycle through each possible T2 position.
            for (thisPosNo in 1:length(theseT1T2Pos)){

                # Identify the actual T2 position corresponding to the
                # position number. For example, the first position number
                # might be the 7th position in the stream.
                thisPos <- theseT1T2Pos[thisPosNo]

                # Add to the pseudo-uniform distribution one unit for every
                # possible error given that T2 position.
                pseudo_uniform[(1-thisPos-minErr+1):(nLetters-thisPos-minErr+1)] = pseudo_uniform[(1-thisPos-minErr+1):(nLetters-thisPos-minErr+1)]+rep(1,nLetters)

            }

            # Cycle through a number of replicates of the fitting
            # procedure with different randomised starting values across
            # the range dictated by the bounds.

            for (thisReplicate in 1:nReplicates){

                # Randomise starting values for each parameter.
                pGuess <- max(c(smallNonZeroNumber, runif(1)))
                muGuess <- (2*muBound*runif(1))-muBound
                sigmaGuess <- sigmaBound*runif(1)+smallNonZeroNumber

                # Compile to feed into the MLE function.
                parameterGuess <- c(pGuess, muGuess, sigmaGuess)
                parameterLowerBound <- c(smallNonZeroNumber, mu_lb_T1T2, sigma_lb_T1T2)
                parameterUpperBound <- c(1, mu_ub_T1T2, sigma_ub_T1T2)

                # Ensure guesses satisfy bounds.
                parameterGuess <- max(c(parameterGuess,parameterLowerBound))
                parameterGuess <- min(c(parameterGuess,parameterUpperBound))

                # Run the MLE function.
                # [currentEstimates, currentCIs] <- mle(theseT1T2Error, 'pdf', pdf_normmixture_single, 'start', parameterGuess, 'lower', parameterLowerBound, 'upper', parameterUpperBound, 'options', options)
                pdf_normmixture_single_par <- function(par)
                {
                    p <- par[1]
                    mu <- par[2]
                    sigma <- par[3]
                    result <- pdf_normmixture_single(theseT1Error, p, mu, sigma)
                    cat("p ", p, " mu ", mu, " sigma ", sigma, " result ", result, "\n")
                    return(-exp(sum(log(result))))
                }                
                fit <- optim(parameterGuess, pdf_normmixture_single_par, lower=parameterLowerBound, upper=parameterUpperBound, control=list(trace=6), method="L-BFGS-B")
                currentEstimates <- fit$par

                # Compute the negative log likelihood of the fitted model.
                thisNegLogLikelihood <- -sum(log(pdf_normmixture_single(theseT1T2Error,currentEstimates[1],currentEstimates[2],currentEstimates[3])))

                # Check whether this is lower than the lowest so far.
                if (minNegLogLikelihood > thisNegLogLikelihood){

                    # If so, store this as the current best estimate.
                    minNegLogLikelihood <- thisNegLogLikelihood
                    bestEstimates <- currentEstimates
                    # bestEstimateCIs <- currentCIs

                }

            }

            # Enter the best estimates into the parameter matrices.
            allT1T2Estimates_byParticipant[thisParticipant,thisLag,] <- bestEstimates
            # allT1T2LowerBounds_byParticipant[thisParticipant,thisLag,] <- bestEstimateCIs[1,]
            # allT1T2UpperBounds_byParticipant[thisParticipant,thisLag,] <- bestEstimateCIs[2,]
            allT1T2MinNegLogLikelihoods_byParticipant[thisParticipant,thisLag] <- minNegLogLikelihood

            # Plot the distributions if required.
            # if (plotFits){
            #     subplot(1,3,3)##ok<UNRCH>
            #     tBars <- hist(theseT1T2Error, xDomain)
            #     tBars <- tBars/sum(tBars)
            #     bar(xDomain,tBars)
            #     hold on
            #     yModelPlot <- pdf_normmixture_single(xDomain, bestEstimates[1], bestEstimates[2], bestEstimates[3])
            #     plot(xDomain,yModelPlot,'r--')
            #     axis square
            #     axis([minErr-1 maxErr+1 0 0.75])
            #     title('T1+T2')
            #     drawnow
            # }

        }

    }

    # Keep track of progress in the command window.
    cat('\n\n')

    # Change directory to store model output.
    # cd([thisPath 'ModelOutput'])

    # Save model output.
    # save(['ModelOutput_' sampleNames{thisSample} '_Single.mat'])

}

# Turn this warning back on.
# warning('on', 'stats:mlecov:NonPosDefHessian')
