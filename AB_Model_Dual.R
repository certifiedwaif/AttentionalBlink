# This script fits a dual-episode model (M2) to error data from
# an attentional blink task.
#
# You should have a folder named 'ModelOutput', where the model output will
# be saved. It will overwrite any saved output with the same name
# ('ModelOutput_ThisSample_Dual.mat'). You should also have a folder named
# 'Data', where the compiled data from each sample should be (called
# something like 'CompiledData_ThisSample.mat'). See the documentation for
# information regarding the format of this data file.
#
# The objective function pdf_Mixture_Single.m should be on the MATLAB path.
#
# After also running the function AB_Model_Single, you should run
# AB_Compare_Models to put everything together and get your final parameter
# estimates.
#
# The script requires the Statistics Toolbox.

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# This section contains some things that will need to be set locally,
# depending on your directory structure and your data.

# Provide the path.
thisPath <- '/Users/experimentalmode/Documents/MATLAB/ABFinal/'

# Provide a name for each sample,
# so files can be read and written with corresponding filenames.
sampleNames <- {'Warwick','MIT','Western','Berkeley','SydneyObject','SydneyWord'}

# Provide some properties of the data for each sample, in order
allNParticipants <- [20 11 12 12 32 31]# Number of participants
allNLetters <- [26 26 26 14 20 20]# Number of items in a stream on each trial
allNTrials <- [100 100 100 210 280 280]# Total number of trials per block across all lags
allNBlocks <- [4 4 4 2 1 1]# Total number of blocks

# Set some model-fitting parameters.
nReplicates <- 100# Number of times to repeat each fit with different starting values
smallNonZeroNumber <- 10^-5# Useful number for when limits can't be exactly zero but can be anything larger
fitMaxIter <- 10^4# Maximum number of fit iterations
fitMaxFunEvals <- 10^4# Maximum number of model evaluations

# Set some parameter bounds. You want these large enough that they span the
# full reasonable range of mean (latency) and SD (precision), but small
# enough to prevent over-fitting to blips in the distributions. These
# values are about right in most cases, but might need some tweaking if
# e.g. you were analysing data with an unusually high or low item rate.
muBound <- 4
sigmaBound <- 4

# Ordinarily you wouldn't want to change these, but you might want to
# provide a different function with a different number of parameters.
nFreeParameters <- 6
pdf_normmixture <- @pdf_Mixture_Dual

# Just for diagnostics. Setting this to 1 will show the fitted
# distributions for every participant at every lag, so it's not practical
# to run it on large datasets.
plotFits <- 0

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

# Declare global variables that need to be accessed by the objective
# function.
global xDomain
global pseudo_uniform){){

# Add folders to the MATLAB path.
addpath(genpath(thisPath))
cd([thisPath 'Data'])

# Turn off this warning to prevent it from filling the command
# window. This is only a problem if you get it on most
# replicates of a fit.
warning('off', 'stats:mlecov:NonPosDefHessian')

# Determine number of samples
nSamples <- numel(sampleNames)

# Cycle through each sample
for (thisSample in 1:nSamples){

    # Load the compiled data file for this sample.
    load(['CompiledData_' sampleNames{thisSample} '.mat'])

    # Extract the relevant task parameters specified above.
    nLetters <- allNLetters(thisSample)
    nTrials <- allNTrials(thisSample)
    nBlocks <- allNBlocks(thisSample)
    nParticipants <- allNParticipants(thisSample)

    # Work out the number of lags.
    listLags <- unique(allLags(:))
    listLags(isnan(listLags)) <- []
    nLags <- numel(listLags)

    # Work out possible positions in the stream for T1.
    listT1Pos <- unique(allT1Pos(:))
    nT1Pos <- numel(listT1Pos)

    # Get the list of T1 errors and extract some properties.
    listT1Errors <- unique(allT1Error(:))# List of unique observed T1 errors
    listT1Errors(isnan(listT1Errors)) <- []# Get rid of NaN values
    nT1Errors <- numel(listT1Errors)# Number of unique observed T1 errors
    minT1Error <- min(listT1Errors)# Lowest (most negative) T1 error
    maxT1Error <- max(listT1Errors)# Highest (most positive) T1 error

    nTrialsPerLag <- (nTrials*nBlocks)/nLags# Calculate number of trials per lag

    # The following output to the command window is just to keep track of
    # what the program is doing.
    fprintf('\n\n%s\n\n', upper(sampleNames{thisSample}))

    # Build empty matrices for storing parameter estimates for each
    # participant, at each lag. Also build matrices to store upper and
    # lower bounds, and minimum negative log likelihoods.

    # Note that this program will also fit the dual-episode model (M2) to
    # T1, even though that won't be used at any point in the analysis
    # proper. However, it can be worthwhile running it as an additional
    # check that T2-selected items don't turn up in the T1 distribution.

    allT1Estimates_byParticipant <- NaN(nParticipants,nLags,nFreeParameters)
    allT1LowerBounds_byParticipant <- NaN(nParticipants,nLags,nFreeParameters)
    allT1UpperBounds_byParticipant <- NaN(nParticipants,nLags,nFreeParameters)
    allT1MinNegLogLikelihoods_byParticipant <- NaN(nParticipants,nLags)

    allT2Estimates_byParticipant <- NaN(nParticipants,nLags,nFreeParameters)
    allT2LowerBounds_byParticipant <- NaN(nParticipants,nLags,nFreeParameters)
    allT2UpperBounds_byParticipant <- NaN(nParticipants,nLags,nFreeParameters)
    allT2MinNegLogLikelihoods_byParticipant <- NaN(nParticipants,nLags)

    allT1T2Estimates_byParticipant <- NaN(nParticipants,nLags,nFreeParameters)
    allT1T2LowerBounds_byParticipant <- NaN(nParticipants,nLags,nFreeParameters)
    allT1T2UpperBounds_byParticipant <- NaN(nParticipants,nLags,nFreeParameters)
    allT1T2MinNegLogLikelihoods_byParticipant <- NaN(nParticipants,nLags)

    # Set fit options.
    options <- statset('MaxIter', fitMaxIter, 'MaxFunEvals', fitMaxFunEvals, 'Display', 'off')

    # Cycle through each participant.
    for (thisParticipant in 1:nParticipants){

        # Keep track of progress in the command window.
        fprintf('\nParticipant %d... ', thisParticipant)

        # Extract the relevant lists of T1 and T2 errors, T1 and T2 stream
        # positions, and corresponding lags.
        T1Error <- squeeze(allT1Error(thisParticipant,:,:))
        T2Error <- squeeze(allT2Error(thisParticipant,:,:))
        T1Pos <- squeeze(allT1Pos(thisParticipant,:,:))
        T2Pos <- squeeze(allT2Pos(thisParticipant,:,:))
        Lags <- squeeze(allLags(thisParticipant,:,:))

        # Cycle through each lag.
        for (thisLag in 1:nLags){

            # Keep track of progress in the command window.
            fprintf('L%d ', thisLag)

            # Find the lag value for this numbered lag. These are usually
            # the same, but this allows for the possibility that some lag
            # values aren't tested.
            thisLagVal <- listLags(thisLag)

            # Go through and replace NaN values. Here, we assume that a NaN
            # value is a guess, and replace it with a random sample from
            # the range of possible errors. You might want to think about
            # this if you have data with lots of missing values, but it
            # makes virtually no difference in these data.

            # Identify T1 trials at this lag.
            hasThisLag <- Lags==thisLagVal
            theseT1Error <- T1Error(hasThisLag)'
            theseT1Pos <- T1Pos(hasThisLag)'

            if (sum(isnan(theseT1Error)) > 0){# If there is at least one NaN

                # Find NaNs.
                replaceCells <- find(isnan(theseT1Error))

                # Find the T1 positions of the trials to be replaced.
                replacePositions <- theseT1Pos(replaceCells)

                # Cycle through each point to be replaced.
                for (thisReplace in 1:numel(replaceCells)){

                    # Replace with a random possible value.
                    theseT1Error(replaceCells(thisReplace)) <- randi(nLetters)-replacePositions(thisReplace)

                }

            }

            # Identify T2 trials at this lag.
            theseT2Error <- T2Error(hasThisLag)'
            theseT2Pos <- T2Pos(hasThisLag)'

            if (sum(isnan(theseT2Error)) > 0){# If there is at least one NaN

                # Find NaNs.
                replaceCells <- find(isnan(theseT2Error))

                # Find the T2 positions of the trials to be replaced.
                replacePositions <- theseT2Pos(replaceCells)

                # Cycle through each point to be replaced.
                for (thisReplace in 1:numel(replaceCells)){

                    # Replace with a random possible value.
                    theseT2Error(replaceCells(thisReplace)) <- randi(nLetters)-replacePositions(thisReplace)

                }

            }

            # Unpack the parameter bounds to feed into the model fitting.
            # These are set near the top of the script, but need to be
            # unpacked for each scenario.

            # Unpack mean (latency) bounds for T1.
            mu1_lb_T1 <- -muBound
            mu1_ub_T1 <- +muBound
            mu2_lb_T1 <- listLags(thisLag)-muBound
            mu2_ub_T1 <- listLags(thisLag)+muBound

            # Unpack SD (precision) bounds for T1.
            sigma1_lb_T1 <- smallNonZeroNumber
            sigma1_ub_T1 <- sigmaBound
            sigma2_lb_T1 <- smallNonZeroNumber
            sigma2_ub_T1 <- sigmaBound

            # Unpack mean (latency) bounds for T2.
            mu1_lb_T2 <- -muBound
            mu1_ub_T2 <- +muBound
            mu2_lb_T2 <- -listLags(thisLag)-muBound
            mu2_ub_T2 <- -listLags(thisLag)+muBound

            # Unpack SD (precision) bounds for T2.
            sigma1_lb_T2 <- smallNonZeroNumber
            sigma1_ub_T2 <- sigmaBound
            sigma2_lb_T2 <- smallNonZeroNumber
            sigma2_ub_T2 <- sigmaBound

            # Unpack mean (latency) bounds for compined T1 & T2.
            mu1_lb_T1T2 <- -muBound
            mu1_ub_T1T2 <- +muBound
            mu2_lb_T1T2 <- -listLags(thisLag)-muBound
            mu2_ub_T1T2 <- -listLags(thisLag)+muBound

            # Unpack SD (precision) bounds for combined T1 & T2.
            sigma1_lb_T1T2 <- smallNonZeroNumber
            sigma1_ub_T1T2 <- sigmaBound
            sigma2_lb_T1T2 <- smallNonZeroNumber
            sigma2_ub_T1T2 <- sigmaBound

            # Fit the model to the T1 distribution.

            # Keep track of the minimum negative log likelihood on each
            # replicate. Start at infinity so the first replicate
            # automatically qualifies as the best candidate up to that
            # point.
            minNegLogLikelihood <- inf

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
            pseudo_uniform in zeros(size(xDomain))){){

            # Cycle through each possible T1 position.
            for (thisPosNo in 1:numel(theseT1Pos)){

                # Identify the actual T1 position corresponding to the
                # position number. For example, the first position number
                # might be the 7th position in the stream.
                thisPos <- theseT1Pos(thisPosNo)

                # Add to the pseudo-uniform distribution one unit for every
                # possible error given that T1 position.
                pseudo_uniform((1-thisPos-minErr+1):(nLetters-thisPos-minErr+1)) in pseudo_uniform((1-thisPos-minErr+1):(nLetters-thisPos-minErr+1))+ones(1,nLetters)){){

            }

            # Cycle through a number of replicates of the fitting
            # procedure with different randomised starting values across
            # the range dictated by the bounds.

            for (thisReplicate in 1:nReplicates){

                # Randomise starting values for each parameter.
                p1Guess <- rand
                mu1Guess <- (2*muBound*rand)-muBound
                sigma1Guess <- sigmaBound*rand
                p2Guess <- rand
                mu2Guess <- listLags(thisLag)+(2*muBound*rand)-muBound
                sigma2Guess <- sigmaBound*rand

                # Ensure guesses satisfy bounds, and round them marginally
                # up or down if necessary.
                p1Guess <- max([smallNonZeroNumber p1Guess])
                mu1Guess <- min([mu1_ub_T1 max([mu1_lb_T1 mu1Guess])])
                sigma1Guess <- min([sigma1_ub_T1 max([sigma1_lb_T1 sigma1Guess])])
                p2Guess <- max([smallNonZeroNumber p2Guess])
                mu2Guess <- min([mu2_ub_T1 max([mu2_lb_T1 mu2Guess])])
                sigma2Guess <- min([sigma2_ub_T1 max([sigma2_lb_T1 sigma2Guess])])

                # Compile to feed into the MLE function.
                parameterGuess <- [p1Guess mu1Guess sigma1Guess p2Guess mu2Guess sigma2Guess]
                parameterLowerBound <- [smallNonZeroNumber mu1_lb_T1 sigma1_lb_T1 smallNonZeroNumber mu2_lb_T1 sigma2_lb_T1]
                parameterUpperBound <- [1 mu1_ub_T1 sigma1_ub_T1 1 mu2_ub_T1 sigma2_ub_T1]

                # Run the MLE function.
                [currentEstimates, currentCIs] <- mle(theseT1Error, 'pdf', pdf_normmixture, 'start', parameterGuess, 'lower', parameterLowerBound, 'upper', parameterUpperBound, 'options', options)

                # Compute the negative log likelihood of the fitted model.
                thisNegLogLikelihood <- -sum(log(pdf_normmixture(theseT1Error,currentEstimates(1),currentEstimates(2),currentEstimates(3),currentEstimates(4),currentEstimates(5),currentEstimates(6))))

                # Check whether this is lower than the lowest so far.
                if (minNegLogLikelihood > thisNegLogLikelihood){

                    # If so, store this as the current best estimate.
                    minNegLogLikelihood <- thisNegLogLikelihood
                    bestEstimates <- currentEstimates
                    bestEstimateCIs <- currentCIs

                }

            }

            # Recover individual efficacy estimates. To ensure the sum of
            # efficacies across the two episodes is no greater than 1, the
            # objective function uses p1 as the sum of the two, and p2 as
            # the relative proportions.
            p1 <- bestEstimates(1)
            p2 <- bestEstimates(4)
            bestEstimates(1) <- p1*(1-p2)
            bestEstimates(4) <- p1*p2

            # Do the same for the confidence intervals on the estimates.
            p1 <- bestEstimateCIs(:,1)
            p2 <- bestEstimateCIs(:,4)
            bestEstimateCIs(:,1) <- p1.*(1-p2)
            bestEstimateCIs(:,4) <- p1.*p2

            # Check whether a switch in the order of parameters is
            # required. When the parameter limits overlap, T1 and and
            # T2 episode estimates can get switched around. We take the
            # later one to be T2.
            if (bestEstimates(2) > bestEstimates(5)){
                bestEstimates <- bestEstimates([4 5 6 1 2 3])
                bestEstimateCIs <- bestEstimateCIs(:,[4 5 6 1 2 3])
            }

            # Enter the best estimates into the parameter matrices.
            allT1Estimates_byParticipant(thisParticipant,thisLag,:) <- bestEstimates
            allT1LowerBounds_byParticipant(thisParticipant,thisLag,:) <- bestEstimateCIs(1,:)
            allT1UpperBounds_byParticipant(thisParticipant,thisLag,:) <- bestEstimateCIs(2,:)
            allT1MinNegLogLikelihoods_byParticipant(thisParticipant,thisLag) <- minNegLogLikelihood

            # Plot the distributions if required.
            if (plotFits){
                pFig <- figure('Color','white','Name',['Participant ' num2str(thisParticipant) ', Lag ' num2str(thisLag)])##ok<UNRCH>
                subplot(1,3,1)
                tBars <- hist(theseT1Error, xDomain)
                tBars <- tBars/sum(tBars)
                bar(xDomain,tBars)
                hold on
                p1 <- bestEstimates(1) + bestEstimates(4)
                p2 <- bestEstimates(4)/p1
                yModelPlot <- pdf_normmixture(xDomain, p1, bestEstimates(2), bestEstimates(3), p2, bestEstimates(5), bestEstimates(6))
                plot(xDomain,yModelPlot,'r--')
                axis square
                axis([minErr-1 maxErr+1 0 0.75])
                title('T1')
                drawnow
            }

            # Fit the model to the T2 distribution.
            minT2Error <- min(theseT2Error)
            maxT2Error <- max(theseT2Error)

            # Keep track of the minimum negative log likelihood on each
            # replicate.
            minNegLogLikelihood <- inf

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
            pseudo_uniform in zeros(size(xDomain))){){

            # Cycle through each possible T2 position.
            for (thisPosNo in 1:numel(theseT2Pos)){

                # Identify the actual T2 position corresponding to the
                # position number. For example, the first position number
                # might be the 7th position in the stream.
                thisPos <- theseT2Pos(thisPosNo)

                # Add to the pseudo-uniform distribution one unit for every
                # possible error given that T2 position.
                pseudo_uniform((1-thisPos-minErr+1):(nLetters-thisPos-minErr+1)) in pseudo_uniform((1-thisPos-minErr+1):(nLetters-thisPos-minErr+1))+ones(1,nLetters)){){

            }

            # Cycle through a number of replicates of the fitting
            # procedure with different randomised starting values across
            # the range dictated by the bounds.

            for (thisReplicate in 1:nReplicates){

                # Randomise starting values for each parameter.
                p1Guess <- rand
                mu1Guess <- (2*muBound*rand)-muBound
                sigma1Guess <- sigmaBound*rand+smallNonZeroNumber
                p2Guess <- rand
                mu2Guess <- -listLags(thisLag)+(2*muBound*rand)-muBound
                sigma2Guess <- sigmaBound*rand+smallNonZeroNumber

                # Ensure guesses satisfy bounds, and round them marginally
                # up or down if necessary.
                p1Guess <- max([smallNonZeroNumber p1Guess])
                mu1Guess <- min([mu1_ub_T2 max([mu1_lb_T2 mu1Guess])])
                sigma1Guess <- min([sigma1_ub_T2 max([sigma1_lb_T2 sigma1Guess])])
                p2Guess <- max([smallNonZeroNumber p2Guess])
                mu2Guess <- min([mu2_ub_T2 max([mu2_lb_T2 mu2Guess])])
                sigma2Guess <- min([sigma2_ub_T2 max([sigma2_lb_T2 sigma2Guess])])

                # Compile to feed into the MLE function.
                parameterGuess <- [p1Guess mu1Guess sigma1Guess p2Guess mu2Guess sigma2Guess]
                parameterLowerBound <- [smallNonZeroNumber mu1_lb_T2 sigma1_lb_T2 smallNonZeroNumber mu2_lb_T2 sigma2_lb_T2]
                parameterUpperBound <- [1 mu1_ub_T2 sigma1_ub_T2 1 mu2_ub_T2 sigma2_ub_T2]

                # Run the MLE function.
                [currentEstimates, currentCIs] <- mle(theseT2Error, 'pdf', pdf_normmixture, 'start', parameterGuess, 'lower', parameterLowerBound, 'upper', parameterUpperBound, 'options', options)

                # Compute the negative log likelihood of the fitted model.
                thisNegLogLikelihood <- -sum(log(pdf_normmixture(theseT2Error,currentEstimates(1),currentEstimates(2),currentEstimates(3),currentEstimates(4),currentEstimates(5),currentEstimates(6))))

                # Check whether this is lower than the lowest so far.
                if (minNegLogLikelihood > thisNegLogLikelihood){

                    # If so, store this as the current best estimate.
                    minNegLogLikelihood <- thisNegLogLikelihood
                    bestEstimates <- currentEstimates
                    bestEstimateCIs <- currentCIs

                }

            }

            # Recover individual efficacy estimates.
            p1 <- bestEstimates(1)
            p2 <- bestEstimates(4)
            bestEstimates(1) <- p1*(1-p2)
            bestEstimates(4) <- p1*p2

            # Recover confidence intervals on the estimates.
            p1 <- bestEstimateCIs(:,1)
            p2 <- bestEstimateCIs(:,4)
            bestEstimateCIs(:,1) <- p1.*(1-p2)
            bestEstimateCIs(:,4) <- p1.*p2

            # Check whether a switch in the order of parameters is
            # required.
            if (bestEstimates(2) < bestEstimates(5)){
                bestEstimates <- bestEstimates([4 5 6 1 2 3])
                bestEstimateCIs <- bestEstimateCIs(:,[4 5 6 1 2 3])
                fprintf('(T2 swap) ')
            }

            # Enter the best estimates into the parameter matrices.
            allT2Estimates_byParticipant(thisParticipant,thisLag,:) <- bestEstimates
            allT2LowerBounds_byParticipant(thisParticipant,thisLag,:) <- bestEstimateCIs(1,:)
            allT2UpperBounds_byParticipant(thisParticipant,thisLag,:) <- bestEstimateCIs(2,:)
            allT2MinNegLogLikelihoods_byParticipant(thisParticipant,thisLag) <- minNegLogLikelihood

            # Plot the distributions if required.
            if (plotFits){
                subplot(1,3,2)##ok<UNRCH>
                tBars <- hist(theseT2Error, xDomain)
                tBars <- tBars/sum(tBars)
                bar(xDomain,tBars)
                hold on
                p1 <- bestEstimates(1) + bestEstimates(4)
                p2 <- bestEstimates(4)/p1
                yModelPlot <- pdf_normmixture(xDomain, p1, bestEstimates(2), bestEstimates(3), p2, bestEstimates(5), bestEstimates(6))
                plot(xDomain,yModelPlot,'r--')
                axis square
                axis([minErr-1 maxErr+1 0 0.75])
                title('T2')
                drawnow
            }

            # Fit the model to the combined T1 + T2 distribution. First,
            # combine the distributions on a common scale (i.e. alter T1
            # errors to reflect position relative to T2.
            theseT1T2Error <- [theseT1Error-listLags(thisLag) theseT2Error]
            minT1T2Error <- min(theseT1T2Error)
            maxT1T2Error <- max(theseT1T2Error)

            # Keep track of the minimum negative log likelihood on each
            # replicate.
            minNegLogLikelihood <- inf

            # Calculate the domain of possible errors (xDomain).
            theseT1T2Pos <- theseT2Pos# Errors are relative to T2
            xPosition <- unique(theseT1T2Pos)
            minX_T1T2 <- min(xPosition)
            maxX_T1T2 <- max(xPosition)
            minErr <- 1-maxX_T1T2-1
            maxErr <- nLetters-minX_T1T2+1
            xDomain <- minErr:maxErr

            # Generate the 'pseudo-uniform' distribution, which is the
            # expected distribution of errors if a random guess was
            # provided on every trial.
            pseudo_uniform in zeros(size(xDomain))){){

            # Cycle through each possible T1 and T2 position.
            for (thisPosNo in 1:numel(theseT1T2Pos)){

                # Identify the actual T1 or T2 position corresponding to the
                # position number. For example, the first position number
                # might be the 7th position in the stream.
                thisPos <- theseT1T2Pos(thisPosNo)

                # Add to the pseudo-uniform distribution one unit for every
                # possible error given that T1 or T2 position.
                pseudo_uniform((1-thisPos-minErr+1):(nLetters-thisPos-minErr+1)) in pseudo_uniform((1-thisPos-minErr+1):(nLetters-thisPos-minErr+1))+ones(1,nLetters)){){

            }

            # Cycle through a number of replicates of the fitting
            # procedure with different randomised starting values across
            # the range dictated by the bounds.

            for (thisReplicate in 1:nReplicates){

                # Randomise starting values for each parameter.
                p1Guess <- rand
                mu1Guess <- (2*muBound*rand)-muBound
                sigma1Guess <- sigmaBound*rand
                p2Guess <- rand
                mu2Guess <- -listLags(thisLag)+(2*muBound*rand)-muBound
                sigma2Guess <- sigmaBound*rand

                # Ensure guesses satisfy bounds, and round them marginally
                # up or down if necessary.
                p1Guess <- max([smallNonZeroNumber p1Guess])
                mu1Guess <- min([mu1_ub_T1T2 max([mu1_lb_T1T2 mu1Guess])])
                sigma1Guess <- min([sigma1_ub_T1T2 max([sigma1_lb_T1T2 sigma1Guess])])
                p2Guess <- max([smallNonZeroNumber p2Guess])
                mu2Guess <- min([mu2_ub_T1T2 max([mu2_lb_T1T2 mu2Guess])])
                sigma2Guess <- min([sigma2_ub_T1T2 max([sigma2_lb_T1T2 sigma2Guess])])

                # Compile to feed into the MLE function.
                parameterGuess <- [p1Guess mu1Guess sigma1Guess p2Guess mu2Guess sigma2Guess]
                parameterLowerBound <- [smallNonZeroNumber mu1_lb_T1T2 sigma1_lb_T1T2 smallNonZeroNumber mu2_lb_T1T2 sigma2_lb_T1T2]
                parameterUpperBound <- [1 mu1_ub_T1T2 sigma1_ub_T1T2 1 mu2_ub_T1T2 sigma2_ub_T1T2]

                # Run the MLE function.
                [currentEstimates, currentCIs] <- mle(theseT1T2Error, 'pdf', pdf_normmixture, 'start', parameterGuess, 'lower', parameterLowerBound, 'upper', parameterUpperBound, 'options', options)

                # Compute the negative log likelihood of the fitted model.
                thisNegLogLikelihood <- -sum(log(pdf_normmixture(theseT1T2Error,currentEstimates(1),currentEstimates(2),currentEstimates(3),currentEstimates(4),currentEstimates(5),currentEstimates(6))))

                # Check whether this is lower than the lowest so far.
                if (minNegLogLikelihood > thisNegLogLikelihood){

                    # If so, store this as the current best estimate.
                    minNegLogLikelihood <- thisNegLogLikelihood
                    bestEstimates <- currentEstimates
                    bestEstimateCIs <- currentCIs

                }

            }

            # Recover individual efficacy estimates.
            p1 <- bestEstimates(1)
            p2 <- bestEstimates(4)
            bestEstimates(1) <- p1*(1-p2)
            bestEstimates(4) <- p1*p2

            # Recover confidence intervals on the estimates.
            p1 <- bestEstimateCIs(:,1)
            p2 <- bestEstimateCIs(:,4)
            bestEstimateCIs(:,1) <- p1.*(1-p2)
            bestEstimateCIs(:,4) <- p1.*p2

            # Check whether a switch in the order of parameters is
            # required.
            if (bestEstimates(2) < bestEstimates(5)){
                bestEstimates <- bestEstimates([4 5 6 1 2 3])
                bestEstimateCIs <- bestEstimateCIs(:,[4 5 6 1 2 3])
                fprintf('(T1+T2 swap) ')
            }

            # Enter the best estimates into the parameter matrices.
            allT1T2Estimates_byParticipant(thisParticipant,thisLag,:) <- bestEstimates
            allT1T2LowerBounds_byParticipant(thisParticipant,thisLag,:) <- bestEstimateCIs(1,:)
            allT1T2UpperBounds_byParticipant(thisParticipant,thisLag,:) <- bestEstimateCIs(2,:)
            allT1T2MinNegLogLikelihoods_byParticipant(thisParticipant,thisLag) <- minNegLogLikelihood

            # Plot the distributions if required.
            if (plotFits){
                subplot(1,3,3)##ok<UNRCH>
                tBars <- hist(theseT1T2Error, xDomain)
                tBars <- tBars/sum(tBars)
                bar(xDomain,tBars)
                hold on
                p1 <- bestEstimates(1) + bestEstimates(4)
                p2 <- bestEstimates(4)/p1
                yModelPlot <- pdf_normmixture(xDomain, p1, bestEstimates(2), bestEstimates(3), p2, bestEstimates(5), bestEstimates(6))
                plot(xDomain,yModelPlot,'r--')
                axis square
                axis([minErr-1 maxErr+1 0 0.75])
                title('T1+T2')
                drawnow
            }

        }

    }

    # Keep track of progress in the command window.
    fprintf('\n\n')

    # Change directory to store model output.
    cd([thisPath 'ModelOutput'])

    # Save model output.
    save(['ModelOutput_' sampleNames{thisSample} '_Dual.mat'])

}

# Turn this warning back on.
warning('on', 'stats:mlecov:NonPosDefHessian')
