% This script takes a dataset from an attentional blink task and generates 
% plots showing the distribution of serial position errors with the models 
% fitted to that data. For illustrative purposes, we combine across 
% participants within the dataset.
%
% The script duplicates some of the functionality of AB_Model_Single,
% AB_Model_Dual, and AB_Compare_Models, but does the fitting to data
% pooled across participants.
%
% You should have a folder named 'Data', where the compiled data from 
% that sample should be. The compiled data should be called something like 
% 'CompiledData_ThisSample.mat'.
%
% The script requires the Econometrics Toolbox, specifically the function
% aicbic(), for calculation of the Bayesian Information Criterion. It also
% requires the Statistics Toolbox, for the nanmean() function, and possibly
% some other things.
%
% The objective functions pdf_Mixture_Single.m and pdf_Mixture_Dual.m
% should be on the MATLAB path.

% -------------------------------------------------------------------------
% M1 (single model fit)
% -------------------------------------------------------------------------
% This section contains some things that will need to be set locally,
% depending on your directory structure and your data.

% Provide the path.
thisPath = '/Users/experimentalmode/Documents/MATLAB/ABFinal/';

% Provide a name for the sample you want to plot, so files can be read with
% corresponding filenames.
sampleName = 'Warwick';

% Provide some properties of the data for the sample.
nParticipants = 20;          % Number of participants
nLetters = 26;               % Number of items in a stream on each trial
nTrials = 100;               % Total number of trials per block across all lags
nBlocks = 4;                 % Total number of blocks
itemRate = 11.125;           % Item rate (items/sec)

% Set some model-fitting parameters.
nReplicates = 100;                          % Number of times to repeat each fit with different starting values
smallNonZeroNumber = 10^-5;                 % Useful number for when limits can't be exactly zero but can be anything larger
fitMaxIter = 10^4;                          % Maximum number of fit iterations
fitMaxFunEvals = 10^4;                      % Maximum number of model evaluations

% Set some parameter bounds. You want these large enough that they span the
% full reasonable range of mean (latency) and SD (precision), but small
% enough to prevent over-fitting to blips in the distributions. These
% values are about right in most cases, but might need some tweaking if
% e.g. you were analysing data with an unusually high or low item rate.
muBound = 4;
sigmaBound = 4;

% Ordinarily you wouldn't want to change these, but you might want to 
% provide a different function with a different number of parameters.
nFreeParameters_Single = 3;
nFreeParameters_Dual = 6;
pdf_normmixture_single = @pdf_Mixture_Single;
pdf_normmixture = @pdf_Mixture_Dual;

% Set the colours used in the plots.
plotColor = {'r','c'};
plotAxes = [-24 24 -.075 .825];
gridPoints = 500; % Number of values to calculate when plotting the model

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% Declare global variables that need to be accessed by the objective
% function.
global xDomain;
global pseudo_uniform;

% Add folders to the MATLAB path.
addpath(genpath(thisPath));
cd([thisPath 'Data']);

% Turn off this warning to prevent it from filling the command
% window. This is only a problem if you get it on most
% replicates of a fit.
warning('off', 'stats:mlecov:NonPosDefHessian');

% Load the compiled data file for this sample.
load(['CompiledData_' sampleName '.mat']);

% Work out the number of lags.
listLags = unique(allLags(:));
listLags(isnan(listLags)) = [];
nLags = numel(listLags);

% Work out possible positions in the stream for T1.
listT1Pos = unique(allT1Pos(:));
nT1Pos = numel(listT1Pos);

% Get the list of T1 errors and extract some properties.
listT1Errors = unique(allT1Error(:));       % List of unique observed T1 errors
listT1Errors(isnan(listT1Errors)) = [];     % Get rid of NaN values
nT1Errors = numel(listT1Errors);            % Number of unique observed T1 errors
minT1Error = min(listT1Errors);             % Lowest (most negative) T1 error
maxT1Error = max(listT1Errors);             % Highest (most positive) T1 error

nTrialsPerLag = (nTrials*nBlocks)/nLags;    % Calculate number of trials per lag

% The following output to the command window is just to keep track of
% what the program is doing.
fprintf('\n\n%s: M1\n\n', upper(sampleName));

% Build empty matrices for storing parameter estimates for each
% participant, at each lag. Also build matrices to store upper and
% lower bounds, and minimum negative log likelihoods.

allT1Estimates_M1 = NaN(nLags,nFreeParameters_Single);
allT1LowerBounds_M1 = NaN(nLags,nFreeParameters_Single);
allT1UpperBounds_M1 = NaN(nLags,nFreeParameters_Single);
allT1MinNegLogLikelihoods_M1 = NaN(nLags,1);

allT2Estimates_M1 = NaN(nLags,nFreeParameters_Single);
allT2LowerBounds_M1 = NaN(nLags,nFreeParameters_Single);
allT2UpperBounds_M1 = NaN(nLags,nFreeParameters_Single);
allT2MinNegLogLikelihoods_M1 = NaN(nLags,1);

allT1T2Estimates_M1 = NaN(nLags,nFreeParameters_Single);
allT1T2LowerBounds_M1 = NaN(nLags,nFreeParameters_Single);
allT1T2UpperBounds_M1 = NaN(nLags,nFreeParameters_Single);
allT1T2MinNegLogLikelihoods_M1 = NaN(nLags,1);

allT1Estimates_M2 = NaN(nLags,nFreeParameters_Dual);
allT1LowerBounds_M2 = NaN(nLags,nFreeParameters_Dual);
allT1UpperBounds_M2 = NaN(nLags,nFreeParameters_Dual);
allT1MinNegLogLikelihoods_M2 = NaN(nLags,1);

allT2Estimates_M2 = NaN(nLags,nFreeParameters_Dual);
allT2LowerBounds_M2 = NaN(nLags,nFreeParameters_Dual);
allT2UpperBounds_M2 = NaN(nLags,nFreeParameters_Dual);
allT2MinNegLogLikelihoods_M2 = NaN(nLags,1);

allT1T2Estimates_M2 = NaN(nLags,nFreeParameters_Dual);
allT1T2LowerBounds_M2 = NaN(nLags,nFreeParameters_Dual);
allT1T2UpperBounds_M2 = NaN(nLags,nFreeParameters_Dual);
allT1T2MinNegLogLikelihoods_M2 = NaN(nLags,1);

% Set fit options.
options = statset('MaxIter', fitMaxIter, 'MaxFunEvals', fitMaxFunEvals, 'Display', 'off');

% Extract the relevant lists of T1 and T2 errors, T1 and T2 stream
% positions, and corresponding lags.
T1Error = allT1Error(:);
T2Error = allT2Error(:);
T1Pos = allT1Pos(:);
T2Pos = allT2Pos(:);
Lags = allLags(:);

% Cycle through each lag.
for thisLag = 1:nLags

    % Keep track of progress in the command window.
    fprintf('L%d ', thisLag);

    % Find the lag value for this numbered lag. These are usually
    % the same, but this allows for the possibility that some lag
    % values aren't tested.
    thisLagVal = listLags(thisLag);

    % Go through and replace NaN values. Here, we assume that a NaN
    % value is a guess, and replace it with a random sample from
    % the range of possible errors. You might want to think about
    % this if you have data with lots of missing values, but it
    % makes virtually no difference in these data.

    % Identify T1 and T2 trials at this lag.
    hasThisLag = Lags==thisLagVal;
    theseT1Error = T1Error(hasThisLag)';
    theseT2Error = T2Error(hasThisLag)';
    theseT1Pos = T1Pos(hasThisLag)';
    theseT2Pos = T2Pos(hasThisLag)';

    if sum(isnan(theseT1Error)) > 0 % If there is at least one NaN

        % Find NaNs.
        replaceCells = find(isnan(theseT1Error));

        % Find the T1 positions of the trials to be replaced.
        replacePositions = theseT1Pos(replaceCells);

        % Cycle through each point to be replaced.
        for thisReplace = 1:numel(replaceCells);

            % Replace with a random possible value.
            theseT1Error(replaceCells(thisReplace)) = randi(nLetters)-replacePositions(thisReplace);

        end

    end

    if sum(isnan(theseT2Error)) > 0 % If there is at least one NaN

        % Find NaNs.
        replaceCells = find(isnan(theseT2Error));

        % Find the T2 positions of the trials to be replaced.
        replacePositions = theseT2Pos(replaceCells);

        % Cycle through each point to be replaced.
        for thisReplace = 1:numel(replaceCells);

            % Replace with a random possible value.
            theseT2Error(replaceCells(thisReplace)) = randi(nLetters)-replacePositions(thisReplace);

        end

    end

    % Combine T1 and T2 distributions on a common scale (i.e. alter
    % T1 errors to reflect position relative to T2.
    theseT1T2Error = [theseT1Error-listLags(thisLag) theseT2Error];

    % Get minimum and maximum error values.
    minT2Error = min(theseT2Error);
    maxT2Error = max(theseT2Error);
    minT1T2Error = min(theseT1T2Error);
    maxT1T2Error = max(theseT1T2Error);

    % Unpack the parameter bounds to feed into the model fitting.
    % These are set near the top of the script, but need to be
    % unpacked for each scenario.

    % Unpack mean (latency) bounds for T1.
    mu_lb_T1 = -muBound;
    mu_ub_T1 = muBound;

    % Unpack SD (precision) bounds for T1.
    sigma_lb_T1 = smallNonZeroNumber;
    sigma_ub_T1 = sigmaBound;

    % Unpack mean (latency) bounds for T2.
    mu_lb_T2 = -muBound;
    mu_ub_T2 = sigmaBound;

    % Unpack SD (precision) bounds for T2.
    sigma_lb_T2 = smallNonZeroNumber;
    sigma_ub_T2 = sigmaBound;

    % Unpack mean (latency) bounds for compined T1 & T2.
    mu_lb_T1T2 = -muBound-listLags(thisLag);
    mu_ub_T1T2 = muBound;

    % Unpack SD (precision) bounds for combined T1 & T2.
    sigma_lb_T1T2 = smallNonZeroNumber;
    sigma_ub_T1T2 = sigmaBound;

    % Fit the model to the T1 distribution.

    % Keep track of the minimum negative log likelihood on each
    % replicate. Start at infinity so the first replicate
    % automatically qualifies as the best candidate up to that
    % point.
    minNegLogLikelihood = inf;

    % Calculate the domain of possible errors (xDomain).
    xPosition = unique(theseT1Pos);
    minX_T1 = min(xPosition);
    maxX_T1 = max(xPosition);
    minErr = 1-maxX_T1-1;
    maxErr = nLetters-minX_T1+1;
    xDomain = minErr:maxErr;

    % Generate the 'pseudo-uniform' distribution, which is the
    % expected distribution of errors if a random guess was
    % provided on every trial. This isn't an actual uniform
    % distribution because the most extreme errors are only
    % possible on trials in which targets appear at their most
    % extreme positions.
    pseudo_uniform = zeros(size(xDomain));

    % Cycle through each possible T1 position.
    for thisPosNo = 1:numel(theseT1Pos)

        % Identify the actual T1 position corresponding to the
        % position number. For example, the first position number
        % might be the 7th position in the stream.
        thisPos = theseT1Pos(thisPosNo);

        % Add to the pseudo-uniform distribution one unit for every
        % possible error given that T1 position.
        pseudo_uniform((1-thisPos-minErr+1):(nLetters-thisPos-minErr+1)) = pseudo_uniform((1-thisPos-minErr+1):(nLetters-thisPos-minErr+1))+ones(1,nLetters);

    end

    % Cycle through a number of replicates of the fitting
    % procedure with different randomised starting values across
    % the range dictated by the bounds.

    for thisReplicate = 1:nReplicates

        % Randomise starting values for each parameter.
        pGuess = max([smallNonZeroNumber rand]);
        muGuess = (2*muBound*rand)-muBound;
        sigmaGuess = sigmaBound*rand+smallNonZeroNumber;

        % Compile to feed into the MLE function.
        parameterGuess = [pGuess muGuess sigmaGuess];
        parameterLowerBound = [smallNonZeroNumber mu_lb_T1 sigma_lb_T1];
        parameterUpperBound = [1 mu_ub_T1 sigma_ub_T1];

        % Ensure guesses satisfy bounds, and round them marginally
        % up or down if necessary.
        parameterGuess = max([parameterGuess;parameterLowerBound]);
        parameterGuess = min([parameterGuess;parameterUpperBound]);

        % Run the MLE function.
        [currentEstimates, currentCIs] = mle(theseT1Error, 'pdf', pdf_normmixture_single, 'start', parameterGuess, 'lower', parameterLowerBound, 'upper', parameterUpperBound, 'options', options);

        % Compute the negative log likelihood of the fitted model.
        thisNegLogLikelihood = -sum(log(pdf_normmixture_single(theseT1Error,currentEstimates(1),currentEstimates(2),currentEstimates(3))));

        % Check whether this is lower than the lowest so far.
        if minNegLogLikelihood > thisNegLogLikelihood

            % If so, store this as the current best estimate.
            minNegLogLikelihood = thisNegLogLikelihood;
            bestEstimates = currentEstimates;
            bestEstimateCIs = currentCIs;

        end

    end

    % Enter the best estimates into the parameter matrices.
    allT1Estimates_M1(thisLag,:) = bestEstimates;
    allT1LowerBounds_M1(thisLag,:) = bestEstimateCIs(1,:);
    allT1UpperBounds_M1(thisLag,:) = bestEstimateCIs(2,:);
    allT1MinNegLogLikelihoods_M1(thisLag) = minNegLogLikelihood;

    % Fit the model to the T2 distribution.

    % Keep track of the minimum negative log likelihood on each
    % replicate.
    minNegLogLikelihood = inf;

    % Calculate the domain of possible errors (xDomain).
    xPosition = unique(theseT2Pos);
    minX_T2 = min(xPosition);
    maxX_T2 = max(xPosition);
    minErr = 1-maxX_T2-1;
    maxErr = nLetters-minX_T2+1;
    xDomain = minErr:maxErr;

    % Generate the 'pseudo-uniform' distribution, which is the
    % expected distribution of errors if a random guess was
    % provided on every trial.
    pseudo_uniform = zeros(size(xDomain));

    % Cycle through each possible T2 position.
    for thisPosNo = 1:numel(theseT2Pos)

        % Identify the actual T2 position corresponding to the
        % position number. For example, the first position number
        % might be the 7th position in the stream.
        thisPos = theseT2Pos(thisPosNo);

        % Add to the pseudo-uniform distribution one unit for every
        % possible error given that T2 position.
        pseudo_uniform((1-thisPos-minErr+1):(nLetters-thisPos-minErr+1)) = pseudo_uniform((1-thisPos-minErr+1):(nLetters-thisPos-minErr+1))+ones(1,nLetters);

    end

    % Cycle through a number of replicates of the fitting
    % procedure with different randomised starting values across
    % the range dictated by the bounds.

    for thisReplicate = 1:nReplicates

        % Randomise starting values for each parameter.
        pGuess = max([smallNonZeroNumber rand]);
        muGuess = (2*muBound*rand)-muBound;
        sigmaGuess = sigmaBound*rand+smallNonZeroNumber;

        % Compile to feed into the MLE function.
        parameterGuess = [pGuess muGuess sigmaGuess];
        parameterLowerBound = [smallNonZeroNumber mu_lb_T2 sigma_lb_T2];
        parameterUpperBound = [1 mu_ub_T2 sigma_ub_T2];

        % Ensure guesses satisfy bounds.
        parameterGuess = max([parameterGuess;parameterLowerBound]);
        parameterGuess = min([parameterGuess;parameterUpperBound]);

        % Run the MLE function.
        [currentEstimates, currentCIs] = mle(theseT2Error, 'pdf', pdf_normmixture_single, 'start', parameterGuess, 'lower', parameterLowerBound, 'upper', parameterUpperBound, 'options', options);

        % Compute the negative log likelihood of the fitted model.
        thisNegLogLikelihood = -sum(log(pdf_normmixture_single(theseT2Error,currentEstimates(1),currentEstimates(2),currentEstimates(3))));

        % Check whether this is lower than the lowest so far.
        if minNegLogLikelihood > thisNegLogLikelihood

            % If so, store this as the current best estimate.
            minNegLogLikelihood = thisNegLogLikelihood;
            bestEstimates = currentEstimates;
            bestEstimateCIs = currentCIs;

        end

    end

    % Enter the best estimates into the parameter matrices.
    allT2Estimates_M1(thisLag,:) = bestEstimates;
    allT2LowerBounds_M1(thisLag,:) = bestEstimateCIs(1,:);
    allT2UpperBounds_M1(thisLag,:) = bestEstimateCIs(2,:);
    allT2MinNegLogLikelihoods_M1(thisLag) = minNegLogLikelihood;

    % Fit the model to the combined T1 + T2 distribution. First,
    % combine the distributions on a common scale (i.e. alter T1
    % errors to reflect position relative to T2.

    % Keep track of the minimum negative log likelihood on each
    % replicate.
    minNegLogLikelihood = inf;

    % Calculate the domain of possible errors (xDomain).
    theseT1T2Pos = [theseT1Pos theseT2Pos];
    xPosition = unique(theseT1T2Pos);
    minX = min([minX_T1-listLags(thisLag) minX_T2]);
    maxX = max([maxX_T1-listLags(thisLag) maxX_T2]);
    minErr = 1-maxX-1;
    maxErr = nLetters-minX+1;
    xDomain = minErr:maxErr;

    % Generate the 'pseudo-uniform' distribution, which is the
    % expected distribution of errors if a random guess was
    % provided on every trial.
    pseudo_uniform = zeros(size(xDomain));

    % Cycle through each possible T2 position.
    for thisPosNo = 1:numel(theseT1T2Pos)

        % Identify the actual T2 position corresponding to the
        % position number. For example, the first position number
        % might be the 7th position in the stream.
        thisPos = theseT1T2Pos(thisPosNo);

        % Add to the pseudo-uniform distribution one unit for every
        % possible error given that T2 position.
        pseudo_uniform((1-thisPos-minErr+1):(nLetters-thisPos-minErr+1)) = pseudo_uniform((1-thisPos-minErr+1):(nLetters-thisPos-minErr+1))+ones(1,nLetters);

    end

    % Cycle through a number of replicates of the fitting
    % procedure with different randomised starting values across
    % the range dictated by the bounds.

    for thisReplicate = 1:nReplicates

        % Randomise starting values for each parameter.
        pGuess = max([smallNonZeroNumber rand]);
        muGuess = (2*muBound*rand)-muBound;
        sigmaGuess = sigmaBound*rand+smallNonZeroNumber;

        % Compile to feed into the MLE function.
        parameterGuess = [pGuess muGuess sigmaGuess];
        parameterLowerBound = [smallNonZeroNumber mu_lb_T1T2 sigma_lb_T1T2];
        parameterUpperBound = [1 mu_ub_T1T2 sigma_ub_T1T2];

        % Ensure guesses satisfy bounds.
        parameterGuess = max([parameterGuess;parameterLowerBound]);
        parameterGuess = min([parameterGuess;parameterUpperBound]);

        % Run the MLE function.
        [currentEstimates, currentCIs] = mle(theseT1T2Error, 'pdf', pdf_normmixture_single, 'start', parameterGuess, 'lower', parameterLowerBound, 'upper', parameterUpperBound, 'options', options);

        % Compute the negative log likelihood of the fitted model.
        thisNegLogLikelihood = -sum(log(pdf_normmixture_single(theseT1T2Error,currentEstimates(1),currentEstimates(2),currentEstimates(3))));

        % Check whether this is lower than the lowest so far.
        if minNegLogLikelihood > thisNegLogLikelihood

            % If so, store this as the current best estimate.
            minNegLogLikelihood = thisNegLogLikelihood;
            bestEstimates = currentEstimates;
            bestEstimateCIs = currentCIs;

        end

    end

    % Enter the best estimates into the parameter matrices.
    allT1T2Estimates_M1(thisLag,:) = bestEstimates;
    allT1T2LowerBounds_M1(thisLag,:) = bestEstimateCIs(1,:);
    allT1T2UpperBounds_M1(thisLag,:) = bestEstimateCIs(2,:);
    allT1T2MinNegLogLikelihoods_M1(thisLag) = minNegLogLikelihood;

end

% The following output to the command window is just to keep track of
% what the program is doing.
fprintf('\n\n%s: M2\n\n', upper(sampleName));

% Cycle through each lag.
for thisLag = 1:nLags

    % Keep track of progress in the command window.
    fprintf('L%d ', thisLag);

    % Find the lag value for this numbered lag. These are usually
    % the same, but this allows for the possibility that some lag
    % values aren't tested.
    thisLagVal = listLags(thisLag);

    % Go through and replace NaN values. Here, we assume that a NaN
    % value is a guess, and replace it with a random sample from
    % the range of possible errors. You might want to think about
    % this if you have data with lots of missing values, but it
    % makes virtually no difference in these data.

    % Identify T1 and T2 trials at this lag.
    hasThisLag = Lags==thisLagVal;
    theseT1Error = T1Error(hasThisLag)';
    theseT2Error = T2Error(hasThisLag)';
    theseT1Pos = T1Pos(hasThisLag)';
    theseT2Pos = T2Pos(hasThisLag)';

    if sum(isnan(theseT1Error)) > 0 % If there is at least one NaN

        % Find NaNs.
        replaceCells = find(isnan(theseT1Error));

        % Find the T1 positions of the trials to be replaced.
        replacePositions = theseT1Pos(replaceCells);

        % Cycle through each point to be replaced.
        for thisReplace = 1:numel(replaceCells);

            % Replace with a random possible value.
            theseT1Error(replaceCells(thisReplace)) = randi(nLetters)-replacePositions(thisReplace);

        end

    end

    if sum(isnan(theseT2Error)) > 0 % If there is at least one NaN

        % Find NaNs.
        replaceCells = find(isnan(theseT2Error));

        % Find the T2 positions of the trials to be replaced.
        replacePositions = theseT2Pos(replaceCells);

        % Cycle through each point to be replaced.
        for thisReplace = 1:numel(replaceCells);

            % Replace with a random possible value.
            theseT2Error(replaceCells(thisReplace)) = randi(nLetters)-replacePositions(thisReplace);

        end

    end

    % Unpack the parameter bounds to feed into the model fitting.
    % These are set near the top of the script, but need to be
    % unpacked for each scenario.

    % Unpack mean (latency) bounds for T1.
    mu1_lb_T1 = -muBound;
    mu1_ub_T1 = +muBound;
    mu2_lb_T1 = listLags(thisLag)-muBound;
    mu2_ub_T1 = listLags(thisLag)+muBound;

    % Unpack SD (precision) bounds for T1.
    sigma1_lb_T1 = smallNonZeroNumber;
    sigma1_ub_T1 = sigmaBound;
    sigma2_lb_T1 = smallNonZeroNumber;
    sigma2_ub_T1 = sigmaBound;

    % Unpack mean (latency) bounds for T2.
    mu1_lb_T2 = -muBound;
    mu1_ub_T2 = +muBound;
    mu2_lb_T2 = -listLags(thisLag)-muBound;
    mu2_ub_T2 = -listLags(thisLag)+muBound;

    % Unpack SD (precision) bounds for T2.
    sigma1_lb_T2 = smallNonZeroNumber;
    sigma1_ub_T2 = sigmaBound;
    sigma2_lb_T2 = smallNonZeroNumber;
    sigma2_ub_T2 = sigmaBound;

    % Unpack mean (latency) bounds for compined T1 & T2.
    mu1_lb_T1T2 = -muBound;
    mu1_ub_T1T2 = +muBound;
    mu2_lb_T1T2 = -listLags(thisLag)-muBound;
    mu2_ub_T1T2 = -listLags(thisLag)+muBound;

    % Unpack SD (precision) bounds for combined T1 & T2.
    sigma1_lb_T1T2 = smallNonZeroNumber;
    sigma1_ub_T1T2 = sigmaBound;
    sigma2_lb_T1T2 = smallNonZeroNumber;
    sigma2_ub_T1T2 = sigmaBound;

    % Fit the model to the T1 distribution.

    % Keep track of the minimum negative log likelihood on each
    % replicate. Start at infinity so the first replicate
    % automatically qualifies as the best candidate up to that
    % point.
    minNegLogLikelihood = inf;

    % Calculate the domain of possible errors (xDomain).
    xPosition = unique(T1Pos);
    minX_T1 = min(xPosition);
    maxX_T1 = max(xPosition);
    minErr = 1-maxX_T1-1;
    maxErr = nLetters-minX_T1+1;
    xDomain = minErr:maxErr;

    % Generate the 'pseudo-uniform' distribution, which is the
    % expected distribution of errors if a random guess was
    % provided on every trial. This isn't an actual uniform
    % distribution because the most extreme errors are only
    % possible on trials in which targets appear at their most
    % extreme positions.
    pseudo_uniform = zeros(size(xDomain));

    % Cycle through each possible T1 position.
    for thisPosNo = 1:numel(T1Pos)

        % Identify the actual T1 position corresponding to the
        % position number. For example, the first position number
        % might be the 7th position in the stream.
        thisPos = T1Pos(thisPosNo);

        % Add to the pseudo-uniform distribution one unit for every
        % possible error given that T1 position.
        pseudo_uniform((1-thisPos-minErr+1):(nLetters-thisPos-minErr+1)) = pseudo_uniform((1-thisPos-minErr+1):(nLetters-thisPos-minErr+1))+ones(1,nLetters);

    end

    % Cycle through a number of replicates of the fitting
    % procedure with different randomised starting values across
    % the range dictated by the bounds.

    for thisReplicate = 1:nReplicates

        % Randomise starting values for each parameter.
        p1Guess = rand;                                 
        mu1Guess = (2*muBound*rand)-muBound;                         
        sigma1Guess = sigmaBound*rand;          
        p2Guess = rand;                                     
        mu2Guess = listLags(thisLag)+(2*muBound*rand)-muBound;  
        sigma2Guess = sigmaBound*rand;

        % Ensure guesses satisfy bounds, and round them marginally
        % up or down if necessary.
        p1Guess = max([smallNonZeroNumber p1Guess]);
        mu1Guess = min([mu1_ub_T1 max([mu1_lb_T1 mu1Guess])]);
        sigma1Guess = min([sigma1_ub_T1 max([sigma1_lb_T1 sigma1Guess])]);
        p2Guess = max([smallNonZeroNumber p2Guess]);
        mu2Guess = min([mu2_ub_T1 max([mu2_lb_T1 mu2Guess])]);
        sigma2Guess = min([sigma2_ub_T1 max([sigma2_lb_T1 sigma2Guess])]);

        % Compile to feed into the MLE function.
        parameterGuess = [p1Guess mu1Guess sigma1Guess p2Guess mu2Guess sigma2Guess];
        parameterLowerBound = [smallNonZeroNumber mu1_lb_T1 sigma1_lb_T1 smallNonZeroNumber mu2_lb_T1 sigma2_lb_T1];
        parameterUpperBound = [1 mu1_ub_T1 sigma1_ub_T1 1 mu2_ub_T1 sigma2_ub_T1];

        % Run the MLE function.
        [currentEstimates, currentCIs] = mle(theseT1Error, 'pdf', pdf_normmixture, 'start', parameterGuess, 'lower', parameterLowerBound, 'upper', parameterUpperBound, 'options', options);

        % Compute the negative log likelihood of the fitted model.
        thisNegLogLikelihood = -sum(log(pdf_normmixture(theseT1Error,currentEstimates(1),currentEstimates(2),currentEstimates(3),currentEstimates(4),currentEstimates(5),currentEstimates(6))));

        % Check whether this is lower than the lowest so far.
        if minNegLogLikelihood > thisNegLogLikelihood

            % If so, store this as the current best estimate.
            minNegLogLikelihood = thisNegLogLikelihood;
            bestEstimates = currentEstimates;
            bestEstimateCIs = currentCIs;

        end

    end

    % Recover individual efficacy estimates. To ensure the sum of
    % efficacies across the two episodes is no greater than 1, the
    % objective function uses p1 as the sum of the two, and p2 as
    % the relative proportions.
    p1 = bestEstimates(1);
    p2 = bestEstimates(4);
    bestEstimates(1) = p1*(1-p2);
    bestEstimates(4) = p1*p2;

    % Do the same for the confidence intervals on the estimates.
    p1 = bestEstimateCIs(:,1);
    p2 = bestEstimateCIs(:,4);
    bestEstimateCIs(:,1) = p1.*(1-p2);
    bestEstimateCIs(:,4) = p1.*p2;

    % Check whether a switch in the order of parameters is
    % required. When the parameter limits overlap, T1 and and
    % T2 episode estimates can get switched around. We take the
    % later one to be T2.
    if bestEstimates(2) > bestEstimates(5)
        bestEstimates = bestEstimates([4 5 6 1 2 3]);
        bestEstimateCIs = bestEstimateCIs(:,[4 5 6 1 2 3]);
    end

    % Enter the best estimates into the parameter matrices.
    allT1Estimates_M2(thisLag,:) = bestEstimates;
    allT1LowerBounds_M2(thisLag,:) = bestEstimateCIs(1,:);
    allT1UpperBounds_M2(thisLag,:) = bestEstimateCIs(2,:);
    allT1MinNegLogLikelihoods_M2(thisLag) = minNegLogLikelihood;

    % Fit the model to the T2 distribution.
    minT2Error = min(theseT2Error);
    maxT2Error = max(theseT2Error);

    % Keep track of the minimum negative log likelihood on each
    % replicate.
    minNegLogLikelihood = inf;

    % Calculate the domain of possible errors (xDomain).
    xPosition = unique(theseT2Pos);
    minX_T2 = min(xPosition);
    maxX_T2 = max(xPosition);
    minErr = 1-maxX_T2-1;
    maxErr = nLetters-minX_T2+1;
    xDomain = minErr:maxErr;

    % Generate the 'pseudo-uniform' distribution, which is the
    % expected distribution of errors if a random guess was
    % provided on every trial.
    pseudo_uniform = zeros(size(xDomain));

    % Cycle through each possible T2 position.
    for thisPosNo = 1:numel(theseT2Pos)

        % Identify the actual T2 position corresponding to the
        % position number. For example, the first position number
        % might be the 7th position in the stream.
        thisPos = theseT2Pos(thisPosNo);

        % Add to the pseudo-uniform distribution one unit for every
        % possible error given that T2 position.
        pseudo_uniform((1-thisPos-minErr+1):(nLetters-thisPos-minErr+1)) = pseudo_uniform((1-thisPos-minErr+1):(nLetters-thisPos-minErr+1))+ones(1,nLetters);

    end

    % Cycle through a number of replicates of the fitting
    % procedure with different randomised starting values across
    % the range dictated by the bounds.

    for thisReplicate = 1:nReplicates

        % Randomise starting values for each parameter.
        p1Guess = rand;
        mu1Guess = (2*muBound*rand)-muBound;              
        sigma1Guess = sigmaBound*rand+smallNonZeroNumber;             
        p2Guess = rand;                                           
        mu2Guess = -listLags(thisLag)+(2*muBound*rand)-muBound;  
        sigma2Guess = sigmaBound*rand+smallNonZeroNumber;

        % Ensure guesses satisfy bounds, and round them marginally
        % up or down if necessary.
        p1Guess = max([smallNonZeroNumber p1Guess]);
        mu1Guess = min([mu1_ub_T2 max([mu1_lb_T2 mu1Guess])]);
        sigma1Guess = min([sigma1_ub_T2 max([sigma1_lb_T2 sigma1Guess])]);
        p2Guess = max([smallNonZeroNumber p2Guess]);
        mu2Guess = min([mu2_ub_T2 max([mu2_lb_T2 mu2Guess])]);
        sigma2Guess = min([sigma2_ub_T2 max([sigma2_lb_T2 sigma2Guess])]);

        % Compile to feed into the MLE function.
        parameterGuess = [p1Guess mu1Guess sigma1Guess p2Guess mu2Guess sigma2Guess];
        parameterLowerBound = [smallNonZeroNumber mu1_lb_T2 sigma1_lb_T2 smallNonZeroNumber mu2_lb_T2 sigma2_lb_T2];
        parameterUpperBound = [1 mu1_ub_T2 sigma1_ub_T2 1 mu2_ub_T2 sigma2_ub_T2];

        % Run the MLE function.
        [currentEstimates, currentCIs] = mle(theseT2Error, 'pdf', pdf_normmixture, 'start', parameterGuess, 'lower', parameterLowerBound, 'upper', parameterUpperBound, 'options', options);

        % Compute the negative log likelihood of the fitted model.
        thisNegLogLikelihood = -sum(log(pdf_normmixture(theseT2Error,currentEstimates(1),currentEstimates(2),currentEstimates(3),currentEstimates(4),currentEstimates(5),currentEstimates(6))));

        % Check whether this is lower than the lowest so far.
        if minNegLogLikelihood > thisNegLogLikelihood

            % If so, store this as the current best estimate.
            minNegLogLikelihood = thisNegLogLikelihood;
            bestEstimates = currentEstimates;
            bestEstimateCIs = currentCIs;

        end

    end

    % Recover individual efficacy estimates.
    p1 = bestEstimates(1);
    p2 = bestEstimates(4);
    bestEstimates(1) = p1*(1-p2);
    bestEstimates(4) = p1*p2;

    % Recover confidence intervals on the estimates.
    p1 = bestEstimateCIs(:,1);
    p2 = bestEstimateCIs(:,4);
    bestEstimateCIs(:,1) = p1.*(1-p2);
    bestEstimateCIs(:,4) = p1.*p2;

    % Check whether a switch in the order of parameters is
    % required.
    if bestEstimates(2) < bestEstimates(5)
        bestEstimates = bestEstimates([4 5 6 1 2 3]);
        bestEstimateCIs = bestEstimateCIs(:,[4 5 6 1 2 3]);
        fprintf('(T2 swap) ');
    end

    % Enter the best estimates into the parameter matrices.
    allT2Estimates_M2(thisLag,:) = bestEstimates;
    allT2LowerBounds_M2(thisLag,:) = bestEstimateCIs(1,:);
    allT2UpperBounds_M2(thisLag,:) = bestEstimateCIs(2,:);
    allT2MinNegLogLikelihoods_M2(thisLag) = minNegLogLikelihood;

    % Fit the model to the combined T1 + T2 distribution. First,
    % combine the distributions on a common scale (i.e. alter T1
    % errors to reflect position relative to T2.
    theseT1T2Error = [theseT1Error-listLags(thisLag) theseT2Error];         
    minT1T2Error = min(theseT1T2Error);
    maxT1T2Error = max(theseT1T2Error);

    % Keep track of the minimum negative log likelihood on each
    % replicate.
    minNegLogLikelihood = inf;

    % Calculate the domain of possible errors (xDomain).
    theseT1T2Pos = [theseT1Pos theseT2Pos];
    xPosition = unique(theseT1T2Pos);
    minX_T1T2 = min([minX_T1-listLags(thisLag) minX_T2]);
    maxX_T1T2 = max([maxX_T1-listLags(thisLag) maxX_T2]);
    minErr = 1-maxX_T1T2-1;
    maxErr = nLetters-minX_T1T2+1;
    xDomain = minErr:maxErr;

    % Generate the 'pseudo-uniform' distribution, which is the
    % expected distribution of errors if a random guess was
    % provided on every trial.
    pseudo_uniform = zeros(size(xDomain));

    % Cycle through each possible T2 position.
    for thisPosNo = 1:numel(theseT2Pos)

        % Identify the actual T2 position corresponding to the
        % position number. For example, the first position number
        % might be the 7th position in the stream.
        thisPos = theseT2Pos(thisPosNo);

        % Add to the pseudo-uniform distribution one unit for every
        % possible error given that T2 position.
        pseudo_uniform((1-thisPos-minErr+1):(nLetters-thisPos-minErr+1)) = pseudo_uniform((1-thisPos-minErr+1):(nLetters-thisPos-minErr+1))+ones(1,nLetters);

    end

    % Cycle through a number of replicates of the fitting
    % procedure with different randomised starting values across
    % the range dictated by the bounds.

    for thisReplicate = 1:nReplicates

        % Randomise starting values for each parameter.
        p1Guess = rand;
        mu1Guess = (2*muBound*rand)-muBound;              
        sigma1Guess = sigmaBound*rand;             
        p2Guess = rand;                                         
        mu2Guess = -listLags(thisLag)+(2*muBound*rand)-muBound;  
        sigma2Guess = sigmaBound*rand;

        % Ensure guesses satisfy bounds, and round them marginally
        % up or down if necessary.
        p1Guess = max([smallNonZeroNumber p1Guess]);
        mu1Guess = min([mu1_ub_T1T2 max([mu1_lb_T1T2 mu1Guess])]);
        sigma1Guess = min([sigma1_ub_T1T2 max([sigma1_lb_T1T2 sigma1Guess])]);
        p2Guess = max([smallNonZeroNumber p2Guess]);
        mu2Guess = min([mu2_ub_T1T2 max([mu2_lb_T1T2 mu2Guess])]);
        sigma2Guess = min([sigma2_ub_T1T2 max([sigma2_lb_T1T2 sigma2Guess])]);

        % Compile to feed into the MLE function.
        parameterGuess = [p1Guess mu1Guess sigma1Guess p2Guess mu2Guess sigma2Guess];
        parameterLowerBound = [smallNonZeroNumber mu1_lb_T1T2 sigma1_lb_T1T2 smallNonZeroNumber mu2_lb_T1T2 sigma2_lb_T1T2];
        parameterUpperBound = [1 mu1_ub_T1T2 sigma1_ub_T1T2 1 mu2_ub_T1T2 sigma2_ub_T1T2];

        % Run the MLE function.
        [currentEstimates, currentCIs] = mle(theseT1T2Error, 'pdf', pdf_normmixture, 'start', parameterGuess, 'lower', parameterLowerBound, 'upper', parameterUpperBound, 'options', options);

        % Compute the negative log likelihood of the fitted model.
        thisNegLogLikelihood = -sum(log(pdf_normmixture(theseT1T2Error,currentEstimates(1),currentEstimates(2),currentEstimates(3),currentEstimates(4),currentEstimates(5),currentEstimates(6))));

        % Check whether this is lower than the lowest so far.
        if minNegLogLikelihood > thisNegLogLikelihood

            % If so, store this as the current best estimate.
            minNegLogLikelihood = thisNegLogLikelihood;
            bestEstimates = currentEstimates;
            bestEstimateCIs = currentCIs;

        end

    end

    % Recover individual efficacy estimates.
    p1 = bestEstimates(1);
    p2 = bestEstimates(4);
    bestEstimates(1) = p1*(1-p2);
    bestEstimates(4) = p1*p2;

    % Recover confidence intervals on the estimates.
    p1 = bestEstimateCIs(:,1);
    p2 = bestEstimateCIs(:,4);
    bestEstimateCIs(:,1) = p1.*(1-p2);
    bestEstimateCIs(:,4) = p1.*p2;

    % Check whether a switch in the order of parameters is
    % required.
    if bestEstimates(2) < bestEstimates(5)
        bestEstimates = bestEstimates([4 5 6 1 2 3]);
        bestEstimateCIs = bestEstimateCIs(:,[4 5 6 1 2 3]);
        fprintf('(T1+T2 swap) ');
    end

    % Enter the best estimates into the parameter matrices.
    allT1T2Estimates_M2(thisLag,:) = bestEstimates;
    allT1T2LowerBounds_M2(thisLag,:) = bestEstimateCIs(1,:);
    allT1T2UpperBounds_M2(thisLag,:) = bestEstimateCIs(2,:);
    allT1T2MinNegLogLikelihoods_M2(thisLag) = minNegLogLikelihood;

end

% Turn this warning back on.
warning('on', 'stats:mlecov:NonPosDefHessian');

% Now, for each lag, we plot the distributions for T1, T2 and T1+T2, along
% with the best-fitting model.

for thisLag = 1:nLags
    
    thisLagVal = listLags(thisLag);
    hasThisLag = Lags==thisLagVal;
    theseT1Error = T1Error(hasThisLag)';
    theseT2Error = T2Error(hasThisLag)';
    theseT1Pos = T1Pos(hasThisLag)';
    theseT2Pos = T2Pos(hasThisLag)';
    
    % Plot T1
    figure('Color','white','Name',[sampleName ': Lag ' num2str(listLags(thisLag)) ' (T1)']); % Create a figure

    % Calculate the domain of possible errors (xDomain).
    xPosition = unique(theseT1Pos);                          % Find the unique positions of T1 in the stream
    minX_T1 = min(xPosition);                           % Find the earliest position of T1
    maxX_T1 = max(xPosition);                           % Find the latest position of T1
    minErr = 1-maxX_T1-1;                               % Find the minimum possible error and subtract 1
    maxErr = nLetters-minX_T1+1;                        % Find the maximum possible error and add 1
    xDomain = minErr:maxErr;                            % Calculate each integer in the error range
    nObservations = numel(theseT1Error);                % Calculate total number of observations
    
    % Generate the 'pseudo-uniform' distribution, which is the
    % expected distribution of errors if a random guess was
    % provided on every trial. This isn't an actual uniform
    % distribution because the most extreme errors are only
    % possible on trials in which targets appear at their most
    % extreme positions.
    pseudo_uniform = zeros(size(xDomain));

    % Cycle through each possible T1 position.
    for thisPosNo = 1:numel(theseT1Pos)

        % Identify the actual T1 position corresponding to the
        % position number. For example, the first position number
        % might be the 7th position in the stream.
        thisPos = theseT1Pos(thisPosNo);

        % Add to the pseudo-uniform distribution one unit for every
        % possible error given that T1 position.
        pseudo_uniform((1-thisPos-minErr+1):(nLetters-thisPos-minErr+1)) = pseudo_uniform((1-thisPos-minErr+1):(nLetters-thisPos-minErr+1))+ones(1,nLetters);

    end
    
    thisHist = hist(theseT1Error,xDomain);              % Get a histogram of T1 error
    thisHist = thisHist/nObservations;                  % Normalise the histogram
    bar(xDomain,thisHist);                              % Plot this as a bar graph
    hold on;                                            % Don't overwrite the axes
	[singleAIC, singleBIC] = aicbic(-allT1MinNegLogLikelihoods_M1(thisLag), nFreeParameters_Single, nObservations); %#ok<ASGLU> % Calculate BIC for the single-episode model
	[dualAIC, dualBIC] = aicbic(-allT1MinNegLogLikelihoods_M2(thisLag), nFreeParameters_Dual, nObservations);       %#ok<ASGLU> % Calculate BIC for the dual-episode model
    xGrid = linspace(minErr,maxErr,gridPoints);                % Grid of x-points at which to calculate the pdf of the model

    bestEstimates = allT1Estimates_M1(thisLag,:);
    yGrid = pdf_normmixture_single(xGrid, bestEstimates(1), bestEstimates(2), bestEstimates(3));
    
    plot(xGrid,yGrid);                                  % Plot the model
    
    % Tidy up the figure
    axis square;
    axis(plotAxes);
    set(gca,'TickDir','out','YMinorTick','on','XMinorTick','on');
    colormap([0 1 1]);
    
    
    % Plot T2
    figure('Color','white','Name',[sampleName ': Lag ' num2str(listLags(thisLag)) ' (T2)']); % Create a figure
    
    % Calculate the domain of possible errors (xDomain).
    xPosition = unique(theseT2Pos);                          % Find the unique positions of T2 in the stream
    minX_T2 = min(xPosition);                           % Find the earliest position of T2
    maxX_T2 = max(xPosition);                           % Find the latest position of T2
    minErr = 1-maxX_T2-1;                               % Find the minimum possible error and subtract 1
    maxErr = nLetters-minX_T2+1;                        % Find the maximum possible error and add 1
    xDomain = minErr:maxErr;                            % Calculate each integer in the error range
    nObservations = numel(theseT2Error);                % Calculate total number of observations

    % Generate the 'pseudo-uniform' distribution, which is the
    % expected distribution of errors if a random guess was
    % provided on every trial.
    pseudo_uniform = zeros(size(xDomain));

    % Cycle through each possible T2 position.
    for thisPosNo = 1:numel(theseT2Pos)

        % Identify the actual T2 position corresponding to the
        % position number. For example, the first position number
        % might be the 7th position in the stream.
        thisPos = theseT2Pos(thisPosNo);

        % Add to the pseudo-uniform distribution one unit for every
        % possible error given that T2 position.
        pseudo_uniform((1-thisPos-minErr+1):(nLetters-thisPos-minErr+1)) = pseudo_uniform((1-thisPos-minErr+1):(nLetters-thisPos-minErr+1))+ones(1,nLetters);

    end
    
    thisHist = hist(theseT2Error,xDomain);              % Get a histogram of T2 error
    thisHist = thisHist/nObservations;                  % Normalise the histogram
    bar(xDomain,thisHist);                              % Plot this as a bar graph
    hold on;                                            % Don't overwrite the axes
	[singleAIC, singleBIC] = aicbic(-allT1T2MinNegLogLikelihoods_M1(thisLag), nFreeParameters_Single, 2*nObservations); %#ok<ASGLU> % Calculate BIC for the single-episode model
	[dualAIC, dualBIC] = aicbic(-allT1T2MinNegLogLikelihoods_M2(thisLag), nFreeParameters_Dual, 2*nObservations);       %#ok<ASGLU> % Calculate BIC for the dual-episode model
    xGrid = linspace(minErr,maxErr,gridPoints);                % Grid of x-points at which to calculate the pdf of the model
    
    if dualBIC > singleBIC
        % Use M1
        bestEstimates = allT2Estimates_M1(thisLag,:);
        yGrid = pdf_normmixture_single(xGrid, bestEstimates(1), bestEstimates(2), bestEstimates(3));
    else
        % Use M2
        bestEstimates = allT2Estimates_M2(thisLag,:);
        p1 = bestEstimates(1) + bestEstimates(4);       % Reconstruct efficacy estimates
        p2 = bestEstimates(4)/p1;
        yGrid = pdf_normmixture(xGrid, p1, bestEstimates(2), bestEstimates(3), p2, bestEstimates(5), bestEstimates(6));
    end
    
    plot(xGrid,yGrid);                                  % Plot the model
    
    % Tidy up the figure
    axis square;
    axis(plotAxes);
    set(gca,'TickDir','out','YMinorTick','on','XMinorTick','on');
    colormap([1 0 0]);
    
    % Plot T1+T2
    figure('Color','white','Name',[sampleName ': Lag ' num2str(listLags(thisLag)) ' (T1+T2)']);
    
    % Combine the distributions on a common scale (i.e. alter T1
    % errors to reflect position relative to T2.
    theseT1T2Error = [theseT1Error-listLags(thisLag); theseT2Error];         
    minT1T2Error = min(theseT1T2Error);
    maxT1T2Error = max(theseT1T2Error);
    
    % Calculate the domain of possible errors (xDomain).
    theseT1T2Pos = theseT2Pos;                          % Errors are relative to T2
    xPosition = unique(theseT1T2Pos);                   % Find the unique positions of T2 in the stream
    minX_T1T2 = min(xPosition);                         % Find the earliest position of T2
    maxX_T1T2 = max(xPosition);                         % Find the latest position of T2
    minErr = 1-maxX_T1T2-1;                             % Find the minimum possible error and subtract 1
    maxErr = nLetters-minX_T1T2+1;                      % Find the maximum possible error and add 1
    xDomain = minErr:maxErr;                            % Calculate each integer in the error range
    nObservations = numel(theseT1T2Error);              % Calculate total number of observations

    % Generate the 'pseudo-uniform' distribution, which is the
    % expected distribution of errors if a random guess was
    % provided on every trial.
    pseudo_uniform = zeros(size(xDomain));

    % Cycle through each possible T2 position.
    for thisPosNo = 1:numel(theseT2Pos)

        % Identify the actual T2 position corresponding to the
        % position number. For example, the first position number
        % might be the 7th position in the stream.
        thisPos = theseT2Pos(thisPosNo);

        % Add to the pseudo-uniform distribution one unit for every
        % possible error given that T2 position.
        pseudo_uniform((1-thisPos-minErr+1):(nLetters-thisPos-minErr+1)) = pseudo_uniform((1-thisPos-minErr+1):(nLetters-thisPos-minErr+1))+ones(1,nLetters);

    end
    
    thisHistT1 = hist(theseT1Error-listLags(thisLag),xDomain);      % Get a histogram of T1 error
    thisHistT1 = thisHistT1/(nObservations/2);                      % Normalise the histogram
    thisHistT2 = hist(theseT2Error,xDomain);                        % Get a histogram of T2 error
    thisHistT2 = thisHistT2/(nObservations/2);                      % Normalise the histogram
    
    bar(xDomain,[thisHistT1; thisHistT2]','stacked');               % Plot this as a bar graph
    hold on;                                                        % Don't overwrite the axes
	[singleAIC, singleBIC] = aicbic(-allT1T2MinNegLogLikelihoods_M1(thisLag), nFreeParameters_Single, 2*nObservations); % Calculate BIC for the single-episode model
	[dualAIC, dualBIC] = aicbic(-allT1T2MinNegLogLikelihoods_M2(thisLag), nFreeParameters_Dual, 2*nObservations);       % Calculate BIC for the dual-episode model
    xGrid = linspace(minErr,maxErr,gridPoints);                     % Grid of x-points at which to calculate the pdf of the model
    
    if dualBIC > singleBIC
        % Use M1
        bestEstimates = allT1T2Estimates_M1(thisLag,:);
        yGrid = pdf_normmixture_single(xGrid, bestEstimates(1), bestEstimates(2), bestEstimates(3));
    else
        % Use M2
        bestEstimates = allT1T2Estimates_M2(thisLag,:);
        p1 = bestEstimates(1) + bestEstimates(4);       % Reconstruct efficacy estimates
        p2 = bestEstimates(4)/p1;
        yGrid = pdf_normmixture(xGrid, p1, bestEstimates(2), bestEstimates(3), p2, bestEstimates(5), bestEstimates(6));
    end
    
    plot(xGrid,2*yGrid);                                % Plot the model
    
    % Tidy up the figure
    axis square;
    axis(plotAxes);
    set(gca,'TickDir','out','YMinorTick','on','XMinorTick','on');
    colormap([0 1 1; 1 0 0]);
    
end