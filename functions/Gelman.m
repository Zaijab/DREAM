function [R_stat] = Gelman(Sequences,MCMCPar)
% Calculates the R-statistic convergence diagnostic
% ----------------------------------------------------
% For more information please refer to: Gelman, A. and D.R. Rubin, 1992. 
% Inference from Iterative Simulation Using Multiple Sequences, 
% Statistical Science, Volume 7, Issue 4, 457-472.
%
% Written by Jasper A. Vrugt
% Los Alamos, August 2007
% ----------------------------------------------------

% Compute the dimensions of Sequences
[n,nrY,m] = size(Sequences);

if (n < 10),
    % Set the R-statistic to a large value
    R_stat = -2 * ones(1,MCMCPar.n);
else
    % Step 1: Determine the sequence means
    meanSeq = mean(Sequences); meanSeq = reshape(meanSeq(:),MCMCPar.n,m)';
    
    % Step 1: Determine the variance between the sequence means 
    B = n * var(meanSeq);
    
    % Step 2: Compute the variance of the various sequences
    for zz = 1:MCMCPar.seq,
        varSeq(zz,:) = var(Sequences(:,:,zz));
    end;
    
    % Step 2: Calculate the average of the within sequence variances
    W = mean(varSeq);
    
    % Step 3: Estimate the target mean
    %mu = mean(meanSeq);
    
    % Step 4: Estimate the target variance
    sigma2 = ((n - 1)/n) * W + (1/n) * B;
    
    % Step 5: Compute the R-statistic
    R_stat = sqrt((m + 1)/m * sigma2 ./ W - (n-1)/m/n);
    
end;