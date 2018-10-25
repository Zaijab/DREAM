function [par_unc,tot_unc] = predMCMC(Pars,RMSE,ModelName,Measurement,MCMCPar,Extra,PredInt,fx);
% Calculates the desired parameter and total prediction uncertainty from the MCMC samples
%
% SYNOPSIS:     [par_unc,tot_unc] = predMCMC(ParSet,RMSE,ModelName,Measurement,MCMCPar,Extra);
%               [par_unc,tot_unc] = predMCMC(ParSet,RMSE,ModelName,Measurement,MCMCPar,Extra,PredInt);
%
% Input:        ParSet = matrix with MCMC samples
%               RMSE = Root Mean Square Error of best solution
%               ModelName = name of the model to be executed
%               Measurement = structure with measurement variables
%               MCMCPar = structure with MCMC variables
%               Extra = structure with additional variables needed for model
%               fx = (optional argument) with model simulations
%               PredInt = (optional) argument with % prediction uncertainty range (e.g. 95%)
%
% Output:       par_unc = matrix with lower + upper bound simulation ranges of parmeter uncertainty
%               tot_unc = matrix with lower + upper bound simulation ranges of total uncertainty

% If PredInt is not defined, automatically assume 95% uncertainty ranges
if nargin < 7,
    PredInt = 95;
end;

% Determine whether fx exists or not
if nargin < 8,
    fx = [];
end;

% How many samples?
N_Pars = size(Pars,1);

% Check whether "fx" exist or not -- if so, then we define Spar directly
if ( isempty(fx) == 0 ),

    % Set Spar to be fx
    Spar = fx(end - N_Pars + 1:end,1:Measurement.N)'; % --> transpose so that each simulation is vertical

else

    % Pre-allocate memory for simulation results
    Spar = zeros(Measurement.N,N_Pars);;

    % Initialize waitbar
    h = waitbar(0,'Running posterior simulations - Please wait...');

    for j = 1:N_Pars,
        % Define evaluation string
        evalstr = ['Spar(:,j) = ',ModelName,'(Pars(j,1:MCMCPar.n),Extra);'];
        % Call model to generate simulated data
        eval(evalstr); 
        % Update the waitbar --> to see simulation progress on screen
        waitbar(j/N_Pars,h);
    end;

    % Close the waitbar
    close(h);

end;

% Define the lower bound of the prediction interval
Lb = floor( (100 - PredInt)/200 * N_Pars);

% Define the upper bound of the prediction interval
Ub = N_Pars - Lb;

% Now add the RMSE to create the total uncertainty (homoscedastic error!)
Smod = Spar + normrnd(0,RMSE,Measurement.N,N_Pars);

% Now sort to get desired ranges
for j = 1:Measurement.N,
    % Sort the model output from the jth parameter combination from low to high
    a = sort(Spar(j,1:N_Pars));
    % And take the desired prediction uncertainty ranges
    par_unc(j,1) = a(Lb); par_unc(j,2) = a(Ub);
    % Same with total uncertainty
    a = sort(Smod(j,1:N_Pars));
    % And take the desired prediction uncertainty ranges
    tot_unc(j,1) = a(Lb); tot_unc(j,2) = a(Ub);
end;