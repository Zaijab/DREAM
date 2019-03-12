% Clear memory
% Clear all;
clc; close all; clear *;

global DREAM_dir EXAMPLE_dir

% Store working directory and subdirectory containing the files needed to run this example
DREAM_dir = pwd; EXAMPLE_dir = [pwd '\example_' num2str(example)];

% Add subdirectory to search path
addpath(EXAMPLE_dir);

% Recommended parameter settings
MCMCPar.seq = 3;                        % Number of Markov chains / sequences
MCMCPar.DEpairs = 1;                    % Number of chain pairs to generate candidate points
MCMCPar.nCR = 3;                        % Number of crossover values used
MCMCPar.k = 10;                         % Thinning parameter for appending X to Z
MCMCPar.parallelUpdate = 0.9;           % Fraction of parallel direction updates
MCMCPar.eps = 5e-2;                     % Perturbation for ergodicity
MCMCPar.steps = inline('MCMCPar.ndraw/(20 * MCMCPar.seq)'); % Steps before calculating convergence diagnostics
MCMCPar.m0 = inline('10 * MCMCPar.n');  % Initial size of matrix Z
MCMCPar.pJumpRate_one = 0.20;           % Probability of selecting a jumprate of 1 --> jump between modes
MCMCPar.pCR = 'Yes';                    % Adaptive tuning of crossover values (Yes or No)
MCMCPar.Restart = 'No';                 % Restart run (Yes or No)
MCMCPar.modout = 'Yes';                 % Return model (function) simulations of samples Yes or No)?
MCMCPar.save = 'Yes';                   % Save output during the run (Yes or No)
MCMCPar.ABC = 'No';                     % Approximate Bayesian Computation or Not?

% -------------------------------------------------------------------------
% If MCMCPar.modout = 'Yes' --> the simulations of the model are stored.
% But this only happens if calibration data vector exists!!
% -------------------------------------------------------------------------

% n-dimensional banana shaped Gaussian distribution

% ---------------------------- Check the following 2 papers ------------------------------- %
%                                                                                           %
%   Vrugt, J.A., C.J.F. ter Braak, C.G.H. Diks, D. Higdon, B.A. Robinson, and J.M. Hyman,   %
%       Accelerating Markov chain Monte Carlo simulation by differential evolution with     %
%       self-adaptive randomized subspace sampling, International Journal of Nonlinear      %
%       Sciences and Numerical Simulation, 10(3), 271-288, 2009.                            %
%                                                                                           %
%   Vrugt, J.A., H.V. Gupta, W. Bouten and S. Sorooshian, A Shuffled Complex Evolution      %
%       Metropolis algorithm for optimization and uncertainty assessment of hydrologic      %
%       model parameters, Water Resour. Res., 39 (8), 1201, doi:10.1029/2002WR001642, 2003. %
%                                                                                           %
% ----------------------------------------------------------------------------------------- %

% Application specific parameter settings
MCMCPar.n = 10;                         % Dimension of the problem (number of parameters to be estimated)
MCMCPar.ndraw = 50000;                  % Maximum number of function evaluations
MCMCPar.T = 1;                          % Each Tth sample is collected in the chains
MCMCPar.prior = 'COV';                  % Latin Hypercube sampling (options, "LHS", "COV" and "PRIOR")
MCMCPar.BoundHandling = 'None';         % Boundary handling (options, "Reflect", "Bound", "Fold", and "None");
MCMCPar.lik = 4;                        % Define the likelihood function --> log-density from model

% Define modelName
ModelName = 'Banana_func';

% Provide information to do initial sampling ("COV")
MCMCPar.mu = zeros(1,MCMCPar.n);        % Provide mean of initial sample
MCMCPar.cov = 5 * eye(MCMCPar.n);       % Initial covariance

% Define the specific properties of the banana function
Extra.cov  = eye(MCMCPar.n); Extra.cov(1,1) = 100;      % Target covariance
Extra.invC = inv(Extra.cov);                            % Inverse of target covariance
Extra.b = 0.1;                                          % For the function

% Run the DREAM_ZS algorithm
[Sequences,X,Z,output,fx] = dream_zs(MCMCPar,ModelName,Extra);

% Create a single matrix with values sampled by chains
ParSet = GenParSet(Sequences,MCMCPar);

%Save output for postprocessing
save GWF1P81R05.mat
