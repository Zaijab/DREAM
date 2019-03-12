% ----------------------------------------------------------------------------------------------%
%                                                                                               %
% DDDDDDDDDDDDDDD    RRRRRRRRRRRRRR     EEEEEEEEEEEEEEEE       AAAAA       MMM             MMM  %
% DDDDDDDDDDDDDDDD   RRRRRRRRRRRRRRR    EEEEEEEEEEEEEEEE       AAAAA       MMMM           MMMM  %
% DDD          DDD   RRR          RRR   EEE                   AAA AAA      MMMMM         MMMMM  %
% DDD          DDD   RRR          RRR   EEE                   AAA AAA      MMMMMM       MMMMMM  %
% DDD          DDD   RRR          RRR   EEE                  AAA   AAA     MMM MMM     MMM MMM  %
% DDD          DDD   RRR          RRR   EEE                  AAA   AAA     MMM  MMM   MMM  MMM  %
% DDD          DDD   RRRRRRRRRRRRRRRR   EEEEEEEEEEEEEEEE    AAA     AAA    MMM   MMM MMM   MMM  %
% DDD          DDD   RRRRRRRRRRRRRRRR   EEEEEEEEEEEEEEEE    AAAAAAAAAAA    MMM    MMMM     MMM  %
% DDD          DDD   RRR          RRR   EEE                AAA       AAA   MMM             MMM  %
% DDD          DDD   RRR          RRR   EEE                AAA       AAA   MMM             MMM  %
% DDD          DDD   RRR          RRR   EEE               AAA         AAA  MMM             MMM  %
% DDD          DDD   RRR          RRR   EEE               AAA         AAA  MMM             MMM  %
% DDDDDDDDDDDDDDDD   RRR          RRR   EEEEEEEEEEEEEEEE  AAA         AAA  MMM             MMM  %
% DDDDDDDDDDDDDDD    RRR          RRR   EEEEEEEEEEEEEEEE  AAA         AAA  MMM             MMM  %
%                                                                                               %
% ----------------------------------------------------------------------------------------------%

% ------------- DREAM with sampling from past and snooker updates: DREAM_ZS --------------------%
%                                                                                               %
% The code presented herein is a Markov Chain Monte Carlo algorithm that runs multiple chains   %
% in parallel for efficient posterior exploration. The algorithm, entitled DREAM_(ZS) is        %
% based on the original DREAM sampling scheme, but uses sampling from an archive of past        %
% states to generate candidate points in each individual chain. Theoy and numerical examples of %
% DREAM_(ZS) have been presented in Vrugt et al. (2009). Details can also be found in           %
% Ter Braak and Vrugt (2008)                                                                    %
%                                                                                               %
% Sampling from past has three main advantages:                                                 %
% (1) Circumvents the requirement of using N = d for posterior exploration. This will speed-up  %
% convergence to a limiting distribution, especially for high-dimensional problems (large d).   %
% (2) Outlier chains do not need explicit consideration. By sampling historical states,         %
% aberrant trajectories an jump directly to the modal region at any time during the             %
% simulation. The N path ways simulated with DREAM_(ZS) therefore maintain detailed balance at  %
% every singe step in the chain.                                                                %
% (3) The transition kernel defining the jumps in each of the chains does not require           %
% information about the current states of the chains. This is of great advantage in a           %
% multi-processor environment where the N candidate points can be generated simultaneously so   %
% that each chain can evolve most efficiently on a different computer. Details of this will be  %
% given in a later publication, which should be ready within the next few months.               %
%                                                                                               %
% DREAM_(ZS) also contains a snooker updater to maximize the diversity of candidate points      %
% and generate jumps beyond parallel direction updates. Finally, DREAM_(ZS) contains subspace   %
% learning in a similar way as DREAM, to maximize the squared jumping distance between two      %
% subsequent points in each chain. This idea has been presented in Vrugt et al. (2008) and      %
% shown to significantly increase the efficiency of posterior exploration. All these options    %
% can be activated from the input file.                                                         %
%                                                                                               %
% DREAM_(ZS) developed by Jasper A. Vrugt and Cajo ter Braak                                    %
%                                                                                               %
% --------------------------------------------------------------------------------------------- % 
%                                                                                               %
% SYNOPSIS: [Sequences,X,Z,output] = dream_zs(MCMCPar,ModelName)                                %
%           [Sequences,X,Z,output] = dream_zs(MCMCPar,ModelName,Extra)                          % 
%           [Sequences,X,Z,output] = dream_zs(MCMCPar,ModelName,Extra,ParRange)                 %
%           [Sequences,X,Z,output] = dream_zs(MCMCPar,ModelName,Extra,ParRange,Measurement)     %
%                                                                                               %
% Input:    MCMCPar = structure with DREAM parameters                                           %
%           ModelName = name of the function                                                    %
%           Extra = optional structure with arguments to be passed to function                  %
%           ParRange = optional structure with parameter ranges                                 %
%           Measurement = optional structure with measurement information                       %
%                                                                                               %
% Output:   Sequences = 3D array with Markov chain evolution                                    %
%           X = final position of chains and correponding density values                        %
%           Z = matrix with thinned sample history                                              %
%           output = structure with convergence properties, acceptance rate, CR values, etc.    %
%                                                                                               %
% The directory \PostProcessing contains a script "PostProcMCMC" that will compute various      %
% posterior statistics (MEAN, STD, MAP, CORR) and create create various plots including,        % 
% marginal posterior parameter distributions, R_stat convergence diagnostic, two-dimensional    %
% correlation plots of the posterior parameter samples, chain convergence plots, and parameter  % 
% and total posterior simulation uncertainty ranges (interval can be specified)                 % 
%                                                                                               %
% --------------------------------------------------------------------------------------------- % 
%                                                                                               %
% This algorithm has been described in:                                                         %
%                                                                                               %
%   Laloy, E., and J.A. Vrugt, High-dimensional posterior exploration of hydrologic models      %
%       using multiple-try DREAM_(ZS) and high-performance computing, Water Resources Research, %
%       48, W01526, doi:10.1029/2011WR010608, 2012.                                             %
%                                                                                               %
%   ter Braak, C.J.F., and J.A. Vrugt, Differential Evolution Markov Chain with snooker updater %
%       and fewer chains, Statistics and Computing, 10.1007/s11222-008-9104-9, 2008             %
%                                                                                               %
%   Vrugt, J.A., E. Laloy, and C.J.F. ter Braak, DiffeRential Evolution Adaptive Metropolis     %
%       with Sampling from the Past and Subspace Updating, SIAM journal on Optimization         %
%                                                                                               %
% For more information please read:                                                             %
%                                                                                               %
%   Vrugt J.A., H.V. Gupta, W. Bouten and S. Sorooshian, A Shuffled Complex Evolution           %
%       Metropolis algorithm for optimization and uncertainty assessment of hydrologic model    %
%       parameters, Water Resour. Res., 39 (8), 1201, doi:10.1029/2002WR001642, 2003.           %
%                                                                                               %
%   ter Braak, C.J.F., A Markov Chain Monte Carlo version of the genetic algorithm Differential %
%       Evolution: easy Bayesian computing for real parameter spaces, Stat. Comput., 16,        %
%       239 - 249, doi:10.1007/s11222-006-8769-1, 2006.                                         %
%                                                                                               %
%   Vrugt, J.A., C.J.F. ter Braak, M.P. Clark, J.M. Hyman, and B.A. Robinson, Treatment of      %
%       input uncertainty in hydrologic modeling: Doing hydrology backward using Markov         %
%       chain Monte Carlo, Water Resour. Res., 44, W00B09, doi:10.1029/2007WR006720, 2008.      %
%                                                                                               %
%   Vrugt, J.A., C.J.F. ter Braak, C.G.H. Diks, D. Higdon, B.A. Robinson, and J.M. Hyman,       %
%       Accelerating Markov chain Monte Carlo simulation by differential evolution with         %
%       self-adaptive randomized subspace sampling, International Journal of Nonlinear Sciences %
%       and Numerical Simulation}, 10(3), 273-290, 2009.                                        %
%                                                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                               %
%     Copyright (C) 2011-2012  the authors                                                      %
%                                                                                               %
%     This program is free software: you can modify it under the terms of the GNU General       %
%     Public License as published by the Free Software Foundation, either version 3 of the      %
%     License, or (at your option) any later version.                                           %
%                                                                                               %
%     This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; %
%     without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. %
%     See the GNU General Public License for more details.                                      %
%                                                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                               %
% Written by Jasper A. Vrugt: jasper@uci.edu                                                    %
%                                                                                               %
% Version 0.5: January 2009                                                                     %
% Version 1.0: April 2011         Maintenance update, explicit treatment of prior distribution  %
% Version 1.1: August 2011        Whittle likelihood function (SPECTRAL ANALYSIS !!)            %
% Version 1.2: April 2012         Simplified code (removed variables) + graphical interface     %
% Version 1.3: June 2012          Simulations stored, new example, and updated likelihood func. %
% Version 1.4: January 2013       Simplification of metrop.m and dream_zs.m                     %
%                                                                                               %
% --------------------------------------------------------------------------------------------- %

% Different test examples
% example 1:  n-dimensional banana shaped Gaussian distribution
% example 2:  n-dimensional Gaussian distribution
% example 3:  n-dimensional multimodal mixture distribution
% example 4:  real-world example using hymod rainfall - runoff model (HYMOD code in MATLAB)
% example 4.1: GW model
% example 5:  real-world example using hymod rainfall - runoff model (HYMOD code in FORTRAN)
% example 6:  rainfall-runoff model with generalized log-likelihood function
% example 7:  HYDRUS-1D soil hydraulic model: using prior information on soil hydraulic parameters
% example 8:  multivariate student t distribution
% example 9:  Rainfall-runoff model with Whittle's likelihood function
% example 10: the use of prior information in a multimodel mixture distrbibution
% example 11: multivariate student t distribution
% example 12: pedometrics problem involving variogram fitting
% example 13: Nash-Cascade example --> heteroscedastic errors
% example 14: ABC inference for hydrologic model
% example 15: ABC inference using 10 bivariate normal distributions
% example 16: Hydrogeophysics example

% Clear memory
% Clear all; 
clc; close all; clear *;

% Which example to run?
example = 4.1;

global DREAM_dir EXAMPLE_dir

% Store working directory and subdirectory containing the files needed to run this example
DREAM_dir = pwd; EXAMPLE_dir = [pwd '\example_' num2str(example)];

% Add subdirectory to search path
addpath(EXAMPLE_dir);

% Recommended parameter settings
MCMCPar.seq = 3;                        % Number of Markov chains / sequences (for high dimensional and highly nonlinear problems, larger values work beter!!)
MCMCPar.DEpairs = 1;                    % Number of chain pairs to generate candidate points
MCMCPar.nCR = 3;                        % Number of crossover values used
MCMCPar.k = 10;                         % Thinning parameter for appending X to Z
MCMCPar.parallelUpdate = 0.9;           % Fraction of parallel direction updates
MCMCPar.eps = 5e-2;                     % Perturbation for ergodicity
MCMCPar.steps = inline('MCMCPar.ndraw/(20 * MCMCPar.seq)'); % Number of steps before calculating convergence diagnostics
MCMCPar.m0 = inline('10 * MCMCPar.n');  % Initial size of matrix Z
MCMCPar.pJumpRate_one = 0.20;           % Probability of selecting a jumprate of 1 --> jump between modes
MCMCPar.pCR = 'Yes';                    % Adaptive tuning of crossover values (Yes or No)
MCMCPar.Restart = 'No';                 % Restart run (Yes or No)
MCMCPar.modout = 'Yes';                 % Return model (function) simulations of samples Yes or No)?
MCMCPar.save = 'Yes';                    % Save output during the run (Yes or No)
MCMCPar.ABC = 'No';                     % Approximate Bayesian Computation or Not?

% -------------------------------------------------------------------------
% If MCMCPar.modout = 'Yes' --> the simulations of the model are stored. 
% But this only happens if calibration data vector exists!!
% -------------------------------------------------------------------------

if example == 4.1,    % GW model 81 parameters
    
    % Problem specific parameter settings
    MCMCPar.seq =5;                                 %Number of Markov chains / sequences (for high dimensional and highly nonlinear problems, larger values work beter!!)
    MCMCPar.n = 81;                                 %Dimension of the problem (number of parameters to be estimated)
    MCMCPar.ndraw = 5e6;                            %Maximum number of function evaluations
    MCMCPar.T = ceil(max(1,MCMCPar.ndraw/1.5e6));   %Each Tth sample is collected in the chains
    MCMCPar.prior = 'LHS';                          %Latin Hypercube sampling (options, "LHS", "COV" and "PRIOR")
    MCMCPar.BoundHandling = 'Bound';                %Boundary handling (options, "Reflect", "Bound", "Fold", and "None"); 
    MCMCPar.modout = 'No';                          %Return model (function) simulations of samples Yes or No)?
    MCMCPar.lik = 3;                                %Define likelihood function -- Sum of Squared Error
    MCMCPar.Best =inf;                              %Best fitness (sum squared error) so far
    % Define modelName
    ModelName = 'Case01';

    % Give the parameter ranges (minimum and maximum values)
    ParRange.minn = 1e-3*ones(1,MCMCPar.n);
    ParRange.maxn = 1e-1*ones(1,MCMCPar.n);

    % Load the data
    load TrueHead.mat
    
    % Define the measured streamflow data
    Measurement.MeasData = TrueHead(:)+e(:);

    % Run the DREAM_ZS algorithm
    [Sequences,X,Z,output,fx] = dream_zs(MCMCPar,ModelName,[],ParRange,Measurement);
end;



if example == 1, % n-dimensional banana shaped Gaussian distribution

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

end;

if example == 2,    % n-dimensional Gaussian distribution

    % ---------------------------- Check the following 2 papers ------------------------------- %
    %                                                                                           %
    %   Vrugt, J.A., C.J.F. ter Braak, C.G.H. Diks, D. Higdon, B.A. Robinson, and J.M. Hyman,   %
    %       Accelerating Markov chain Monte Carlo simulation by differential evolution with     %
    %       self-adaptive randomized subspace sampling, International Journal of Nonlinear      %
    %       Sciences and Numerical Simulation, 10(3), 271-288, 2009.                            %
    %                                                                                           %
    %   Ter Braak, C.J.F., and J.A. Vrugt, Differential Evolution Markov Chain with snooker     %
    %       updater and fewer chains, Statistics and Computing, 10.1007/s11222-008-9104-9,      %
    %       2008.                                                                               %
    %                                                                                           %
    % ----------------------------------------------------------------------------------------- %

    % Problem specific parameter settings
    MCMCPar.n = 2;                          % Dimension of the problem (number of parameters to be estimated)
    MCMCPar.ndraw = 5e4;                    % Maximum number of function evaluations
    MCMCPar.T = 1;                          % Each Tth sample is collected in the chains (This parameter was not part of MS)
    MCMCPar.prior = 'LHS';                  % Latin Hypercube sampling (options, "LHS", "COV" and "PRIOR")
    MCMCPar.BoundHandling = 'None';         % Boundary handling (options, "Reflect", "Bound", "Fold", and "None");
    MCMCPar.lik = 4;                        % Model returns log-density

    % Define modelName
    ModelName = 'normalfunc';

    % Give the parameter ranges (minimum and maximum values) (mostly overdispersed)
    ParRange.minn = [-5 * ones(1,MCMCPar.n)]; ParRange.maxn = [15 * ones(1,MCMCPar.n)];

    % Give the parameter ranges (minimum and maximum values) (underdispersed)
    % ParRange.minn = [9.9 * ones(1,MCMCPar.n)]; ParRange.maxn = [10 * ones(1,MCMCPar.n)];

    % ------ Define covariance matrix of target distribution --------------

%     % Construct the d x d covariance matrix
%     A = 0.5 * eye(MCMCPar.n) + 0.5 * ones(MCMCPar.n);
%     % Rescale to variance-covariance matrix of interest
%     for i = 1:MCMCPar.n
%         for j = 1:MCMCPar.n
%             C(i,j) = A(i,j) * sqrt(i * j);
%         end
%     end

     C = 1.1 * eye(MCMCPar.n) + 0.9 * ones(MCMCPar.n);


    % Store inverse of covariance for "normalfunc"
    Extra.invC = inv(C);
    
    % ---------------------------------------------------------------------

    % Run the DREAM_ZS algorithm   
    [Sequences,X,Z,output,fx] = dream_zs(MCMCPar,ModelName,Extra,ParRange);

end;

if example == 3,    % n-dimensional multimodal mixture distribution

    % ---------------------------- Check the following two papers ----------------------------- %
    %                                                                                           %
    %   Vrugt, J.A., C.J.F. ter Braak, C.G.H. Diks, D. Higdon, B.A. Robinson, and J.M. Hyman,   %
    %       Accelerating Markov chain Monte Carlo simulation by differential evolution with     %
    %       self-adaptive randomized subspace sampling, International Journal of Nonlinear      %
    %       Sciences and Numerical Simulation, 10(3), 271-288, 2009.                            %
    %                                                                                           %
    %   Ter Braak, C.J.F., and J.A. Vrugt, Differential Evolution Markov Chain with snooker     %
    %       updater and fewer chains, Statistics and Computing, 10.1007/s11222-008-9104-9,      %
    %       2008.                                                                               %
    %                                                                                           %
    % ----------------------------------------------------------------------------------------- %

    % Problem specific parameter settings
    MCMCPar.n = 10;                         % Dimension of the problem (number of parameters to be estimated)
    MCMCPar.ndraw = 200000;                 % Maximum number of function evaluations
    MCMCPar.T = 10;                         % Each Tth sample is collected in the chains
    MCMCPar.prior = 'COV';                  % Latin Hypercube sampling (options, "LHS", "COV" and "PRIOR")
    MCMCPar.BoundHandling = 'None';         % Boundary handling (options, "Reflect", "Bound", "Fold", and "None");
    MCMCPar.lik = 1;                        % Model returns directly the probability density

    % Define modelName
    ModelName = 'mixturemodel';

    % Provide information to do covariance sampling ("COV")
    MCMCPar.mu = zeros(1,MCMCPar.n);        % Provide mean of initial sample
    MCMCPar.cov = eye(MCMCPar.n);           % Initial covariance

    % Run the DREAM_ZS algorithm   
    [Sequences,X,Z,output,fx] = dream_zs(MCMCPar,ModelName);

end;

if example == 4,    % HYMOD rainfall - runoff model (coded in MATLAB)

    % ---------------------------- Check the following 3 papers ------------------------------- %
    %                                                                                           %
    %   Vrugt, J.A., C.J.F. ter Braak, M.P. Clark, J.M. Hyman, and B.A. Robinson, Treatment of  %
    %      input uncertainty in hydrologic modeling: Doing hydrology backward with Markov chain %
    %      Monte Carlo simulation, Water Resources Research, 44, W00B09,                        %
    %      doi:10.1029/2007WR006720, 2008.                                                      %
    %                                                                                           %
    %   Vrugt, J.A., C.J.F. ter Braak, H.V. Gupta, and B.A. Robinson, Equifinality of formal    %
    %       (DREAM) and informal (GLUE) Bayesian approaches in hydrologic modeling?, Stochastic %
    %       Environmental Research and Risk Assessment, 23(7), 1011-1026, 				        %
    %       doi:10.1007/s00477-008-0274-y, 2009                                                 %
    %                                                                                           %
    %   Vrugt, J.A., H.V. Gupta, W. Bouten and S. Sorooshian, A Shuffled Complex Evolution      %
    %       Metropolis algorithm for optimization and uncertainty assessment of hydrologic      %
    %       model parameters, Water Resour. Res., 39 (8), 1201, doi:10.1029/2002WR001642, 2003. %
    %                                                                                           %
    % ----------------------------------------------------------------------------------------- %

    % Problem specific parameter settings
    MCMCPar.n = 5;                          % Dimension of the problem (number of parameters to be estimated)
    MCMCPar.ndraw = 5e3;                  % Maximum number of function evaluations
    MCMCPar.T = 1;                          % Each Tth sample is collected in the chains
    MCMCPar.prior = 'LHS';                  % Latin Hypercube sampling (options, "LHS", "COV" and "PRIOR")
    MCMCPar.BoundHandling = 'Reflect';      % Boundary handling (options, "Reflect", "Bound", "Fold", and "None"); 
    MCMCPar.modout = 'Yes';                 % Return model (function) simulations of samples Yes or No)?
    MCMCPar.lik = 3;                        % Define likelihood function -- Sum of Squared Error
    %MCMCPar.lik == 1, Model returns posterior density
    %MCMCPar.lik == 2, Log-density function
    %MCMCPar.lik == 3, Model returns vector of predictions
    %MCMCPar.lik == 4, Model returns log posterior density
    %MCMCPar.lik == 8, % Generalized log likelihood (GL)
    

    % Define modelName
    ModelName = 'hymodMATLAB';

    % Give the parameter ranges (minimum and maximum values)
    ParRange.minn = [1.0 0.10 0.10 0.00 0.10]; ParRange.maxn = [500 2.00 0.99 0.10 0.99];

    % Load the Leaf River data
    load bound.txt;
    % Then read the boundary conditions -- only do two years
    Extra.MaxT = 795;
    % Define the PET, Measured Streamflow and Precipitation.
    Extra.E = bound(1:Extra.MaxT,5); Extra.P = sum(bound(1:Extra.MaxT,6:9),2);
    % Area factor to translate HYMOD output in mm/d to m3/s (calibration data); (area Leaf River is 1944 km2)
    Extra.F = 1944 * (1000 * 1000 ) / (1000 * 60 * 60 * 24);
    % Define the measured streamflow data
    Measurement.MeasData = bound(65:Extra.MaxT,4);
    % We need to specify the Measurement error of the data in Measurement.Sigma
    % With option 3, Measurement.Sigma is integrated out the likelihoon function
    % With any other option, Sigma needs to be defined

    % We can estimate the measurement error directly if we use temporal differencing
    % The function MeasError provides an estimate of error versus flow level
    % out = MeasError(Measurement.MeasData;
    % For the Leaf River watershed this results in a heteroscedastic error
    % that is about 10% of the actual measured discharge value, thus
    % You can check this by plotting out(:,1) versus out(:,2)
    % Measurement.Sigma = 0.1*Measurement.MeasData; % --> option 2

    % Run the DREAM_ZS algorithm
    [Sequences,X,Z,output,fx] = dream_zs(MCMCPar,ModelName,Extra,ParRange,Measurement);
    
end;

if example == 5,    % HYMOD rainfall - runoff model (coded in FORTRAN)

    % Problem specific parameter settings
    MCMCPar.n = 5;                          % Dimension of the problem (number of parameters to be estimated)
    MCMCPar.ndraw = 10000;                  % Maximum number of function evaluations
    MCMCPar.T = 1;                          % Each Tth sample is collected in the chains
    MCMCPar.prior = 'LHS';                  % Latin Hypercube sampling (options, "LHS", "COV" and "PRIOR")
    MCMCPar.BoundHandling = 'Reflect';      % Boundary handling (options, "Reflect", "Bound", "Fold", and "None");
    MCMCPar.lik = 3;                        % Define likelihood function -- Sum of Squared Error

    % Define modelName
    ModelName = 'hymodFORTRAN';

    % Give the parameter ranges (minimum and maximum values)
    ParRange.minn = [1.0 0.10 0.10 0.00 0.10]; ParRange.maxn = [500 2.00 0.99 0.10 0.99];

    % Leaf River data -- forcing conditions not needed --> externally loaded by FORTRAN executable
    load bound.txt;

    % Then read the boundary conditions -- only do two years
    Extra.MaxT = 795;

    % Area factor to translate HYMOD output in mm/d to m3/s (calibration data); (area Leaf River is 1944 km2)
    Extra.F = 1944 * (1000 * 1000 ) / (1000 * 60 * 60 * 24);

    % Define the measured streamflow data -- use to compute likelihood function
    Measurement.MeasData = bound(65:Extra.MaxT,4);

    % Run the DREAM_ZS algorithm
    [Sequences,X,Z,output,fx] = dream_zs(MCMCPar,ModelName,Extra,ParRange,Measurement);
    
end;

if example == 6,    % Rainfall-runoff model with generalized log-likelihood

    % ---------------------------- Check the following 2 papers ------------------------------- %
    %                                                                                           %
    % G. Schoups, J.A. Vrugt, F. Fenicia, and N.C. van de Giesen (2010), Corruption of          %
    %     accuracy and efficiency of Markov Chain Monte Carlo simulation by inaccurate          %
    %     numerical implementation of conceptual hydrologic models, Water Resources             %
    %     Research, 46, W10530, doi:10.1029/2009WR008648.                                       %
    %                                                                                           %
    % G. Schoups, and J.A. Vrugt (2010), A formal likelihood function for parameter and         %
    %     predictive inference of hydrologic models with correlated, heteroscedastic and        %
    %     non-Gaussian errors, Water Resources Research, 46, W10531, doi:10.1029/2009WR008933.  %
    %                                                                                           %
    % ----------------------------------------------------------------------------------------- %

    % Problem specific parameter settings
    MCMCPar.n = 11;                         % Dimension of the problem
    MCMCPar.ndraw = 15000;                  % Maximum number of function evaluations
    MCMCPar.T = 1;                          % Each Tth sample is collected in the chains
    MCMCPar.prior = 'LHS';                  % Latin Hypercube sampling (options, "LHS", "COV" and "PRIOR")
    MCMCPar.BoundHandling = 'Reflect';      % Boundary handling (options, "Reflect", "Bound", "Fold", and "None");
    MCMCPar.lik = 8;                        % Generalized likelihood function

    % Define modelName
    ModelName = 'hmodel';

    % Give the parameter ranges (minimum and maximum values)
    %parno:       1     2     3     4     5     6     7      8    9     10    11   12   13   14   15   16   17   18   19   20  21
    %parname:     fA    Imax  Smax  Qsmax alE   alS   alF    Pf   Kfast Kslow std0 std1 beta xi   mu1  phi1 phi2 phi3 phi4 K   lambda
    Extra.fpar = [1     0     100   10    100   1e-6  1e-6   0    2     70    0.1  0    0    1    0    0    0    0    0    0   1];
    % Minimum parameter values
    parmin =     [1     0     10    0     1e-6  1e-6 -10     0    0     0     0    0   -1    0.1  0    0    0    0    0    0   0.1 ];
    % maximum parameter values
    parmax =     [1     10    1000  100   100   1e-6  10     0    10    150   1    1    1    10   100  1    1    1    1    1   1];
    % Select the parameters to be sampled
    Extra.idx_vpar = [2 3 4 5 7 9 10 11 12 13 16];
    
    % And define the associated Parameter ranges
    ParRange.minn = parmin(Extra.idx_vpar); ParRange.maxn = parmax(Extra.idx_vpar);

    % Load the French Broad data
    daily_data = load('03451500.dly');

    % First two years are warm-up
    Extra.idx = [731:size(daily_data,1)]';

    % Define the PET, Measured Streamflow and Precipitation.
    Extra.Precip = daily_data(:,4); Extra.Ep = daily_data(:,5);

    % Define the measured streamflow data
    Measurement.MeasData = daily_data(Extra.idx,6);

    % Run the DREAM_ZS algorithm
    [Sequences,X,Z,output,fx] = dream_zs(MCMCPar,ModelName,Extra,ParRange,Measurement);
    
end;

if example == 7,	% HYDRUS-1D soil hydraulic model: using prior information on soil hydraulic parameters

    % -------------------------------- Check the following paper ------------------------------ %
    %                                                                                           %
    %   B. Scharnagl, J.A. Vrugt, H. Vereecken, and M. Herbst (2011), Bayesian inverse          %
    %	modeling of soil water dynamics at the field scale: using prior information             %
    %	on soil hydraulic properties, Hydrology and Earth System Sciences.                      %
    %                                                                                           %
    % ----------------------------------------------------------------------------------------- %

    % Problem specific parameter settings
    MCMCPar.n = 7;                          % Dimension of the problem
    MCMCPar.ndraw = 5000;                   % Maximum number of function evaluations
    MCMCPar.T = 1;                          % Each Tth sample is collected in the chains
    MCMCPar.prior = 'PRIOR';                % Latin Hypercube sampling (options, "LHS", "COV" and "PRIOR")
    MCMCPar.BoundHandling = 'Reflect';      % Boundary handling (options, "Reflect", "Bound", "Fold", and "None");
    MCMCPar.lik = 4;                        % The model directly returns the log-density

    % Define model name
    ModelName = 'HYDRUS';

    % Define feasible parameter space (minimum and maximum values)
    %				1		2		3				4			5			6		7
    %				[thetar	thetas	log10(alpha)	log10(n)	log10(Ks)	L		hLB
    ParRange.minn =	[0.0430	0.4090	-2.5528			0.1790		-2.2366		-5.4900	-250];
    ParRange.maxn =	[0.0910 0.4810	-2.0706			0.2670		-0.0800		6.2700	-50];

    % Provide observational data and data needed to modify the initial and boundary conditions
    [Extra] = ProvideData;

    % Specify the prior distributions for the various parameters
    MCMCPar.prior_marginal = {  'normrnd(0.0670,0.0060)',...
                                'normrnd(0.4450,0.0090)',...
                                'normrnd(-2.310,0.0600)',...
                                'normrnd(0.2230,0.0110)',...
                                'normrnd(-1.160,0.2700)',...
                                'normrnd(0.3900,1.4700)',...
                                'unifrnd(-250,-50)'};

    % Run the DREAM_ZS algorithm
    [Sequences,X,Z,output,fx] = dream_zs(MCMCPar,ModelName,Extra,ParRange);

end;

if example == 8,    % Simple 1D mixture distribution -- Approximate Bayesian Computation

    % ---------------------------- Check the following 2 papers ------------------------------- %
    %                                                                                           %
    % M. Sadegh, and J.A. Vrugt (2013), Approximate Bayesian computation using DREAM: Theory,   %
    %     Numerical Implementation and Case Studies, Water Resources Research, In Prep.         %
    %                                                                                           %
    % B.M. Turner, and P.B. Sederberg (2013), Approximate Bayesian computation with             %
    %     differential evolution, Journal of Mathematical Psychology, In Press.                 %
    %                                                                                           %
    % ----------------------------------------------------------------------------------------- %

    % Problem specific parameter settings
    MCMCPar.n = 1;                          % Dimension of the problem
    MCMCPar.ndraw = 50000;                  % Maximum number of function evaluations
    MCMCPar.T = 1;                          % Each Tth sample is collected in the chains
    MCMCPar.prior = 'LHS';                  % Latin Hypercube sampling (options, "LHS", "COV" and "PRIOR")
    MCMCPar.BoundHandling = 'Reflect';      % Boundary handling (options, "Reflect", "Bound", "Fold", and "None");
    MCMCPar.lik = 10;                       % ABC likelihood function (standardized)
    MCMCPar.ABC = 'Yes';                    % Specify that we perform ABC
    MCMCPar.delta = 0.025;                  % Delta of the noisy ABC implementation
    MCMCPar.rho = inline('X - Y');          % Define the distance - this case the difference

    % Define modelName
    ModelName = 'ABC_func';

    % Give the parameter ranges (minimum and maximum values)
    ParRange.minn = [-10]; ParRange.maxn = [10];

    % Define Measurement.MeasData --> "Y" in paper
    Measurement.MeasData = 0;
    
    % Run the DREAM_ZS algorithm
    [Sequences,X,Z,output,fx] = dream_zs(MCMCPar,ModelName,[],ParRange,Measurement);

end;

if example == 9,    % Rainfall-runoff model with Whittle's log-likelihood (spectral analysis)

    % ---------------------------- Check the following 2 papers ------------------------------- %
    %                                                                                           %
    % G. Schoups, J.A. Vrugt, F. Fenicia, and N.C. van de Giesen (2010), Corruption of          %
    %     accuracy and efficiency of Markov Chain Monte Carlo simulation by inaccurate          %
    %     numerical implementation of conceptual hydrologic models, Water Resources             %
    %     Research, 46, W10530, doi:10.1029/2009WR008648.                                       %
    %                                                                                           %
    % G. Schoups, and J.A. Vrugt (2010), A formal likelihood function for parameter and         %
    %     predictive inference of hydrologic models with correlated, heteroscedastic and        %
    %     non-Gaussian errors, Water Resources Research, 46, W10531, doi:10.1029/2009WR008933.  %
    %                                                                                           %
    % ----------------------------------------------------------------------------------------- %

    % and for Whittle's likelihood function and application:

    % ----------------------------------------------------------------------------------------- %
    %                                                                                           %
    % A. Montanari, and E. Toth (2007), Calibration of hydrological models in the spectral      %
    % domain: An opportunity for scarcely gauged basins?, Water Resources Research, 43, W05434, %
    % doi:10.1029/2006WR005184.                                                                 %
    %                                                                                           %
    % ----------------------------------------------------------------------------------------- %

    % Problem specific parameter settings
    MCMCPar.n = 7;                          % Dimension of the problem
    MCMCPar.ndraw = 10000;                  % Maximum number of function evaluations
    MCMCPar.T = 1;                          % Each Tth sample is collected in the chains
    MCMCPar.prior = 'LHS';                  % Latin Hypercube sampling (options, "LHS", "COV" and "PRIOR")
    MCMCPar.BoundHandling = 'Reflect';      % Boundary handling (options, "Reflect", "Bound", "Fold", and "None");
    MCMCPar.lik = 9;                        % Whittle's likelood function -- spectral analysis of data

    % Define modelName
    ModelName = 'hmodel';

    % Give the parameter ranges (minimum and maximum values)
    %parno:       1     2     3     4     5     6     7      8    9     10    11   12   13   14   15   16   17   18   19   20  21
    %parname:     fA    Imax  Smax  Qsmax alE   alS   alF    Pf   Kfast Kslow std0 std1 beta xi   mu1  phi1 phi2 phi3 phi4 K   lambda
    Extra.fpar = [1     0     100   10    100   1e-6  1e-6   0    2     70    0.1  0    0    1    0    0    0    0    0    0   1];
    parmin =     [1     0     10    0     1e-6  1e-6 -10     0    0     0     0    0   -1    0.1  0    0    0    0    0    0   0.1 ];
    parmax =     [1     10    1000  100   100   1e-6  10     0    10    150   1    1    1    10   100  1    1    1    1    1   1];
    Extra.idx_vpar = [2 3 4 5 7 9 10];
    % Define parameter ranges
    ParRange.minn = parmin(Extra.idx_vpar); ParRange.maxn = parmax(Extra.idx_vpar);

    % Load the French Broad data
    daily_data = load('03451500.dly');

    % First two years are warm-up
    Extra.idx = [731:size(daily_data,1)]';

    % Define the PET, Measured Streamflow and Precipitation.
    Extra.Precip = daily_data(:,4); Extra.Ep = daily_data(:,5);

    % Define the measured streamflow data
    Measurement.MeasData = daily_data(Extra.idx,6);

    % Run the DREAM_ZS algorithm
    [Sequences,X,Z,output,fx] = dream_zs(MCMCPar,ModelName,Extra,ParRange,Measurement);

end;

if example == 10,    % the use of prior information in a multimodel mixture distrbibution

    % ---------------------------- Check the following two papers ----------------------------- %
    %                                                                                           %
    %   Vrugt, J.A., C.J.F. ter Braak, C.G.H. Diks, D. Higdon, B.A. Robinson, and J.M. Hyman,   %
    %       Accelerating Markov chain Monte Carlo simulation by differential evolution with     %
    %       self-adaptive randomized subspace sampling, International Journal of Nonlinear      %
    %       Sciences and Numerical Simulation, 10(3), 271-288, 2009.                            %
    %                                                                                           %
    %   Ter Braak, C.J.F., and J.A. Vrugt, Differential Evolution Markov Chain with snooker     %
    %       updater and fewer chains, Statistics and Computing, 10.1007/s11222-008-9104-9,      %
    %       2008.                                                                               %
    %                                                                                           %
    % ----------------------------------------------------------------------------------------- %

    % Problem specific parameter settings
    MCMCPar.n = 2;                          % Dimension of the problem (number of parameters to be estimated)
    MCMCPar.ndraw = 10000 ;                 % Maximum number of function evaluations
    MCMCPar.T = 1;                          % Each Tth sample is collected in the chains
    MCMCPar.prior = 'PRIOR';                % Latin Hypercube sampling (options, "LHS", "COV" and "PRIOR")
    MCMCPar.BoundHandling = 'None';         % Boundary handling (options, "Reflect", "Bound", "Fold", and "None");
    MCMCPar.lik = 1;                        % Model returns directly the probability density

    % Define modelName
    ModelName = 'mixturemodel';

    % Provide information to do covariance sampling ("COV")
    MCMCPar.mu = zeros(1,MCMCPar.n);        % Provide mean of initial sample
    MCMCPar.cov = eye(MCMCPar.n);           % Initial covariance

    % Specify the prior distributions for both parameters
    MCMCPar.prior_marginal = {  'normrnd(-5,0.1)',...
                                'normrnd(-5,0.1)',...
                                };
    % So the mixture models has two modes at -5 and 5; with the specified prior
    % distribution the mode around 5 should disappear. You can compare the
    % theoretical distribution with the DREAM(ZS) results by plotting the
    % target distribution and adding in the density derived from the posterior samples. 
    
    % Run the DREAM_ZS algorithm   
    [Sequences,X,Z,output,fx] = dream_zs(MCMCPar,ModelName);

end;

if example == 11,    % multivariate student t distribution with 60 degrees of freedom

    % ---------------------------- Check the following paper ---------------------------------- %
    %                                                                                           %
    %   Ter Braak, C.J.F., and J.A. Vrugt, Differential Evolution Markov Chain with snooker     %
    %       updater and fewer chains, Statistics and Computing, 10.1007/s11222-008-9104-9,      %
    %       2008.                                                                               %
    %                                                                                           %
    % ----------------------------------------------------------------------------------------- %

    % Problem specific parameter settings
    MCMCPar.n = 25;                         % Dimension of the problem (number of parameters to be estimated)
    MCMCPar.ndraw = 100000;                 % Maximum number of function evaluations
    MCMCPar.T = 10;                         % Each Tth sample is collected in the chains
    MCMCPar.prior = 'LHS';                  % Latin Hypercube sampling (options, "LHS", "COV" and "PRIOR")
    MCMCPar.BoundHandling = 'None';         % Boundary handling (options, "Reflect", "Bound", "Fold", and "None");
    MCMCPar.lik = 1;                        % Model returns density
    
    % Define modelName
    ModelName = 'multi_student';
    
    % Give the parameter ranges (minimum and maximum values)
    ParRange.minn = [-5 * ones(1,MCMCPar.n)]; ParRange.maxn = [15 * ones(1,MCMCPar.n)];

    % ---------------------- Define covariance matrix ---------------------
    
    % Construct the dxd correlation matrix
    Extra.corr = 0.5 * eye(MCMCPar.n) + 0.5 * ones(MCMCPar.n);
    
    % Define the example specific properties used to compute output
    Extra.mu = zeros(1,MCMCPar.n);          
    
    % How many degrees of freedom of student distribution used as target function?
    Extra.df = 60;

    % Make sure C is a valid covariance matrix
    [Extra.R,err] = cholcov(Extra.corr,0);
    
    % Define logSqrtDetC
    Extra.logSqrtDetC = sum(log(diag(Extra.R)));
    
    % Define dimensionality
    Extra.d = MCMCPar.n;

    % ---------------------------------------------------------------------

    % Run the DREAM_ZS algorithm   
    [Sequences,X,Z,output,fx] = dream_zs(MCMCPar,ModelName,Extra,ParRange);
    
end;

if example == 12,    % pedometrics problem involving variogram fitting

    % ---------------------------- Check the following paper ---------------------------------- %
    %                                                                                           %
    %   Minasny, B., J.A. Vrugt, and A.B. McBratney, Confronting uncertainty in model-based     %
    %       geostatistics using Markov chain Monte Carlo simulation, Geoderma, 163, 150-622,    %
    %       doi:10.1016/j.geoderma.2011.03.011.                                                 %
    %                                                                                           %   
    % ----------------------------------------------------------------------------------------- %

    % Problem specific parameter settings
    MCMCPar.n = 5;                          % Dimension of the problem (number of parameters to be estimated)
    MCMCPar.ndraw = 10000;                  % Maximum number of function evaluations
    MCMCPar.T = 1;                          % Each Tth sample is collected in the chains
    MCMCPar.prior = 'LHS';                  % Latin Hypercube sampling (options, "LHS", "COV" and "PRIOR")
    MCMCPar.BoundHandling = 'Reflect';      % Boundary handling (options, "Reflect", "Bound", "Fold", and "None");
    MCMCPar.lik = 4;                        % Define likelihood function -- log-likelihood function

    % Define modelName
    ModelName = 'blpmodel';

    % Give the parameter ranges (minimum and maximum values)
    ParRange.minn = [0 0.00 0.00 0.00 0]; ParRange.maxn = [100 100 1000 1000 20];

    % Load the data
    load 'forest.txt'; data = forest; x = data(:,1); y = data(:,2); z = data(:,3); 
    
    % ---------------------------------------------------------------------
    % Create matrix of coordinates
    X = [x , y];  
    % Create trend vector
    Extra.M = ones(length(z),1);
    % distance matrix between all observations
    Extra.H = distmat(X);   
    % Define np
    Extra.np = 1;
    % Define nt
    Extra.nt = 5;
    % Define Extra.N
    Extra.N = length(z);
    % Define z
    Extra.z = z;
    % ---------------------------------------------------------------------
    
    % Run the DREAM_ZS algorithm
    [Sequences,X,Z,output,fx] = dream_zs(MCMCPar,ModelName,Extra,ParRange);
    
    % Create a single matrix with values sampled by chains
    ParSet = GenParSet(Sequences,MCMCPar);

    % Postprocess the results to generate some fitting results
    postproc_variogram
    
end;

if example == 13,    % Nash-Cascade series of reservoirs

    % ---------------------------- Check the following 2 papers ------------------------------- %
    %                                                                                           %
    %   Nash, J.E., A unit hydrograph study with particular reference to British catchments,    %
    %      Proceedings - Institution of Civil Engineers, 17, 249-282, 1960.                     % 
    %                                                                                           %
    %   Nash, J.E., J.V. Sutcliffe, River flow forecasting through conceptual models part I -   %
    %      A discussion of principles, Journal of Hydrology, 10(3), 282-290, 1970.              %
    %                                                                                           %   
    % ----------------------------------------------------------------------------------------- %

    % Problem specific parameter settings
    MCMCPar.n = 1;                          % Dimension of the problem
    MCMCPar.ndraw = 10000;                  % Maximum number of function evaluations
    MCMCPar.T = 1;                          % Each Tth sample is collected in the chains
    MCMCPar.prior = 'LHS';                  % Latin Hypercube sampling (options, "LHS", "COV" and "PRIOR")
    MCMCPar.BoundHandling = 'Reflect';      % Boundary handling (options, "Reflect", "Bound", "Fold", and "None");
    MCMCPar.lik = 2;                        % Define likelihood function -- Sum of Squared Error
    MCMCPar.modout = 'Yes';                 % Return model (function) simulations of samples Yes or No)?

    % Define modelName
    ModelName = 'Nash_Cascade';

    % Give the parameter ranges (minimum and maximum values)
    ParRange.minn = [  1  ]; 
    ParRange.maxn = [ 100 ];

    % Load the French Broad data
    daily_data = load('03451500.dly');

    % Define maximum time. 
    Extra.maxT = 365;

    % Define the precipitation
    Extra.Precip = daily_data(1:Extra.maxT,4); 

    % Area factor to translate Nash-Cascade output from mm/d to m3/s 
    Extra.F = 767 * (1000 * 1000 ) / (1000 * 60 * 60 * 24);

    % Create the synthetic time series
    [S] = Nash_Cascade(2,Extra);

    % Now define the heteroscedastic measurement error
    Measurement.Sigma = max(0.2 * S , 1e-2);
    
    % And perturb the model simulation with a heteroscedastic measurement error
    Measurement.MeasData = normrnd(S,Measurement.Sigma);
    
    % Run the DREAM_ZS algorithm
    [Sequences,X,Z,output,fx] = dream_zs(MCMCPar,ModelName,Extra,ParRange,Measurement);

end;

if example == 14, % Rainfall runoff modeling problem using Schoups and Vrugt model (2010)
 
    % ---------------------------- Check the following 3 papers ------------------------------- %
    %                                                                                           %
    % J.A. Vrugt, and M. Sadegh, Towards diagnostic model calibration and evaluation:           %
    %     Appproximate Bayesian computation, Water Resources Research, 2012, In Review.         %
    %                                                                                           %
    % M. Sadegh, and J.A. Vrugt (2013), Approximate Bayesian computation using DREAM: Theory,   %
    %     Numerical Implementation and Case Studies, Water Resources Research, In Prep.         %
    %                                                                                           %
    % G. Schoups, and J.A. Vrugt (2010), A formal likelihood function for parameter and         %
    %     predictive inference of hydrologic models with correlated, heteroscedastic and        %
    %     non-Gaussian errors, Water Resources Research, 46, W10531, doi:10.1029/2009WR008933.  %
    %                                                                                           %
    % ----------------------------------------------------------------------------------------- %

    % Problem specific parameter settings
    MCMCPar.n = 7;                          % Dimension of the problem
    MCMCPar.ndraw = 25000;                  % Maximum number of function evaluations
    MCMCPar.T = 1;                          % Each Tth sample is collected in the chains
    MCMCPar.prior = 'LHS';                  % Latin Hypercube sampling (options, "LHS", "COV" and "PRIOR")
    MCMCPar.BoundHandling = 'Reflect';      % Boundary handling (options, "Reflect", "Bound", "Fold", and "None");
    MCMCPar.lik = 10;                       % ABC likelihood function (standardized)
    MCMCPar.ABC = 'Yes';                    % Specify that we perform ABC
    MCMCPar.delta = 0.025;                  % Delta of the noisy ABC implementation
    MCMCPar.rho = inline('X - Y');          % Define the distance - in this case the difference
    
    % Give the parameter ranges (minimum and maximum values)
    %parno         1     2      3      4     5      6      7    
    %parname:     Imax  Smax  Qsmax   alE   alF   Kfast  Kslow  

    % Minimum parameter values
    parmin =     [ 0     10     0    1e-6   -10     0      0    ];
    % maximum parameter values
    parmax =     [ 10   1000   100   100     10     10    150   ];
    
    % Select the parameters to be sampled
    Extra.idx_vpar = [1:7]; 
    
    % And define the associated Parameter ranges
    ParRange.minn = parmin(Extra.idx_vpar); ParRange.maxn = parmax(Extra.idx_vpar);

    % Load the French Broad data
    daily_data = load('03451500.dly');

    % First two years are warm-up
    Extra.idx = [731:size(daily_data,1)]';

    % Define the PET, Measured Streamflow and Precipitation.
    Extra.P = daily_data(:,4); Extra.Ep = daily_data(:,5);

    % Define the measured streamflow data
    MeasData = daily_data(Extra.idx,6);

    % Now calculate the summary metrics
    [Measurement.MeasData] = CalcMetrics(MeasData,Extra.P(Extra.idx))';
     
    % Define modelname
    ModelName = 'rainfall_runoff';

    % Run the DREAM_ZS algorithm
    [Sequences,X,Z,output,fx] = dream_zs(MCMCPar,ModelName,Extra,ParRange,Measurement);
    
end;

if example == 15,    % 10 bivariate normal distributions -- Approximate Bayesian Computation

    % ---------------------------- Check the following 2 papers ------------------------------- %
    %                                                                                           %
    % M. Sadegh, and J.A. Vrugt (2013), Approximate Bayesian computation using DREAM: Theory,   %
    %     Numerical Implementation and Case Studies, Water Resources Research, In Prep.         %
    %                                                                                           %
    % B.M. Turner, and P.B. Sederberg (2013), Approximate Bayesian computation with             %
    %     differential evolution, Journal of Mathematical Psychology, In Press.                 %
    %                                                                                           %
    % ----------------------------------------------------------------------------------------- %

    %How many bivariate normal distributions?
    Extra.Npairs = 10;
    
    % Problem specific parameter settings
    MCMCPar.n = 2 * Extra.Npairs;           % Dimension of the problem
    MCMCPar.ndraw = 30000;                  % Maximum number of function evaluations
    MCMCPar.T = 1;                          % Each Tth sample is collected in the chains
    MCMCPar.prior = 'LHS';                  % Latin Hypercube sampling (options, "LHS", "COV" and "PRIOR")
    MCMCPar.BoundHandling = 'Fold';         % Boundary handling (options, "Reflect", "Bound", "Fold", and "None");
    MCMCPar.lik = 10;                       % ABC likelihood function (standardized)
    MCMCPar.ABC = 'Yes';                    % Specify that we perform ABC
    MCMCPar.delta = 0.025;                  % Delta of the noisy ABC implementation
    
    % Define distance function -- in this case RMSE
    MCMCPar.rho = inline(' sqrt( 1 / 20 * sum((X - Y).^2)) ');  
    
    % Define modelName
    ModelName = 'ABC_binormal';

    % Give the parameter ranges (minimum and maximum values)
    ParRange.minn = [zeros(1,2*Extra.Npairs) ]; ParRange.maxn = [10 * ones(1,2*Extra.Npairs) ];

    % Lets create the data -- first create mu
    x = LHS(ParRange.minn,ParRange.maxn,1000); idx = round(unifrnd(0.5,1000.49)); x = x(idx,1:MCMCPar.n-1);
    
    % Now create data
    Measurement.MeasData = ABC_binormal(x,Extra);
        
    % Run the DREAM_ZS algorithm
    [Sequences,X,Z,output,fx] = dream_zs(MCMCPar,ModelName,Extra,ParRange,Measurement);

end;

if example == 16,   % Crosshole GPR slowness distribution based on a straight-ray approximation using the discrete cosine transform. 
                    % The problem is simplified compared to the paper cited below as it considers a problem in which the true model 
                    % represent a model with the same dimension as the inverse parameterization and it uses straight rays.

    % ### Results can be visualized by the function visualize_results.m ###

    % ---------------------------- Check the following paper ---------------------------------- %
    %                                                                                           %
    %  N. Linde, and J.A. Vrugt, Spatially distributed soil moisture from traveltime            %
    %       observations of crosshole ground penetrating radar using Markov chain Monte Carlo,  %
    %       Vadose Zone Journal, 2013                                                           %
    %                                                                                           %
    % ----------------------------------------------------------------------------------------- %
       
    Extra.parx = 8;                         % Inversion parameters in x-direction (DCT order)
    Extra.parz = 8;                         % Inversion parameters in z-direction (DCT order)
                                            % 64 dimensions in total
    
    % Problem specific parameter settings
    MCMCPar.n = Extra.parx * Extra.parz;    % Dimension of the problem (number of parameters to be estimated)
    MCMCPar.ndraw = 300000;                 % Maximum number of function evaluations
    MCMCPar.T = 25;                         % Each Tth sample is collected in the chains
    MCMCPar.prior = 'COV';                  % Latin Hypercube sampling (options, "LHS", "COV" and "PRIOR")
    MCMCPar.BoundHandling = 'None';         % Boundary handling (options, "Reflect", "Bound", "Fold", and "None");
    MCMCPar.lik = 4;                        % The model directly returns the log-density
    %MCMCPar.save='Yes';

    % Define distance function -- in this case RMSE
    MCMCPar.rho = inline(' sqrt( 1 / 20 * sum((X - Y).^2)) ');  
   
    % Define modelName
    ModelName = 'DCT_GPR';
    Extra.error = 0.5;                      % Standard deviation of Gaussian traveltime error

    % Give the parameter ranges (minimum and maximum values)
    ParRange.minn(1) = 30 * 0.7696; ParRange.maxn(1) = 30 * 1.301;              % Corresponds to 1/0.05 to 1/0.17 ns/m in logarithmic units

    % Set-up forward kernel
    % The x-axis is varying the fastest
    x = 0 : 0.1 : 3;                        % Boundaries of uniform x-grid (3 by 3 m grid); 0.1 m discretization
    z = 0 : 0.1:3;                          % Boundaries of uniform z-grid

    sourcex = 0.01;                         % x-position of source
    sourcez = 0.05 : 0.1 : 3.;              % z-positions of sources
    receiverx = 2.99;                       % x-position of receivers
    receiverz = 0.05 : .1 : 3;              % z-positions of receivers
    nsource = length(sourcez); nreceiver = length(receiverz);

    % Calculate acquisition geometry (multiple-offset gather)
    for j = 1 : nsource
        for i = 1 : nreceiver
            data( ( j - 1 ) * nreceiver + i , :) = [sourcex sourcez(j) receiverx receiverz(i)];
        end
    end

    % Calculate forward modeling kernel, Courtesy Dr. James Irving, UNIL
    Extra.J = tomokernel_straight(data,x,z); % Distance of ray-segment in each cell for each ray

    % Grid-cells in horizontal and vertical direction
    Extra.dimhor = length(x) - 1; Extra.dimver = length(z) - 1;

    Extra.por = 0.36;                       % Porosity field
    nz = 30;                                % Original discretization of water saturation model
    nx = 40;
    wcon = load('Sw.dat');                  % Load water saturation model: Courtesy of Dr. Michael Kowalsky (LBNL), Kowalsky et al. (2005; WRR)
    count = 0;
    for k = 1 : nz
        for i = 6 : nx - 5                  % Make model 3 by 3 m
            wcont( k , i - 5) = wcon( nz - k + 1 , i);
        end
    end
    wcont = wcont * Extra.por;              % Transform saturation into water content
    wcont = dct2(wcont);                    % Transform water content model into DCT space

    wtrunc = zeros(Extra.dimver,Extra.dimhor);
    for j = 1 : Extra.parz
        for i = 1 : Extra.parx
            wtrunc(j,i) = wcont(j,i);       % Truncate at same order as inverse parameterization
        end
    end
    wtrunc = idct2(wtrunc);                 % Do inverse DCT
    wtrunc = wtrunc';

    % Translate water content model into slowness using the Refractive Index/CRIM model
    Extra.pw = 81;                          % Permittivity of water
    Extra.pa = 1;                           % Permittivity of air
    Extra.ps = 5;                           % Permittivity of mineral grains
    slowtrue = wtrunc(:) * sqrt(Extra.pw) + (Extra.por-wtrunc(:)) * sqrt(Extra.pa) + (1 - Extra.por) * sqrt(Extra.ps); % sqrt of effective permittivity
    slowtrue = slowtrue / 0.3; % True slowness field

    % Simulate data for true slowness model
    % Add Gaussian uncorrelated noise with a standard deviation of 1 ns.
    Extra.datasim = Extra.J * slowtrue + Extra.error * randn(nsource * nreceiver,1); % Simulate data
    Extra.model_DCT = zeros(Extra.dimver,Extra.dimhor); % Matrix in which proposed model is assigned

    % Scale DCT coefficients such that all models are possible
    dum = zeros(Extra.dimver,Extra.dimhor);
    count = 1;
    for i = 1 : Extra.parz
        for j = 1 : Extra.parx
            if (i>1 | j>1)
                count = count + 1;
                dum(i,j) = 1;
                dummy = idct2(dum);
                ParRange.minn(count) = -1.7 / max(max(abs(dummy)));  % 0.2657
                ParRange.maxn(count) =  1.7 / max(max(abs(dummy)));  % 0.2657
                dum(i,j) = 0;
            end
        end
    end
    % ---------------------------------------------------------------------

    % Provide information to do initial sampling ("COV") --> The initial
    % chain positions are concentrated to the middle of the parameter ranges
    % This will speed up the calculation -- but cannot be done in practice!

    % Define mu and covariance of multinormal distribution used with "COV"
    MCMCPar.mu = ParRange.minn + 0.5 * ( ParRange.maxn - ParRange.minn );   % Provide mean of initial sample
    MCMCPar.cov = 0.001 * diag( ParRange.maxn - ParRange.minn );            % Initial covariance

    % Run the DREAM_ZS algorithm
    [Sequences,X,Z,output,fx] = dream_zs(MCMCPar,ModelName,Extra,ParRange);

end;

% Create a single matrix with values sampled by chains
ParSet = GenParSet(Sequences,MCMCPar);

%Save output for postprocessing
save GWF1P81R05.mat

% --------------------------------------------------------------------------------------------- %
% ------------------------------------ POSTPROCESSING ----------------------------------------- %
%                                                                                               %
% For postprocessing of results --> please go to directory \PostProcessing and run the          %  
% "postprocMCMC" script. This will compute various statistic and create a number of different   % 
% plots, including the R_stat convergence diagnostic, marginal posterior parameter              % 
% distributions, two-dimensional correlation plots of the posterior parameter samples, and      %  
% parameter and total uncertainty posterior simulation uncertainty ranges.                      %
%                                                                                               %
% --------------------------------------------------------------------------------------------- %
% --------------------------------------------------------------------------------------------- %