clear all
close all
clc
% (mu/mu_w, lambda)-CMA-ES
% URL: http://www.lri.fr/~hansen/purecmaes.m
% References: See end of file. Last change: October, 21, 2010
% --------------------  Initialization --------------------------------
%Problem setup files
strfitnessfct = 'Case01';       %Opt 0.1,0.02,0.06
global TrueHead %#ok<NUSED>
global CMARange
global ParaRange
global rangemin
global rangemax
global t
global SMC
global MSMC
global MCMC
global MCMCr
global fSMC
global fMSMC
global fMCMC

load TrueHead
rangemin=-3;
rangemax=3;
CMARange=[rangemin;rangemax];
ParaRange=[1e-3;1e-1];
DisplayUnit=1;

%Optimization Setup
Learning =0;
N = 81;                          % number of objective variables/problem dimension
sigma = 0.5;                     % coordinate wise standard deviation (step size)
stopfitness = 0;              % stop if fitness < stopfitness (minimization)
lambda=4+floor(3*log(N));        % population size, offspring number 4+floor(3*log(N));
stopeval =3e6;                   % stop after stopeval number of function evaluations 1e1*N^4;
NumRun=1;                        % Number of Optimization runs for statistical analysis
OptRun=1;

% --------------------  Initialization --------------------------------
xmean = rand(N,1);                      % objective variables initial point

% Strategy parameter setting: Selection
mu = lambda/2;                          % number of parents/points for recombination
weights = log(mu+1/2)-log(1:mu)';       % muXone array for weighted recombination
mu = floor(mu);
weights = weights/sum(weights);         % normalize recombination weights array
mueff=sum(weights)^2/sum(weights.^2);   % variance-effectiveness of sum w_i x_i

% Strategy parameter setting: Adaptation
cc = (4 + mueff/N) / (N+4 + 2*mueff/N); % time constant for cumulation for C
cs = (mueff+2) / (N+mueff+5);           % t-const for cumulation for sigma control

if Learning ==0
    c1 = 2 / ((N+1.3)^2+mueff);             % learning rate for rank-one update of C
    cmu = 2 * (mueff-2+1/mueff) / ((N+2)^2+mueff);      % and for rank-mu update
else
    mueff=1;
    M=(N+2)/3;
    c1=min(1,lambda/6)/(M+(2*sqrt(M))+(mueff/N));
    LRate=(0.3+mueff-2+(1/mueff))/(M+(4*sqrt(M))+(mueff/2));
    cmu=min(c1-1,LRate);
end

damps = 1 + 2*max(0, sqrt((mueff-1)/(N+1))-1) + cs; % damping for sigma usually close to 1

% Initialize dynamic (internal) strategy parameters and constants
pc = zeros(N,1); ps = zeros(N,1);   % evolution paths for C and sigma
B = eye(N,N);                       % B defines the coordinate system
D = ones(N,1);                      % diagonal D defines the scaling
C = B * diag(D.^2) * B';            % covariance matrix C
invsqrtC = B * diag(D.^-1) * B';    % C^-1/2
eigeneval = 0;                      % track update of B and D
chiN=N^0.5*(1-1/(4*N)+1/(21*N^2));  % expectation of ||N(0,I)|| == norm(randn(N,1))
%  out.dat = []; out.datx = [];  % for plotting output
outdat=NaN(ceil(stopeval/lambda),N+2);
GCount=0;

%Block-wise sampler initialization for t=1
t=1;                                          %Sampling steps
xMC=NaN(1,N,lambda);                          %Estimands (non-scaled)
SMC=NaN(floor(stopeval/lambda),N,lambda);     %Seq MC
MCMCr=NaN(floor(stopeval/lambda),N,lambda);     %Seq MC
fSMC=NaN(floor(stopeval/lambda),lambda);      %f(x)
for k=1:lambda
    MCMCr(1,:,k)=xmean + sigma * B * (D .* randn(N,1));     %First set of estimands
    SMC(1,:,k) = interp1(CMARange,ParaRange,MCMCr(1,:,k));
    HC=reshape(SMC(1,:,k),9,9);
    fSMC(1,k)=MODFLOW(HC,k);                              %Objective function call
end
MSMC= SMC;      %Meteropolis seq MC
MCMC= SMC;      %MCMC
MCMCr= SMC;     %Meteropolis seq MC
fMSMC=fSMC;
fMCMC=fSMC;



% -------------------- Generation Loop --------------------------------
counteval = 0;  % the next 40 lines contain the 20 lines of interesting code
CurrentFitness=inf;
while counteval < stopeval && CurrentFitness >= stopfitness;
    % Generate and evaluate lambda offspring
    tic
    arx=zeros(N,lambda);
    arfitness=zeros(1,lambda);
    t=t+1;
    for k=1:lambda,
        e=sigma * B * (D .* randn(N,1));
        arx(:,k) = xmean + e;   % m + sig * Normal(0,C)
        arfitness(k) = feval(strfitnessfct, arx(:,k),k,e);      % Objective function call
        counteval = counteval+1;
    end
    
    % Sort by fitness and compute weighted mean into xmean
    [arfitness, arindex] = sort(arfitness);  % minimization
    xold = xmean;
    xmean = arx(:,arindex(1:mu)) * weights;  % recombination, new mean value
    xmin = arx(:, arindex(1)); % Return best point of last iteration.
    
    % Cumulation: Update evolution paths
    ps = (1-cs) * ps ...
        + sqrt(cs*(2-cs)*mueff) * invsqrtC * (xmean-xold) / sigma;
    hsig = sum(ps.^2)/(1-(1-cs)^(2*counteval/lambda))/N < 2 + 4/(N+1);
    pc = (1-cc) * pc ...
        + hsig * sqrt(cc*(2-cc)*mueff) * (xmean-xold) / sigma;
    
    % Adapt covariance matrix C
    artmp = (1/sigma) * (arx(:,arindex(1:mu)) - repmat(xold,1,mu));  % mu difference vectors
    C = (1-c1-cmu) * C ...                      % regard old matrix
        + c1 * (pc * pc' ...                    % plus rank one update
        + (1-hsig) * cc*(2-cc) * C) ...         % minor correction if hsig==0
        + cmu * artmp * diag(weights) * artmp'; % plus rank mu update
    
    % Adapt step size sigma
    sigma = sigma * exp((cs/damps)*(norm(ps)/chiN - 1));
    
    % Update B and D from C
    if counteval - eigeneval > lambda/(c1+cmu)/N/10 % to achieve O(N^2)
        eigeneval = counteval;
        C = triu(C) + triu(C,1)';                   % enforce symmetry
        if Learning==0
            [B,D] = eig(C);                         % eigen decomposition, B==normalized eigenvectors
        else
            D=diag(diag(C));
        end
        D = sqrt(diag(D));                          % D contains standard deviations now
        invsqrtC = B * diag(D.^-1) * B';
    end
    CurrentFitness=arfitness(1);
    Generation_Time=toc;
    
    % Displaying and saving results
    more off;  % turn pagination off in Octave
    GCount= GCount+1;
    if round(counteval/(lambda*DisplayUnit))== counteval/(lambda*DisplayUnit)
        disp([num2str(counteval/lambda) ' : ' num2str(arfitness(1)) ' : ' num2str(Generation_Time)])
    end
    outdat(GCount,:)=[arfitness(1)  Generation_Time xmin'];
end



