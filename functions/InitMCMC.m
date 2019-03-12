function [MCMCPar,pCR,lCR,CR,Iter,iteration,T,fx,m_func,Sequences,Z,Table_JumpRate,iloc,output] = InitMCMC(MCMCPar,Measurement)
% Initializes important variables for use in the algorithm

% Set MCMCPar.m
MCMCPar.m = MCMCPar.m0; 

% Define the crossover values as geometrical series
MCMCPar.CR = cumsum((1/MCMCPar.nCR) * ones(1,MCMCPar.nCR));

% Derive the number of elements in the output file
Nelem = floor(MCMCPar.ndraw/MCMCPar.seq) + 1;

% Initialize output information -- N_CR  
output.CR = zeros(floor(Nelem/MCMCPar.steps),MCMCPar.nCR+1);

% Initialize output information -- AR
output.AR = zeros(floor(Nelem/MCMCPar.steps),2); output.AR(1,1:2) = [MCMCPar.seq -1];

% Initialize output information -- R statistic
output.R_stat = zeros(floor(Nelem/MCMCPar.steps),MCMCPar.n+1);

% Calculate multinomial probabilities of each of the nCR CR values
pCR = (1/MCMCPar.nCR) * ones(1,MCMCPar.nCR);

% Calculate the actual CR values based on p
[CR,lCR] = GenCR(MCMCPar,pCR); 

% Check whether we store all the model simulations or not
if strcmp(MCMCPar.modout,'Yes');
    if Measurement.N > 0
        % Define matrix with model simulations
        fx = zeros(Measurement.N,floor(MCMCPar.ndraw/MCMCPar.T));
        % Set m_func
        m_func = MCMCPar.seq;
    else
        fx = []; m_func = [];
    end;
else
    fx = []; m_func = [];
end;

% Initialize Sequences with zeros
Sequences = zeros(floor(Nelem/MCMCPar.T),MCMCPar.n+2,MCMCPar.seq);

% Define Z 
Z = zeros(floor(MCMCPar.m0 + MCMCPar.seq * (MCMCPar.ndraw - MCMCPar.m0) / (MCMCPar.seq * MCMCPar.k)),MCMCPar.n+2);  %something wrong (remove Par.seq)

% Generate the Table with JumpRates (dependent on number of dimensions and number of pairs
for zz = 1:MCMCPar.DEpairs,
    Table_JumpRate(:,zz) = 2.38./sqrt(2 * zz * [1:MCMCPar.n]'); 
end;

% Initialize Iter 
Iter = MCMCPar.seq; iloc = 1; iteration = 2; T = 0;

% Change MCMCPar.steps to make sure to get nice iteration numbers in first loop
MCMCPar.steps = MCMCPar.steps - 1;