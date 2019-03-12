function [S_mod,Y] = rainfall_runoff(theta,data);
% Rainfall runoff model 

%% Assign parameters
data.Imax  = theta(1);      % interception storage capacity (mm)
data.Sumax = theta(2);      % unsaturated zone storage capacity (mm)
data.Qsmax = theta(3);      % maximum percolation rate (mm/d)
data.aE    = theta(4);      % evaporation coefficient
data.aF    = theta(5);      % runoff coefficient
data.aS    = 1e-6;          % percolation coefficient
data.Kf    = theta(6);      % fast-flow response time (d)
data.Ks    = theta(7);      % slow-flow response time (d)

%% Integration options
options.InitialStep = 1;                 % initial time-step (d)
options.MaxStep     = 1;                 % maximum time-step (d)
options.MinStep     = 1e-5;              % minimum time-step (d)
options.RelTol      = 1e-3;              % relative tolerance
options.AbsTol      = 1e-3*ones(5,1);    % absolute tolerances (mm)
options.Order       = 2;                 % 2nd order accurate method (Heun)

%% Initial conditions
y0 = 1e-5 * ones(5,1); 

%% Running time
tout = [ 0 : size(data.P,1) ];

%% Run model C
y = crr_model(tout,y0,data,options);

% Now compute discharge
Y = diff(y(5,1:end))'; Y = Y(data.idx);

% Now calculate the summary metrics
[S_mod] = CalcMetrics(Y,data.P(data.idx))';