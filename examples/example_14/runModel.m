%% Example script to compile and run conceptual rainfall-runoff model
clear all;

%% Generate mex file (need to do this only once)
mex crr_model.c;

%% Mopex data
mopexdata = load('03451500.dly');
idx = find(mopexdata(:,1)>1959 & mopexdata(:,1)<1999); 
n = 365;
tout = 0:n;

%% Forcing and parameters
data.P     = mopexdata(idx(1:n),4)';    % daily rainfall (mm/d)
data.Ep    = mopexdata(idx(1:n),5)';    % daily evaporation (mm/d)
data.Imax  = 3;        % interception storage capacity (mm)
data.Sumax = 1000;      % unsaturated zone storage capacity (mm)
data.Qsmax = 10;        % maximum percolation rate (mm/d)
data.aE    = 50;        % evaporation coefficient
data.aF    = -1;        % runoff coefficient
data.aS    = 1e-6;      % percolation coefficient
data.Kf    = 3;         % fast-flow response time (d)
data.Ks    = 70;        % slow-flow response time (d)

%% Integration options
options.InitialStep = 1;                 % initial time-step (d)
options.MaxStep     = 1;                 % maximum time-step (d)
options.MinStep     = 1e-6;              % minimum time-step (d)
options.RelTol      = 1e-3;              % relative tolerance
options.AbsTol      = 1e-3*ones(5,1);    % absolute tolerances (mm)
options.Order       = 2;                 % 2nd order accurate method (Heun)

%% Initial conditions
y0 = 1e-6*ones(5,1);

%% Run model C
tic;
y = crr_model(tout,y0,data,options);
cpu_C = toc

%figure;  plot(tout(2:end),y(1,2:end),'k');      % interception storage
%hold on; plot(tout(2:end),y(2,2:end),'b');      % unsat zone storage
%hold on; plot(tout(2:end),y(3,2:end),'r');      % fast-flow zone storage
%hold on; plot(tout(2:end),y(4,2:end),'g');      % slow-flow zone storage
figure; plot(tout(2:end),y(5,2:end)-y(5,1:end-1),'k'); % discharge
hold on;
%% Run model Matlab
tic;
ns = n; u0 = y0'; y = [];
for s = 1:ns
    [dummy,u] = ode23(@(t,u) mtlb_fRhs(t,u,s,data),[0 1],u0,options);
    u0 = u(end,:);
    y = [y; u0];
end
y = [y0 y'];
cpu_matlab = toc

%figure;  plot(tout(2:end),y(1,2:end),'k');      % interception storage
%hold on; plot(tout(2:end),y(2,2:end),'b');      % unsat zone storage
%hold on; plot(tout(2:end),y(3,2:end),'r');      % fast-flow zone storage
%hold on; plot(tout(2:end),y(4,2:end),'g');      % slow-flow zone storage
%figure; 
plot(tout(2:end),y(5,2:end)-y(5,1:end-1),'r.'); % discharge
