function [SimRR] = hymodMATLAB(Pars,Extra);
% Runs the HYMOD model

% Define the rainfall
PET = Extra.E; Precip = Extra.P; MaxT = Extra.MaxT;
% Define the parameters
cmax = Pars(1); bexp = Pars(2); alpha = Pars(3); Rs = Pars(4); Rq = Pars(5);

% Set the initial states
x_loss = 0.0;

% Initialize slow tank state
x_slow = 0; % --> works ok if calibration data starts with low discharge

% Initialize state(s) of quick tank(s)
x_quick(1:3,1) = 0; outflow = [];

% Now loop over the forcing data
for t = 1:MaxT,
    
    % Assign precipitation and evapotranspiration
    Pval = Precip(t,1); PETval = PET(t,1);
    
    % Compute excess precipitation and evaporation
    [ER1,ER2,x_loss] = excess(x_loss,cmax,bexp,Pval,PETval);
    
    % Calculate total effective rainfall
    ET = ER1 + ER2;
    
    % Now partition ER between quick and slow flow reservoirs
    UQ = alpha*ET; US = (1-alpha)*ET;
  
    % Route slow flow component with single linear reservoir
    [x_slow,QS] = linres(x_slow,US,Rs);
    
    % Route quick flow component with linear reservoirs
    inflow = UQ; 
    
    for k = 1:3,
        % Linear reservoir
        [x_quick(k),outflow] = linres(x_quick(k),inflow,Rq); inflow = outflow;
    end;

    % Compute total flow for timestep
    output(t,1) = (QS + outflow);
    
end;

SimRR = Extra.F * output(65:MaxT,1);