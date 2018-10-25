function [S_mod] = CalcMetrics(Y,P);
% Calculate the summary metrics of the discharge data

% Determine the size of the data
N = size(Y,1);

% Derive FDC
FDC = CalcFDC(Y);

% The second column is the exceedance probability
p = FDC(:,2);

% Now define the streamflow
x = FDC(:,1);

% Now fit a retention function
pars = [1 1];

% Fit the FDC
pars = fminsearch(@(pars) VG_Ret(pars,x,p),pars);

% Define the FDC metrics
S_mod(3) = pars(1); S_mod(4) = pars(2);

% Now derive the annual runoff coefficient
S_mod(1) = sum(Y)/sum(P);

% Now derive the annual baseflow coefficient
yb(1) = 0.25 * Y(1); phi = 0.925;

% Now loop
for j = 2:N;
    yb(j,1) = min( phi * yb(j-1,1) + (1/2) * (1 - phi) * ( Y(j,1) + Y(j-1,1) ) , Y(j,1));
end;

% Now calculate the annual baseflow coefficient
S_mod(2) = sum(yb)/sum(Y);