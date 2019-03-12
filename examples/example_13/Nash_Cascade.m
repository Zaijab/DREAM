function [SimRR] = Nash_Cascade(k,Extra);
% Nash-Cascade unit hydrograph -- three linear reservoirs in series

if k < 1,
    % Write to screen
    disp('Recession constant < 1 day --> numerical errors possible')
end;

% Define help matrix
A = zeros(Extra.maxT,Extra.maxT);

% Define number of linear reservoirs
n = 3;

% Define Time
Time = [ 1 : Extra.maxT ];

% Calculate unit hydrograph
IUH = 1 / ( k * Gamma(n) )  * ( Time / k ).^(n-1) .* exp( -Time / k ); 
    
% Now loop over time
for t = 1 : Extra.maxT,
    
    % Define idx
    idx = [ 1 : Extra.maxT - (t - 1) ];
    
    % Calculate flow
    A(t,t:Extra.maxT) = Extra.Precip(t) * IUH(idx);
    
end

% Now determine total flow
SimRR = Extra.F * sum(A)';