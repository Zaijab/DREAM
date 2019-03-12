function [S] = Banshp(x,Extra);
% Hyperbolic shaped posterior probability distribution

% Each line is one parameter combination
x(:,2) = x(:,2) + Extra.bpar*x(:,1).^2 - 100*Extra.bpar;

% Compute the SSR
S = -0.5 * (x * Extra.imat) * x';	