function [log_p] = Banana_func(x,Extra);
% Hyperbolic shaped posterior probability distribution

% Define the nonlinearity of the banana shaped distribution
b = 0.1; 

% Each line is one parameter combination
x(2) = x(2) + b * x(1).^2 - 100 * b;

% Calculate the log density
log_p = -0.5 * x * Extra.invC * x';