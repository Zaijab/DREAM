function y = exppdf(x,mu);
% Exponential distribution
if x < 0,
    y = 0;
else
    % Overwrite of existing function
    y = mu * exp(-mu * x);
end;