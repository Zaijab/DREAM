function [log_p] = normalfunc(x,Extra)
% Multivariate normal distribution with specified covariance matrix

% Calculate the log density
log_p = -0.5 * x * Extra.invC * x';