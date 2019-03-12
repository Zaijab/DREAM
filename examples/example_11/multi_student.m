function [p] = multi_student(x,Extra);
% Multivariate Student t distribution with 60 degrees of freedom

% Extra.df defines the degrees of freedom
% Extra.d defines the dimensionality
% Extra.R defines the cholesky of the correlation matrix

% Normalize x with R
Z = x / Extra.R;

% Define logNumer
logNumer = -((Extra.df + Extra.d)/2) .* log(1+sum(Z.^2, 2)./Extra.df);

% Define logDenom
logDenom = Extra.logSqrtDetC + (Extra.d/2) * log(Extra.df*pi);

% Define multivariate pdf
p = exp(gammaln((Extra.df + Extra.d)/2) - gammaln(Extra.df/2) + logNumer - logDenom);