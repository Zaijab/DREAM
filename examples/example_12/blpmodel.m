function [L] = blpmodel(p,Extra);
% p contains the parameters of the linear model and variogram function best linear predictor

% parameter of the linear model
beta = p(1:Extra.np);

% the linear model
S = Extra.M * beta';

% parameter of the variogram
theta = p ( Extra.np + 1 : Extra.nt );   

% compute covariance structure of the data
C = cova_matern(theta(1),theta(2),theta(3),theta(4),Extra.H);    

% Compute Log likelihood
dC = det(C) + realmin;

% Calculate the L1 norm
L1 = ( Extra.z - S )'* inv(C) *( Extra.z - S );

% Calculate the log-likelihood
L  = -0.5 * ( log(dC) + L1 );