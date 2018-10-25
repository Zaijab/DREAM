function [z0, zs]=predict_blp(x,y,z,x0,y0,c0,c1,range,nu,beta,M,m)
% Best Linear Unbiased Predictor
% z(x0) = m(x) beta + err(x)
% Input:
%   x,y,z   : data
%   x0,y0   : coord to be predicted
%   c0,c1,range,nu  : Matern paremeter
%   [optional]
%   beta    : parameter of trend model
%   M       : design matrix for x,y
%   m       : design matrix for x0,y0
% Output:
%   z0      : predicted
%   zs      : MSE/variance of predicted

n=length(x);
n0=length(x0);
if(nargin<=10),
    M=ones(n,1);                % M=trend matrix
    m=ones(n0,1);               % m(x0)
end

% between data points
X = [x,y];
H = distmat(X, X);    % distance
C = cova_matern(c0,c1,range,nu,H);    % covariance matrix
CI = inv(C);

%[L,D,E,pneg]=mchol(C);
%U = (L*D^(1/2))';
%UI=inv(U);
%CI = UI*UI';

if(nargin<10),
    W = M' * CI * M;
    beta = inv(W) * M' * CI * z;   % parameter of model
end

for i=1:n0
    xp = x0(i);
    yp = y0(i);
    mp = m(i,:)';
    [z0(i,:),zs(i,:)] = blp(xp,yp,x,y,z,c0,c1,range,nu,M,beta,mp,CI);
end
