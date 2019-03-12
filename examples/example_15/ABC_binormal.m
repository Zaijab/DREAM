function [S] = ABC_binormal(x,Extra);
% Generate samples from bivariate normal distribution

% Check whether delta has been specified
if size(x,2) == 2 * Extra.Npairs,
    % Set delta to zero
    x = [x 0];
end;

% Define the sigma of the bivariate distribution
sigma = 0.01^2 * eye(2); 

% Now reorganize the x-vector so that we get bivariate mu values
idx_1 = [1 : 2 : 2 * Extra.Npairs - 1]; idx_2 = [2 : 2 : 2 * Extra.Npairs];

% Define mu
mu = [x(idx_1)' x(idx_2)'];

% Now loop over each bivariate normal distribution
for i = 1:Extra.Npairs,
    % Sample 50 points from each bivariate normal distribution
    X(i,1:2) = mean ( mvnrnd ( mu(i,1:2) , sigma, 50 ) + normrnd( 0 , x(end) , 50 , 2 ) );
end

% Now return vector with model simulated values
S = X(:);