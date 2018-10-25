function [y] = linrnd(a,b);
% Function for updraft velocity

% Assumption is that pdf(x) is as follows
% "y(x) = max(a + b * x,0)"

% Check the number of input arguments
if nargin < 3,
    n = 1; % default value
end;

% Solve for which x-value the function is zero
x_0 = -a/b; 

% Define the pdf
if x_0 < 0,
    x = linspace(x_0,0,1000);
else
    x = linspace(0,x_0,1000);
end;

% Now define the pdf
pdf = a + b * x;

% Now create the cdf
cdf = cumsum(pdf)./sum(pdf);

% Make sure the values of "cdf" are distinct
idx = find(diff(cdf) == 0) + 1;

% How many elements?
N = prod(size(idx));

% Then update these values with a very small value
cdf(idx) = [1:N] * 1e-10 + cdf(idx);

% Now create "n" random values of weights consistent with the prior cdf
w = rand(n,1);

if w <= x(2),
    
    % Just set the cdf value
    y = cdf(1);
    
else

    % And calculate the corresponding x values
    y = interp1(cdf,x,w);
    
end;