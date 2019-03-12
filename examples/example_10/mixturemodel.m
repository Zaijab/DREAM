function [p] = mixturemodel(x,Extra);
% Example of mixture model with two multinormal distributions centered at -5 and 5. 

% Calculate probability density at given value of x
p = 1/3 * mvnpdf(x, -5 , eye(size(x,2)) ) + 2/3 * mvnpdf(x, 5 , eye(size(x,2)) );