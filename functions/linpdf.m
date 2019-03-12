function [pdf] = linpdf(x,a,b);
% Function for updraft velocity

% Assumption is that pdf(x) is as follows
% "y(x) = max(a + b * x,0)"

% Now define the pdf
pdf = max(a + b * x,0);