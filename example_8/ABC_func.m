function [X] = ABC_func(x,Extra);
% Simple one-dimensional example to illustrate ABC

% Example function
if rand < 1/2
    X = normrnd(x(1),1/10); 
else
    X = normrnd(x(1),1);
end;